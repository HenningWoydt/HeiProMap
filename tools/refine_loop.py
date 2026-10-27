#!/usr/bin/env python3
"""Refinement loop around the wse-spmv scheduler.

M0 is a Waver-HeiProMap mapping (waver-heipromap-spectral, as in wse-spmv/run_all_mappers.py) or a given one. In every
round the current best mapping gets CHILDREN random variants: a few swaps of the contents of two nearby PEs and a few
single-row moves to a neighboring PE, drawn near the PEs of the critical chain of the best mapping's schedule. All
variants are scheduled with spmv-map-eval; the best of the parent and the children becomes the new best (so the
makespan never rises). Everything is stored in the run directory:

  config.json            the settings of the run
  log.jsonl              one line per evaluated mapping (round 0 = M0, then every child); makespan is the objective
                         over the scheduler seeds, makespans the value of each seed
  holdout.json           makespans of M0 and the final best mapping with scheduler seeds that were not optimized
  best_rNN.u32           the full best mapping after round NN (raw little-endian u32, row -> PE), only when it changed
  rNN/child_CC.npz       child CC of round NN as a difference to the best mapping of the round before (vertices, pe)
  rNN/chain.json         the critical chain of the best mapping that was the parent of round NN
  final_viz/             the schedule dump (--dump-viz) of the final best mapping

  python3 tools/refine_loop.py --rounds 10 --children 10 [--start partition.txt] [--out DIR]
  python3 tools/refine_loop.py --instance data/small100/instance.json ...   (smaller instance, see make_small_instance.py)
  python3 tools/refine_loop.py --reconstruct RUN_DIR ROUND CHILD OUT.u32     (rebuild a stored child)

A raw .u32 mapping can be read with numpy.fromfile(path, dtype='<u4'); entry r is the PE (row * 316 + col) of matrix
row r (poisson3d 215 x 215 x 216, row = i + 215 * (j + 215 * k)).
"""
import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
HEI = HERE.parent
WSE = Path(os.path.expanduser("~/Projects/wse-spmv"))

# The instance (set by set_instance): the default is the full one of wse-spmv/run_all_mappers.py.
FULL = {"nx": 215, "ny": 215, "nz": 216, "cols": 316, "rows": 316,
        "graph": "/tmp/heipromap_input.graph", "topology": "/tmp/heipromap_topology.graph",
        "base": "/tmp/heipromap_base.wspvmap"}
NX = NY = NZ = N = COLS = ROWS = P = CAP = 0


def set_instance(inst):
    """Sizes of the instance; CAP = ceil(1.03 * N / P) is the balance limit of Waver-HeiProMap."""
    global NX, NY, NZ, N, COLS, ROWS, P, CAP, INSTANCE
    INSTANCE = inst
    NX, NY, NZ, COLS, ROWS = inst["nx"], inst["ny"], inst["nz"], inst["cols"], inst["rows"]
    N = NX * NY * NZ
    P = COLS * ROWS
    CAP = -(-103 * N // (100 * P))


def eval_template(ram, inst):
    """The evaluator command of run_all_mappers.py with the sizes of the instance."""
    t = list(ram.COMMAND_TEMPLATE)
    for flag, key in (("--nx", "nx"), ("--ny", "ny"), ("--nz", "nz"), ("--cols", "cols"), ("--rows", "rows")):
        t[t.index(flag) + 1] = str(inst[key])
    return t


def load_wse():
    """The command template and helpers of run_all_mappers.py (its paths are relative to the wse-spmv root)."""
    os.chdir(WSE)
    sys.path.insert(0, str(WSE))
    argv, sys.argv = sys.argv, [sys.argv[0]]
    import run_all_mappers as ram
    sys.argv = argv
    return ram


def neighbors(u):
    i = u % NX
    j = (u // NX) % NY
    k = u // (NX * NY)
    r = []
    if i > 0: r.append(u - 1)
    if i < NX - 1: r.append(u + 1)
    if j > 0: r.append(u - NX)
    if j < NY - 1: r.append(u + NX)
    if k > 0: r.append(u - NX * NY)
    if k < NZ - 1: r.append(u + NX * NY)
    return r


def read_mapping(path):
    path = str(path)
    if path.endswith(".u32"):
        return np.fromfile(path, dtype="<u4")
    return np.fromfile(path, dtype=np.uint32, sep=" ")  # one PE id per line


# ------------------------------------------------------------------------------------------------------------------
# evaluation
# ------------------------------------------------------------------------------------------------------------------

class Evaluator:
    """Schedules mappings with spmv-map-eval. The objective of a mapping is a statistic of its makespans over several
    scheduler seeds: the makespan of a single seed varies by several percent for the same mapping."""

    def __init__(self, ram, out, tail_repair, tail_passes, template, base, objective="median"):
        self.ram = ram
        self.template = template
        self.base = base
        self.objective = objective
        self.work = out / "work"
        self.work.mkdir(exist_ok=True)
        self.pack = out / "wspvmap_pack"
        if not self.pack.exists():
            subprocess.run(["gcc", "-O2", "-o", str(self.pack), str(HERE / "wspvmap_pack.c")], check=True)
        self.flags = []
        if tail_repair != "off":
            self.flags = ["--tail-repair", tail_repair, "--tail-repair-passes", str(tail_passes)]

    def _pack(self, part, tag):
        raw = self.work / f"{tag}.u32"
        mapping = self.work / f"{tag}.wspvmap1"
        part.astype("<u4").tofile(raw)
        subprocess.run([str(self.pack), self.base, str(raw), str(mapping)], check=True)
        raw.unlink()
        return mapping

    def _schedule(self, mapping, seed, dump_dir=None):
        cmd = self.template + ["--map-seed", "0", "--mapping-in", str(mapping), "--schedule-seed", str(seed)]
        cmd += self.flags
        if dump_dir is not None:
            cmd += ["--dump-viz", str(dump_dir), "--dump-viz-detail-limit", "200000000"]
        r = subprocess.run(self.ram.limited(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            d = json.loads(r.stdout.strip().split("\n")[-1])
            sc, ph = d["schedule"], d["physical"]
            return {"makespan": sc["makespan"], "res_lb": sc["resource_lower_bound"],
                    "hop_volume": ph["hop_volume"], "max_link": ph["max_directed_link_load"],
                    "chain_message": sc["critical_chain"]["message"], "chain_length": sc["critical_chain"]["messages"]}
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(f"scheduling {mapping} failed (exit {r.returncode}): {e}\n{r.stderr[-500:]}")

    def run(self, part, tag, seeds, want_chain=False):
        """Metrics of the mapping: `makespan` is the objective over the seeds, `makespans` the value of every seed.
        With want_chain, `chain_pes` is the union of the PEs of the critical chains of the schedules of all seeds."""
        mapping = self._pack(part, tag)
        results, chain = [], set()
        for seed in seeds:
            dump = self.work / f"{tag}_dump{seed}" if want_chain else None
            r = self._schedule(mapping, seed, dump)
            if want_chain:
                chain.update(chain_pes(dump, r["chain_message"])[0])
                shutil.rmtree(dump)
            results.append(r)
        mapping.unlink()
        makespans = [r["makespan"] for r in results]
        stat = {"median": float(np.median(makespans)), "mean": float(np.mean(makespans)), "max": max(makespans)}
        out = dict(results[0])
        out.update(makespan=stat[self.objective], makespans=makespans, chain_pes=sorted(chain))
        return out

    def dump(self, part, tag, seed, dump_dir):
        """The schedule dump (--dump-viz) of the mapping for one seed."""
        mapping = self._pack(part, tag)
        self._schedule(mapping, seed, dump_dir)
        mapping.unlink()


def chain_pes(dump_dir, last_message):
    """PEs (row * 316 + col) and messages of the critical chain: the binding predecessors of the last message."""
    meta = json.load(open(Path(dump_dir) / "schedule-meta.json"))
    n = meta["num_messages"]
    buf = gzip.open(Path(dump_dir) / meta["files"]["messages"]).read()
    src = np.frombuffer(buf, "<u4", n, 0)
    dst = np.frombuffer(buf, "<u4", n, 4 * n)
    binding = np.frombuffer(buf, "<u4", n, 12 * n)
    chain, cur = [], int(last_message)
    while cur != 0xFFFFFFFF and cur not in chain:
        chain.append(cur)
        cur = int(binding[cur])
    chain.reverse()
    pes = sorted({int(src[c]) for c in chain} | {int(dst[c]) for c in chain})
    return pes, chain


# ------------------------------------------------------------------------------------------------------------------
# random variants
# ------------------------------------------------------------------------------------------------------------------

def region_of(pes, margin):
    """The PEs within `margin` (per axis) of a PE of the list."""
    out = set()
    for p in pes:
        x, y = p % COLS, p // COLS
        for dy in range(-margin, margin + 1):
            for dx in range(-margin, margin + 1):
                if 0 <= x + dx < COLS and 0 <= y + dy < ROWS:
                    out.add((y + dy) * COLS + x + dx)
    return np.array(sorted(out), dtype=np.int64)


def make_child(parent, order, starts, counts0, region, rng, n_swaps, n_moves, radius):
    child = parent.copy()
    counts = counts0.copy()
    used = set()
    swaps = tries = 0
    while swaps < n_swaps and tries < 200 * n_swaps:
        tries += 1
        a = int(rng.choice(region))
        bx = a % COLS + int(rng.integers(-radius, radius + 1))
        by = a // COLS + int(rng.integers(-radius, radius + 1))
        if not (0 <= bx < COLS and 0 <= by < ROWS):
            continue
        b = by * COLS + bx
        if b == a or a in used or b in used:
            continue
        child[order[starts[a]:starts[a + 1]]] = b  # the contents of two PEs trade places
        child[order[starts[b]:starts[b + 1]]] = a
        counts[a], counts[b] = counts[b], counts[a]
        used.update((a, b))
        swaps += 1
    moves = tries = 0
    while moves < n_moves and tries < 400 * n_moves:
        tries += 1
        a = int(rng.choice(region))
        rows = order[starts[a]:starts[a + 1]]
        u = int(rows[rng.integers(len(rows))])
        cur = int(child[u])
        others = [int(child[v]) for v in neighbors(u) if child[v] != cur]
        if not others:
            continue
        target = others[int(rng.integers(len(others)))]
        if counts[target] >= CAP or counts[cur] <= 1:
            continue
        child[u] = target
        counts[cur] -= 1
        counts[target] += 1
        moves += 1
    return child, swaps, moves


# ------------------------------------------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------------------------------------------

def reconstruct(run_dir, rnd, child, out_path):
    run_dir = Path(run_dir)
    entry = None
    for line in open(run_dir / "log.jsonl"):
        e = json.loads(line)
        if e["round"] == rnd and e["child"] == child:
            entry = e
    if entry is None:
        sys.exit("no such mapping in the log")
    part = read_mapping(run_dir / entry["parent"]).copy()
    delta = np.load(run_dir / f"r{rnd:02d}" / f"child_{child:02d}.npz")
    part[delta["vertices"]] = delta["pe"]
    part.astype("<u4").tofile(out_path)
    print(f"wrote {out_path} (round {rnd}, child {child}, makespan {entry['makespan']})")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="run directory (default data/refine_loop/run_<time> in the HeiProMap repo)")
    ap.add_argument("--start", help="starting mapping (text, one PE per line, or .u32); default: run Waver-HeiProMap")
    ap.add_argument("--rounds", type=int, default=10)
    ap.add_argument("--children", type=int, default=10)
    ap.add_argument("--swaps", type=int, default=3, help="PE content swaps per child")
    ap.add_argument("--moves", type=int, default=30, help="single-row moves per child")
    ap.add_argument("--radius", type=int, default=3, help="max distance (per axis) of a swap")
    ap.add_argument("--margin", type=int, default=2, help="the changes are drawn within this distance of the chain")
    ap.add_argument("--jobs", type=int, default=2, help="concurrent evaluations (about 3.3 GB each at full size)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tail-repair", default="off", help="scheduler tail repair mode (off, gap-walk, ...)")
    ap.add_argument("--tail-repair-passes", type=int, default=40)
    ap.add_argument("--seeds", type=int, default=4, help="scheduler seeds 0..SEEDS-1 per mapping")
    ap.add_argument("--objective", choices=("median", "mean", "max"), default="median",
                    help="statistic of the makespans of the seeds that is minimized")
    ap.add_argument("--holdout", type=int, default=12,
                    help="at the end, M0 and the best mapping are scheduled with this many unseen seeds (1000...)")
    ap.add_argument("--waver-args", default="", help="extra arguments for the Waver-HeiProMap run that makes M0")
    ap.add_argument("--instance", help="instance.json of a smaller instance (default: the full 215 x 215 x 216 one)")
    ap.add_argument("--reconstruct", nargs=4, metavar=("RUN_DIR", "ROUND", "CHILD", "OUT"))
    args = ap.parse_args()
    inst = json.load(open(args.instance)) if args.instance else dict(FULL)
    set_instance(inst)

    if args.reconstruct:
        run_dir, rnd, child, out_path = args.reconstruct
        reconstruct(run_dir, int(rnd), int(child), out_path)
        return

    out = Path(args.out).resolve() if args.out else HEI / "data" / "refine_loop" / time.strftime("run_%Y%m%d_%H%M%S")
    out.mkdir(parents=True, exist_ok=True)
    if args.start:
        args.start = str(Path(args.start).resolve())  # the script changes into the wse-spmv directory
    json.dump(vars(args), open(out / "config.json", "w"), indent=1)
    ram = load_wse()
    ev = Evaluator(ram, out, args.tail_repair, args.tail_repair_passes, eval_template(ram, inst), inst["base"],
                   args.objective)
    seeds = list(range(args.seeds))
    log = open(out / "log.jsonl", "a")

    def record(entry):
        log.write(json.dumps(entry) + "\n")
        log.flush()

    # ---- M0
    if args.start:
        best = read_mapping(args.start)
    else:
        partition_txt = out / "M0_partition.txt"
        cmd = [str(HEI / "build" / "Waver-HeiProMap"),
               "--graph", inst["graph"], "--topology", inst["topology"],
               "--mapping", str(partition_txt), "--imbalance", "0.03", "--config", "fast", "--seed", "0",
               "--threads", "1", "--kaminpar-threads", "1", "--kaminpar-preset", "default",
               "--alpha-send", "116", "--alpha-recv", "111", "--alpha-link", "2", "--beta-send", "0",
               "--beta-recv", "8.0", "--beta-link", "1.16",
               "--prepartitioner", "spectral3d", "--mode", "spectral"] + args.waver_args.split()
        print("making M0 with Waver-HeiProMap ...", flush=True)
        subprocess.run(ram.limited(cmd), check=True, stdout=subprocess.DEVNULL)
        best = read_mapping(partition_txt)
        partition_txt.unlink()
    if len(best) != N:
        sys.exit(f"the mapping has {len(best)} entries, expected {N}")
    best = best.astype(np.uint32)

    best_file = "best_r00.u32"
    best.astype("<u4").tofile(out / best_file)
    t0 = time.time()
    res = ev.run(best, "M0", seeds, want_chain=True)
    best_make = res["makespan"]
    chain = res["chain_pes"]
    record({"round": 0, "child": 0, "parent": best_file, **res, "accepted": True, "seconds": time.time() - t0})
    print(f"M0: {args.objective} makespan {best_make} over seeds {seeds}: {res['makespans']} (res.LB {res['res_lb']}, "
          f"{len(chain)} PEs on the chains)", flush=True)

    rng = np.random.default_rng(args.seed)
    history = [best_make]
    for rnd in range(1, args.rounds + 1):
        t_round = time.time()
        (out / f"r{rnd:02d}").mkdir(exist_ok=True)
        json.dump({"parent": best_file, "chain_pes": chain}, open(out / f"r{rnd:02d}" / "chain.json", "w"))
        order = np.argsort(best, kind="stable")
        counts0 = np.bincount(best, minlength=P).astype(np.int64)
        starts = np.concatenate([[0], np.cumsum(counts0)])
        region = region_of(chain, args.margin)

        children = []
        for c in range(1, args.children + 1):
            child, swaps, moves = make_child(best, order, starts, counts0, region, rng,
                                             args.swaps, args.moves, args.radius)
            diff = np.nonzero(child != best)[0].astype(np.uint32)
            np.savez(out / f"r{rnd:02d}" / f"child_{c:02d}.npz", vertices=diff, pe=child[diff])
            children.append((c, child, swaps, moves, len(diff)))

        def evaluate(item):
            c, child, swaps, moves, changed = item
            t = time.time()
            r = ev.run(child, f"r{rnd:02d}c{c:02d}", seeds)
            return c, child, swaps, moves, changed, r, time.time() - t

        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            results = list(pool.map(evaluate, children))

        winner = min(results, key=lambda x: x[5]["makespan"])
        accepted_child = winner[0] if winner[5]["makespan"] < best_make else None
        for c, child, swaps, moves, changed, r, secs in results:
            record({"round": rnd, "child": c, "parent": best_file, **r, "swaps": swaps, "moves": moves,
                    "changed_rows": changed, "accepted": c == accepted_child, "seconds": secs})
        spread = sorted(x[5]["makespan"] for x in results)
        line = (f"round {rnd}: children {spread[0]} .. {spread[-1]} (median {spread[len(spread) // 2]}), "
                f"parent {best_make}")
        if accepted_child is not None:
            _, child, _, _, _, r, _ = winner
            best, best_make = child, r["makespan"]
            best_file = f"best_r{rnd:02d}.u32"
            best.astype("<u4").tofile(out / best_file)
            chain = ev.run(best, f"r{rnd:02d}best", seeds, want_chain=True)["chain_pes"]
            line += f" -> new best {best_make} (child {accepted_child})"
        else:
            line += " -> no improvement"
        history.append(best_make)
        print(f"{line}  [{time.time() - t_round:.0f}s]", flush=True)

    # schedule dump of the final best mapping, for the analysis
    ev.dump(best, "final", 0, out / "final_viz")

    # the objective was minimized over the seeds 0..SEEDS-1, so check the result on seeds that were never used
    if args.holdout > 0:
        m0 = read_mapping(out / "best_r00.u32").astype(np.uint32)
        held = list(range(1000, 1000 + args.holdout))
        def held_out(item):
            name, part = item
            return name, ev.run(part, "holdout_" + name, held)["makespans"]
        with ThreadPoolExecutor(max_workers=min(2, args.jobs)) as pool:
            hold = dict(pool.map(held_out, [("M0", m0), ("best", best)]))
        json.dump({"seeds": held, **hold}, open(out / "holdout.json", "w"))
        for name, ms in hold.items():
            print(f"held-out seeds {held[0]}..{held[-1]}, {name}: min {min(ms)} median {int(np.median(ms))} "
                  f"mean {int(np.mean(ms))} max {max(ms)}")

    shutil.rmtree(out / "work", ignore_errors=True)
    print(f"{args.objective} makespan per round: {history}")
    print(f"final best: {best_file} ({best_make}); everything is in {out}")


if __name__ == "__main__":
    main()
