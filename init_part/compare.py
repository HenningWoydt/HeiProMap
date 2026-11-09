#!/usr/bin/env python3
"""Compare flow-bisection against the baseline (multisection, kaffpa-strong) on the
coarsest graphs in data/coarsest. Runs are serial (timing). Baseline results are cached in
init_part/baseline_cache.json, so only the first call pays for them (use --refresh to redo).

  python3 init_part/compare.py [--seeds 3] [--graphs wave,144] [--extra "--flow-bisection-seeds 8"]
"""
import argparse, json, math, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build" / "InitPart"
GRAPH_DIR = ROOT / "data" / "coarsest"
CACHE = ROOT / "init_part" / "baseline_cache.json"
COMMON = ["--hierarchy", "4:8:6", "--distance", "1:10:100", "--imbalance", "0.03", "--config", "strong", "--threads", "1"]
BASE = ["--partitioning", "multisection"]
NEW = ["--partitioning", "flow-bisection"]


def run(graph, seed, alg_args, extra):
    cmd = [str(BIN), "--graph", str(graph), "--seed", str(seed)] + COMMON + alg_args + extra
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    def g(key):
        m = re.search(rf"^{re.escape(key)}\s*:\s*(\S+)", out, re.M)
        return float(m.group(1)) if m else math.nan
    return {"qap": g("QAP"), "ms": g("Time (ms)"), "empty": g("#empty partitions"),
            "oload": g("#oload partitions"), "maxw": g("max block w"), "lmax": g("Lmax")}


def mean(xs): return sum(xs) / len(xs)
def geomean(xs): return math.exp(sum(math.log(x) for x in xs) / len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--graphs", default="", help="comma separated names (default: all in data/coarsest)")
    ap.add_argument("--extra", default="", help="extra InitPart args for the new partitioner")
    ap.add_argument("--refresh", action="store_true", help="re-run the baseline")
    a = ap.parse_args()
    extra = a.extra.split()

    names = [n for n in a.graphs.split(",") if n] or sorted({p.name.rsplit("_s", 1)[0] for p in GRAPH_DIR.glob("*_s0.graph")})
    cache = {} if a.refresh or not CACHE.exists() else json.loads(CACHE.read_text())

    rows = []
    for name in names:
        graph = GRAPH_DIR / f"{name}_s0.graph"  # coarsening is deterministic: one file per graph
        key = lambda s: f"{name}|{s}"
        base, new = [], []
        for s in range(a.seeds):
            if key(s) not in cache:
                cache[key(s)] = run(graph, s, BASE, [])
                CACHE.write_text(json.dumps(cache, indent=1))
            base.append(cache[key(s)])
            new.append(run(graph, s, NEW, extra))
            print(f"  {name} seed {s}: base {base[-1]['qap']:.0f} ({base[-1]['ms']:.0f} ms)  new {new[-1]['qap']:.0f} ({new[-1]['ms']:.0f} ms)", file=sys.stderr, flush=True)
        rows.append((name, base, new))

    # per graph: mean over the seeds; total row: geometric mean over the graphs
    # balance = max block weight / Lmax (<= 1 is feasible), oload = number of overloaded blocks
    hdr = f"{'graph':15} {'base QAP':>10} {'new QAP':>10} {'ratio':>6} | {'base ms':>8} {'new ms':>8} {'speedup':>7} | {'base bal':>8} {'new bal':>8} {'base ol':>7} {'new ol':>6}"
    print("\n" + hdr)
    print("-" * len(hdr))
    cols = {k: [] for k in ("bq", "nq", "qr", "bt", "nt", "sp", "bb", "nb")}
    for name, base, new in rows:
        bq, nq = mean([r["qap"] for r in base]), mean([r["qap"] for r in new])
        bt, nt = mean([r["ms"] for r in base]), mean([r["ms"] for r in new])
        bb = mean([r["maxw"] / r["lmax"] for r in base]); nb = mean([r["maxw"] / r["lmax"] for r in new])
        bo = mean([r["oload"] for r in base]); no = mean([r["oload"] for r in new])
        for k, v in zip(cols, (bq, nq, nq / bq, bt, nt, bt / nt, bb, nb)): cols[k].append(v)
        print(f"{name:15} {bq:10.0f} {nq:10.0f} {nq/bq:6.3f} | {bt:8.0f} {nt:8.0f} {bt/nt:6.2f}x | {bb:8.3f} {nb:8.3f} {bo:7.1f} {no:6.1f}")
    print("-" * len(hdr))
    G = {k: geomean(v) for k, v in cols.items()}
    print(f"{'geomean':15} {G['bq']:10.0f} {G['nq']:10.0f} {G['qr']:6.3f} | {G['bt']:8.0f} {G['nt']:8.0f} {G['sp']:6.2f}x | {G['bb']:8.3f} {G['nb']:8.3f}")
    print(f"(mean over {a.seeds} seeds per graph, geomean over {len(rows)} graphs; ratio < 1: new has the better QAP; speedup > 1: new is faster; bal = max block weight / Lmax)")


if __name__ == "__main__":
    main()
