#!/usr/bin/env python3
"""Compares Waver-HeiProMap configurations on a small instance with the wse-spmv scheduler.

The makespan of one mapping varies by several percent with the scheduler seed, so every configuration is run with
MAPPER_SEEDS Waver seeds and every resulting mapping is scheduled with SCHED_SEEDS scheduler seeds; the value of a
mapping is the median over its scheduler seeds. The configurations are NAME=ARGS pairs (extra Waver arguments).

  python3 tools/bench_small.py --instance data/small100/instance.json \
      "fast=" "eco=--config eco" "smooth31=--spectral-big-rounds 1 --spectral-big-depth 8 --spectral-big-n 50000"
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import refine_loop as rl  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("configs", nargs="+", help="NAME=ARGS, extra arguments of Waver-HeiProMap (the first one is the baseline)")
    ap.add_argument("--instance", default=str(rl.HEI / "data" / "small100" / "instance.json"))
    ap.add_argument("--mapper-seeds", type=int, default=4)
    ap.add_argument("--sched-seeds", type=int, default=12)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--out", help="directory for the results (default data/bench_small/<time>)")
    ap.add_argument("--tail-repair", default="off")
    ap.add_argument("--tail-repair-passes", type=int, default=40)
    args = ap.parse_args()

    inst = json.load(open(args.instance))
    rl.set_instance(inst)
    out = Path(args.out).resolve() if args.out else rl.HEI / "data" / "bench_small" / time.strftime("%Y%m%d_%H%M%S")
    out.mkdir(parents=True, exist_ok=True)
    ram = rl.load_wse()
    ev = rl.Evaluator(ram, out, args.tail_repair, args.tail_repair_passes, rl.eval_template(ram, inst), inst["base"])
    sched_seeds = list(range(1000, 1000 + args.sched_seeds))
    configs = [(c.split("=", 1)[0], c.split("=", 1)[1].split()) for c in args.configs]

    def task(job):
        name, extra, seed = job
        tag = f"{name}_s{seed}"
        partition = out / "work" / f"{tag}.txt"
        cmd = [str(rl.HEI / "build" / "Waver-HeiProMap"),
               "--graph", inst["graph"], "--topology", inst["topology"],
               "--mapping", str(partition), "--imbalance", "0.03", "--config", "fast", "--seed", str(seed),
               "--threads", "1", "--kaminpar-threads", "1", "--kaminpar-preset", "default",
               "--alpha-send", "116", "--alpha-recv", "111", "--alpha-link", "2", "--beta-send", "0",
               "--beta-recv", "8.0", "--beta-link", "1.16",
               "--prepartitioner", "spectral3d", "--mode", "spectral"] + extra
        t = time.time()
        subprocess.run(ram.limited(cmd), check=True, stdout=subprocess.DEVNULL)
        solve = time.time() - t
        part = rl.read_mapping(partition).astype(np.uint32)
        partition.unlink()
        res = ev.run(part, tag, sched_seeds)
        return name, seed, solve, res

    (out / "work").mkdir(exist_ok=True)
    jobs = [(n, e, s) for n, e in configs for s in range(args.mapper_seeds)]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(task, jobs))
    shutil.rmtree(out / "work", ignore_errors=True)
    json.dump([{"config": n, "mapper_seed": s, "solve_seconds": t, "makespans": r["makespans"], "res_lb": r["res_lb"],
                "hop_volume": r["hop_volume"], "max_link": r["max_link"]} for n, s, t, r in results],
              open(out / "results.json", "w"), indent=1)

    base = configs[0][0]
    medians = {n: np.array([np.median(r["makespans"]) for nn, s, t, r in results if nn == n]) for n, _ in configs}
    pooled = {n: np.concatenate([r["makespans"] for nn, s, t, r in results if nn == n]) for n, _ in configs}
    print(f"{'config':24s} {'solve s':>8s} {'median of medians':>18s} {'per mapper seed (median over sched seeds)':>46s} "
          f"{'pooled mean':>12s} {'vs ' + base:>12s}")
    for n, _ in configs:
        solve = np.mean([t for nn, s, t, r in results if nn == n])
        diff = 100.0 * (np.mean(pooled[n]) / np.mean(pooled[base]) - 1.0)
        # standard error of the difference of the pooled means (all runs treated as independent)
        se = 100.0 * np.sqrt(np.var(pooled[n], ddof=1) / len(pooled[n]) + np.var(pooled[base], ddof=1) / len(pooled[base])) \
            / np.mean(pooled[base])
        per_seed = " ".join(f"{int(m)}" for m in medians[n])
        print(f"{n:24s} {solve:8.1f} {int(np.median(medians[n])):18d} {per_seed:>46s} {np.mean(pooled[n]):12.0f} "
              f"{diff:+6.1f}% +-{se:.1f}")
    print(f"results in {out}")


if __name__ == "__main__":
    main()
