#!/usr/bin/env python3
"""Sweep one flow-bisection flag. Quality phase in parallel (seeds), timing phase serial (first seed).

  scripts/bench_flow_param.py --flag --flow-bisection-init-fraction --values 0,0.05,0.1,0.2 [--initial-c 8]
"""
import argparse, copy, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_flow_bisection as b

ap = argparse.ArgumentParser()
ap.add_argument("--flag", required=True)
ap.add_argument("--values", required=True)
ap.add_argument("--configs", default="fast,eco")
ap.add_argument("--graphs", default=",".join(b.DEFAULT_GRAPHS))
ap.add_argument("--seeds", default="0,1,2,3,4")
ap.add_argument("--initial-c", type=int, default=8)
ap.add_argument("--jobs", type=int, default=32)
a = ap.parse_args()
a.hierarchy, a.distance, a.imbalance, a.threads = "4:8:6", "1:10:100", "0.03", 1
a.memlimit, a.no_memlimit, a.out = "100G", False, b.ROOT / "data" / "bench_flow_bisection"
(a.out / "mappings").mkdir(parents=True, exist_ok=True)
vals, configs, graphs = a.values.split(","), a.configs.split(","), a.graphs.split(",")
seeds = [int(s) for s in a.seeds.split(",")]

def job(g, c, v, s):
    aa = copy.copy(a)
    aa.extra = [a.flag, v]
    aa.out = a.out
    return (g, c, "flow", s, aa), v

def run(x):
    j, v = x
    # separate mapping files per value
    j[4].out = a.out
    r = b.run_one(j)
    r["value"] = v
    return r

jobs = [job(g, c, v, s) for c in configs for g in graphs for v in vals for s in seeds]
print(f"quality: {len(jobs)} runs in parallel", flush=True)
with ThreadPoolExecutor(a.jobs) as ex:
    qrows = list(ex.map(run, jobs))
tjobs = [job(g, c, v, seeds[0]) for c in configs for g in graphs for v in vals]
print(f"timing: {len(tjobs)} runs serial", flush=True)
trows = [run(x) for x in tjobs]

def m(rows, c, g, v, k):
    return b.mean([r[k] for r in rows if (r["config"], r["graph"], r["value"]) == (c, g, v) and r["ok"]])

print(f"\n{a.flag} (c={a.initial_c}); QAP mean over seeds {seeds} relative to value {vals[0]}; time/init serial seed {seeds[0]}\n")
print(f"{'config':<7}{'value':>8}{'QAP rel':>9}{'init QAP rel':>13}{'oload':>7}{'init ms':>9}{'total s':>9}")
for c in configs:
    for v in vals:
        qr = [m(qrows, c, g, v, "final_qap") / m(qrows, c, g, vals[0], "final_qap") for g in graphs]
        ir = [m(qrows, c, g, v, "init_qap") / m(qrows, c, g, vals[0], "init_qap") for g in graphs]
        print(f"{c:<7}{v:>8}{b.geomean(qr):>9.3f}{b.geomean(ir):>13.3f}"
              f"{b.mean([m(qrows, c, g, v, 'init_overloaded') for g in graphs]):>7.1f}"
              f"{b.mean([m(trows, c, g, v, 'init_part_ms') for g in graphs]):>9.0f}"
              f"{b.mean([m(trows, c, g, v, 'time_s') for g in graphs]):>9.2f}")
    print()
