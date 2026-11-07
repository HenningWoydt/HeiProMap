#!/usr/bin/env python3
"""Quick serial benchmark: default HeiProMap with 1 thread, regular refiners vs --force-parallel-suite.

  scripts/bench_parallel_suite.py [--graphs ..] [--configs fast,eco,strong] [--seeds 0]
"""
import argparse
import copy
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_flow_bisection as b  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--graphs", default=",".join(b.DEFAULT_GRAPHS))
ap.add_argument("--configs", default="fast,eco,strong")
ap.add_argument("--seeds", default="0,1,2,3,4")
ap.add_argument("--hierarchy", default="4:8:6")
ap.add_argument("--distance", default="1:10:100")
ap.add_argument("--imbalance", default="0.03")
ap.add_argument("--memlimit", default="100G")
ap.add_argument("--no-memlimit", action="store_true")
ap.add_argument("--out", type=Path, default=b.ROOT / "data" / "bench_flow_bisection")
a = ap.parse_args()
a.threads, a.initial_c, a.default_partitioning = 1, 0, True
(a.out / "mappings").mkdir(parents=True, exist_ok=True)
graphs, configs, seeds = a.graphs.split(","), a.configs.split(","), [int(s) for s in a.seeds.split(",")]
variants = {"regular": [], "parallel-suite": ["--force-parallel-suite", "true"]}

rows = []
for c in configs:
    for g in graphs:
        for s in seeds:
            for name, extra in variants.items():
                print(f"{c} {g} seed {s} {name}", flush=True)
                aa = copy.copy(a)
                aa.extra = extra
                # run_one takes alg only for the file name / partitioning flag; "default" keeps defaults
                r = b.run_one((g, c, "default-" + name, s, aa))
                r["variant"] = name
                rows.append(r)
with open(a.out / "parallel_suite.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=b.FIELDS + ["variant"])
    w.writeheader()
    w.writerows(rows)

def m(g, c, v, k):
    return b.mean([r[k] for r in rows if (r["graph"], r["config"], r["variant"]) == (g, c, v) and r["ok"]])

print("\nfinal QAP / total time (s), 1 thread; ratio = parallel-suite / regular\n")
print(f"{'config':<8}{'graph':<14}{'QAP reg':>11}{'QAP suite':>11}{'ratio':>7}{'t reg':>8}{'t suite':>9}{'ratio':>7}")
for c in configs:
    qr, tr = [], []
    for g in graphs:
        a1, a2, t1, t2 = m(g, c, "regular", "final_qap"), m(g, c, "parallel-suite", "final_qap"), \
            m(g, c, "regular", "time_s"), m(g, c, "parallel-suite", "time_s")
        qr.append(a2 / a1); tr.append(t2 / t1)
        print(f"{c:<8}{g:<14}{a1:>11.0f}{a2:>11.0f}{a2 / a1:>7.3f}{t1:>8.2f}{t2:>9.2f}{t2 / t1:>7.2f}")
    print(f"{c:<8}{'GEOMEAN':<14}{'':>11}{'':>11}{b.geomean(qr):>7.3f}{'':>8}{'':>9}{b.geomean(tr):>7.2f}\n")
