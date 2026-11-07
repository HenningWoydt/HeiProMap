#!/usr/bin/env python3
"""Parallel quality sweep of flow-bisection over initial_c (timings are NOT meaningful).

Compares final QAP against the cached serial runs of the current algorithm (sweep_cache.json,
written by bench_c_sweep.py) at --ref-c. Rows go to <out>/flow_c_sweep.csv.

  scripts/bench_flow_c_sweep.py [--flow-c 4,8,16,32,64] [--ref-c 64] [--seeds 0] [--jobs 32]
"""
import argparse
import copy
import shlex
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_flow_bisection as b  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", default=",".join(b.DEFAULT_GRAPHS))
    ap.add_argument("--configs", default="fast,eco,strong")
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--flow-c", default="4,8,16,32,64")
    ap.add_argument("--ref-c", type=int, default=64, help="initial_c of the current algorithm to compare to")
    ap.add_argument("--jobs", type=int, default=32)
    ap.add_argument("--extra", default="", help="extra HeiProMap flags, e.g. \"--flow-bisection-init-fraction 0\"")
    ap.add_argument("--tag", default="", help="suffix of the csv file name")
    ap.add_argument("--hierarchy", default="4:8:6")
    ap.add_argument("--distance", default="1:10:100")
    ap.add_argument("--imbalance", default="0.03")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--memlimit", default="100G")
    ap.add_argument("--no-memlimit", action="store_true")
    ap.add_argument("--out", type=Path, default=b.ROOT / "data" / "bench_flow_bisection")
    a = ap.parse_args()
    (a.out / "mappings").mkdir(parents=True, exist_ok=True)
    graphs, configs = a.graphs.split(","), a.configs.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    cs = [int(c) for c in a.flow_c.split(",")]

    jobs = []
    for c in configs:
        for g in graphs:
            for s in seeds:
                for ic in cs:
                    aa = copy.copy(a)
                    aa.initial_c = ic
                    aa.extra = shlex.split(a.extra)
                    # separate mapping file per c
                    aa.out = a.out
                    jobs.append((g, c, "flow", s, aa))
    print(f"{len(jobs)} runs, {a.jobs} in parallel", flush=True)

    def run(job):
        row = b.run_one(job)
        row["initial_c"] = job[4].initial_c
        return row

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rows = list(ex.map(run, jobs))
    with open(a.out / f"flow_c_sweep{a.tag}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=b.FIELDS + ["initial_c"])
        w.writeheader()
        w.writerows(rows)

    cache_file = a.out / "sweep_cache.json"
    cache = json.loads(cache_file.read_text()) if cache_file.exists() else {}

    def cur(g, c, s, k):
        key = "|".join(map(str, [g, c, "current", s, a.ref_c, a.hierarchy, a.distance, a.imbalance, a.threads]))
        return cache[key][k] if key in cache else float("nan")

    def flow(g, c, ic, k):
        return b.mean([r[k] for r in rows if (r["graph"], r["config"], r["initial_c"]) == (g, c, ic)
                       and r["ok"]])

    for c in configs:
        print(f"\n=== {c}: flow final QAP per initial_c, and ratio to current c={a.ref_c} (seeds {seeds}) ===")
        print(f"{'graph':<14}{'current':>12}" + "".join(f"{'flow c=' + str(x):>14}" for x in cs))
        ratios = {x: [] for x in cs}
        for g in graphs:
            ref = b.mean([cur(g, c, s, "final_qap") for s in seeds])
            line = f"{g:<14}{ref:>12.0f}"
            for x in cs:
                q = flow(g, c, x, "final_qap")
                ratios[x].append(q / ref)
                line += f"{q:>14.0f}"
            print(line)
        print(f"{'ratio geomean':<14}{'':>12}" + "".join(f"{b.geomean(ratios[x]):>14.3f}" for x in cs))
        print(f"{'init QAP mean':<14}{'':>12}" + "".join(
            f"{b.mean([flow(g, c, x, 'init_qap') for g in graphs]):>14.0f}" for x in cs))
        print(f"{'oload blocks':<14}{'':>12}" + "".join(
            f"{b.mean([flow(g, c, x, 'init_overloaded') for g in graphs]):>14.1f}" for x in cs))
    print(f"\nsaved: {a.out / ('flow_c_sweep' + a.tag + '.csv')}")


if __name__ == "__main__":
    main()
