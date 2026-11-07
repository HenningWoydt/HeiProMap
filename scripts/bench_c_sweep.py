#!/usr/bin/env python3
"""Serial sweep: flow-bisection with one initial_c vs the current algorithm for several initial_c.

Everything runs one after the other (timings are meaningful). Every finished run is stored in
<out>/sweep_cache.json (resumable; reruns skip finished runs) and all rows go to <out>/sweep_runs.csv.

  scripts/bench_c_sweep.py [--flow-c 8] [--current-c 4,8,16,32,64] [--seeds 0]
                           [--graphs ..] [--configs fast,eco,strong] [--refresh]
"""
import argparse
import copy
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_flow_bisection as b  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", default=",".join(b.DEFAULT_GRAPHS))
    ap.add_argument("--configs", default="fast,eco,strong")
    ap.add_argument("--seeds", default="0,1,2,3,4")
    ap.add_argument("--flow-c", type=int, default=8)
    ap.add_argument("--current-c", default="4,8,16,32,64")
    ap.add_argument("--hierarchy", default="4:8:6")
    ap.add_argument("--distance", default="1:10:100")
    ap.add_argument("--imbalance", default="0.03")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--memlimit", default="100G")
    ap.add_argument("--no-memlimit", action="store_true")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--out", type=Path, default=b.ROOT / "data" / "bench_flow_bisection")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "mappings").mkdir(exist_ok=True)
    graphs, configs = a.graphs.split(","), a.configs.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    cs_current = [int(c) for c in a.current_c.split(",")]

    cache_file = a.out / "sweep_cache.json"
    cache = {} if a.refresh or not cache_file.exists() else json.loads(cache_file.read_text())
    jobs = [(g, c, "flow", s, a.flow_c) for c in configs for g in graphs for s in seeds]
    jobs += [(g, c, "current", s, ic) for c in configs for g in graphs for s in seeds for ic in cs_current]

    rows = []
    for i, (g, c, alg, s, ic) in enumerate(jobs, 1):
        key = "|".join(map(str, [g, c, alg, s, ic, a.hierarchy, a.distance, a.imbalance, a.threads]))
        if key not in cache:
            print(f"[{i}/{len(jobs)}] {alg} c={ic} {c} {g} seed {s}", flush=True)
            aa = copy.copy(a)
            aa.initial_c = ic
            row = b.run_one((g, c, alg, s, aa))
            if not row["ok"]:
                continue
            row["initial_c"] = ic
            cache[key] = row
            cache_file.write_text(json.dumps(cache, indent=1))
        rows.append(cache[key])

    fields = b.FIELDS + ["initial_c"]
    with open(a.out / "sweep_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    def get(g, c, alg, ic, k):
        v = [r[k] for r in rows if (r["graph"], r["config"], r["alg"], r["initial_c"]) == (g, c, alg, ic)
             and r["seed"] in seeds]
        return b.mean(v)

    for c in configs:
        print(f"\n=== config {c}: final QAP / total time (s), seeds {seeds} ===")
        head = f"{'graph':<14}{'flow c=' + str(a.flow_c):>20}" + "".join(f"{'cur c=' + str(x):>20}" for x in cs_current)
        print(head)
        gm_q = {x: [] for x in cs_current}
        gm_t = {x: [] for x in cs_current}
        for g in graphs:
            fq, ft = get(g, c, "flow", a.flow_c, "final_qap"), get(g, c, "flow", a.flow_c, "time_s")
            line = f"{g:<14}{fq:>12.0f}{ft:>7.1f}s "
            for x in cs_current:
                q, t = get(g, c, "current", x, "final_qap"), get(g, c, "current", x, "time_s")
                gm_q[x].append(fq / q)
                gm_t[x].append(ft / t)
                line += f"{q:>12.0f}{t:>7.1f}s "
            print(line)
        print(f"{'flow/cur QAP':<14}{'':>20}" + "".join(f"{b.geomean(gm_q[x]):>20.3f}" for x in cs_current))
        print(f"{'flow/cur time':<14}{'':>20}" + "".join(f"{b.geomean(gm_t[x]):>20.2f}" for x in cs_current))
    print(f"\nsaved: {a.out / 'sweep_runs.csv'}  (cache: {cache_file})")


if __name__ == "__main__":
    main()
