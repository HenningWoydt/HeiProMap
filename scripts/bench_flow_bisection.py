#!/usr/bin/env python3
"""Small benchmark suite: flow-bisection vs the current initial partitioning (multisection).

Two phases, for the configs fast/eco/strong:
  1. QUALITY  - all runs in parallel (default 32 single-threaded runs at once). Timings are
                meaningless here, only the QAP counts. A config passes if the geometric mean of
                final-QAP(flow) / final-QAP(current) is <= 1 + tolerance (default 5%).
  2. TIMING   - only for passing configs: runs one after the other (nothing else running).
                Results of the current algorithm are cached in <out>/cache_current.json and reused.

Usage:
  scripts/bench_flow_bisection.py [--graphs a,b,..] [--configs fast,eco,strong] [--seeds 0,1,2]
      [--timing-seeds 0] [--hierarchy 4:8:6] [--distance 1:10:100] [--threads 1] [--initial-c N]
      [--jobs 32] [--tolerance 0.05] [--quality-only] [--force-timing] [--refresh-cache]
Note: the cache of the current algorithm does not notice code changes of HeiProMap itself, use
--refresh-cache after changing anything that affects the current algorithm.
"""
import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GRAPH_DIR = Path(os.environ.get("GRAPH_DIR", ROOT.parent / "ProMapRepo" / "data" / "mapping"))
DEFAULT_GRAPHS = ["fe_ocean", "cont-300.mtx", "598a", "wave", "m14b"]
ALGS = {"current": "multisection", "flow": "flow-bisection"}
ANSI = re.compile(r"\x1b\[[0-9;]*m")
FIELDS = ["graph", "config", "alg", "seed", "ok", "init_qap", "final_qap", "init_part_ms",
          "init_overloaded", "imb", "time_s"]


def run_one(job):
    graph, config, alg, seed, a = job
    mapping = a.out / "mappings" / f"{graph}_{config}_{alg}_{seed}.txt"
    cmd = [str(ROOT / "build" / "HeiProMap"), "--graph", str(GRAPH_DIR / f"{graph}.graph"),
           "--mapping", str(mapping), "--hierarchy", a.hierarchy, "--distance", a.distance,
           "--imbalance", a.imbalance, "--config", config, "--seed", str(seed),
           "--threads", str(a.threads)]
    if alg in ALGS and not getattr(a, "default_partitioning", False):
        cmd += ["--partitioning", ALGS[alg]]
    if a.initial_c:
        cmd += ["--initial-c", str(a.initial_c)]
    cmd += getattr(a, "extra", [])
    if not a.no_memlimit:
        cmd = ["systemd-run", "--user", "--scope", "-q", "-p", f"MemoryMax={a.memlimit}"] + cmd
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = ANSI.sub("", p.stdout)
    row = {"graph": graph, "config": config, "alg": alg, "seed": seed, "ok": p.returncode == 0}

    def grab(pattern):
        m = re.search(pattern, out)
        return float(m.group(1)) if m else math.nan

    row["init_qap"] = grab(r"Init\. QAP\s*:\s*([\d.]+)")
    row["final_qap"] = grab(r"Final QAP\s*:\s*([\d.]+)")
    row["init_part_ms"] = grab(r"Init\. Part\.\s*:\s*([\d.]+)")
    row["init_overloaded"] = grab(r"Init\. #oload partitions\s*:\s*([\d.]+)")
    row["time_s"] = grab(r"Total Time in HeiProMap\.cpp \(s\):\s*([\d.]+)")
    m = re.search(r"^\s*0\s*\|[^|]*\|[^|]*\|[^|]*\|[^|]*\|[^|]*\|\s*([\d.]+)", out, re.M)
    row["imb"] = float(m.group(1)) if m else math.nan
    if not row["ok"] or math.isnan(row["final_qap"]):
        row["ok"] = False
        sys.stderr.write(f"FAILED {graph} {config} {alg} seed={seed}: {(p.stdout + p.stderr)[-200:]}\n")
    return row


def mean(xs):
    xs = [x for x in xs if not math.isnan(x)]
    return sum(xs) / len(xs) if xs else math.nan


def geomean(xs):
    xs = [x for x in xs if x > 0 and not math.isnan(x)]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else math.nan


def pick(rows, graph, config, alg, key, seeds=None):
    return mean([r[key] for r in rows if r["graph"] == graph and r["config"] == config
                 and r["alg"] == alg and r["ok"] and (seeds is None or r["seed"] in seeds)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", default=",".join(DEFAULT_GRAPHS))
    ap.add_argument("--configs", default="fast,eco,strong")
    ap.add_argument("--seeds", default="0,1,2,3,4", help="seeds of the quality phase")
    ap.add_argument("--timing-seeds", default="0,1,2,3,4", help="seeds of the (serial) timing phase")
    ap.add_argument("--hierarchy", default="4:8:6")
    ap.add_argument("--distance", default="1:10:100")
    ap.add_argument("--imbalance", default="0.03")
    ap.add_argument("--threads", type=int, default=1, help="threads of each HeiProMap run")
    ap.add_argument("--initial-c", type=int, default=0, help="0 = config default")
    ap.add_argument("--jobs", type=int, default=32, help="parallel runs in the quality phase")
    ap.add_argument("--tolerance", type=float, default=0.05)
    ap.add_argument("--quality-only", action="store_true")
    ap.add_argument("--force-timing", action="store_true", help="time all configs even if quality is bad")
    ap.add_argument("--refresh-cache", action="store_true")
    ap.add_argument("--memlimit", default="100G")
    ap.add_argument("--no-memlimit", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "bench_flow_bisection")
    a = ap.parse_args()
    (a.out / "mappings").mkdir(parents=True, exist_ok=True)
    graphs, configs = a.graphs.split(","), a.configs.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    tseeds = [int(s) for s in a.timing_seeds.split(",")]

    # cache of serial runs of the current algorithm
    cache_file = a.out / "cache_current.json"
    cache = {} if a.refresh_cache or not cache_file.exists() else json.loads(cache_file.read_text())

    def key(r):
        return "|".join(map(str, [r["graph"], r["config"], r["seed"], a.hierarchy, a.distance,
                                  a.imbalance, a.threads, a.initial_c]))

    def save_cache():
        cache_file.write_text(json.dumps(cache, indent=1))

    # ---------------- phase 1: quality ----------------
    rows, todo = [], []
    for g in graphs:
        for c in configs:
            for s in seeds:
                todo.append((g, c, "flow", s, a))
                r = dict(graph=g, config=c, seed=s)
                if key(r) in cache:
                    rows.append(cache[key(r)])
                else:
                    todo.append((g, c, "current", s, a))
    print(f"PHASE 1 (quality): {len(todo)} runs ({len(rows)} current runs cached), {a.jobs} in parallel", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rows += list(ex.map(run_one, todo))

    passed = []
    print(f"\nfinal QAP, mean over seeds {seeds}; ratio = flow / current\n")
    print(f"{'config':<8}{'graph':<14}{'current':>12}{'flow':>12}{'ratio':>8}{'init ratio':>11}")
    for c in configs:
        ratios, iratios = [], []
        for g in graphs:
            qc, qf = pick(rows, g, c, "current", "final_qap"), pick(rows, g, c, "flow", "final_qap")
            ic, if_ = pick(rows, g, c, "current", "init_qap"), pick(rows, g, c, "flow", "init_qap")
            ratios.append(qf / qc)
            iratios.append(if_ / ic)
            print(f"{c:<8}{g:<14}{qc:>12.0f}{qf:>12.0f}{qf / qc:>8.3f}{if_ / ic:>11.3f}")
        gm = geomean(ratios)
        ok = gm <= 1 + a.tolerance
        if ok:
            passed.append(c)
        print(f"{c:<8}{'GEOMEAN':<14}{'':>12}{'':>12}{gm:>8.3f}{geomean(iratios):>11.3f}   "
              f"{'PASS' if ok else 'FAIL'} (limit {1 + a.tolerance:.2f})\n")
    with open(a.out / "quality_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    if a.quality_only:
        return
    timing_configs = configs if a.force_timing else passed
    if not timing_configs:
        print("no config passed the quality gate, skipping the timing phase (--force-timing to override)")
        return

    # ---------------- phase 2: timing, strictly serial ----------------
    print(f"PHASE 2 (timing, serial): configs {timing_configs}, seeds {tseeds}", flush=True)
    trows = []
    for c in timing_configs:
        for g in graphs:
            for s in tseeds:
                r = dict(graph=g, config=c, seed=s)
                if key(r) in cache:
                    trows.append(cache[key(r)])
                else:
                    print(f"  current {c} {g} seed {s} ...", flush=True)
                    row = run_one((g, c, "current", s, a))
                    if row["ok"]:
                        cache[key(r)] = row
                        save_cache()
                    trows.append(row)
                print(f"  flow    {c} {g} seed {s} ...", flush=True)
                trows.append(run_one((g, c, "flow", s, a)))

    print(f"\ntiming (serial), seeds {tseeds}; ratio = flow / current\n")
    print(f"{'config':<8}{'graph':<14}{'init ms cur':>12}{'init ms flow':>13}{'total s cur':>12}"
          f"{'total s flow':>13}{'ratio':>7}{'oload cur':>10}{'oload flow':>11}")
    for c in timing_configs:
        tr = []
        for g in graphs:
            p = lambda alg, k: pick(trows, g, c, alg, k, tseeds)
            tr.append(p("flow", "time_s") / p("current", "time_s"))
            print(f"{c:<8}{g:<14}{p('current', 'init_part_ms'):>12.0f}{p('flow', 'init_part_ms'):>13.0f}"
                  f"{p('current', 'time_s'):>12.2f}{p('flow', 'time_s'):>13.2f}{tr[-1]:>7.2f}"
                  f"{p('current', 'init_overloaded'):>10.0f}{p('flow', 'init_overloaded'):>11.0f}")
        print(f"{c:<8}{'GEOMEAN':<14}{'':>12}{'':>13}{'':>12}{'':>13}{geomean(tr):>7.2f}\n")
    with open(a.out / "timing_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(trows)
    print(f"results in {a.out}")


if __name__ == "__main__":
    main()
