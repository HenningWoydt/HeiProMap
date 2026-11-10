#!/usr/bin/env python3
"""Full strong pipeline: multisection (cached) vs flow-bisection, serial.

  python3 init_part/compare_full.py [--seeds 3] [--graphs wave,144] [--extra "--flow-bisection-seeds 8"] [--refresh]
Prints mean final QAP and total time per graph and a geomean row.
"""
import argparse, json, math, re, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build" / "HeiProMap"
DATA = ROOT.parent / "ProMapRepo" / "data" / "mapping"
FILES = {"144": "144.graph", "2cubes_sphere": "2cubes_sphere.mtx.graph", "598a": "598a.graph",
         "G2_circuit": "G2_circuit.mtx.graph", "cont-300": "cont-300.mtx.graph", "cop20k_A": "cop20k_A.mtx.graph",
         "fe_ocean": "fe_ocean.graph", "fe_rotor": "fe_rotor.graph", "thermomech_TC": "thermomech_TC.mtx.graph",
         "wave": "wave.graph"}
for _n in ["Dubcova3", "bmwcra_1", "boneS01", "cfd2", "shipsec5"]:
    FILES[_n] = _n + ".mtx.graph"
for _n in ["m14b", "auto", "afshell9", "ecology2", "thermal2", "G3_circuit", "nlr"]:
    FILES[_n] = _n + ".graph"
COMMON = ["--hierarchy", "4:8:6", "--distance", "1:10:100", "--imbalance", "0.03", "--threads", "1"]


def run(graph, seed, alg, extra, profile="strong"):
    with tempfile.NamedTemporaryFile(suffix=".txt") as tmp:
        cmd = [str(BIN), "--graph", str(DATA / FILES[graph]), "--mapping", tmp.name, "--seed", str(seed),
               "--partitioning", alg, "--config", profile] + COMMON + extra
        out = subprocess.run(cmd, capture_output=True, text=True).stdout
    q = re.search(r"Final QAP\s*:\s*(\d+)", out)
    t = re.search(r"Total solve time\s*:\s*([\d.]+)", out)
    i = re.search(r"Init\. Part\.\s*:\s*([\d.]+)", out)
    if not q or not t:
        sys.exit(f"failed: {' '.join(cmd)}\n{out[-500:]}")
    return {"qap": int(q.group(1)), "ms": float(t.group(1)), "init_ms": float(i.group(1)) if i else 0.0}


def geo(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--graphs", default="")
    ap.add_argument("--extra", default="")
    ap.add_argument("--profile", default="strong")
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    CACHE = Path(__file__).resolve().parent / ("full_cache.json" if a.profile == "strong" else f"full_cache_{a.profile}.json")
    graphs = a.graphs.split(",") if a.graphs else list(FILES)
    extra = a.extra.split()
    cache = {} if a.refresh or not CACHE.exists() else json.loads(CACHE.read_text())
    rows = []
    for g in graphs:
        b, n = [], []
        for s in range(a.seeds):
            key = f"{g}|{s}"
            if key not in cache:
                cache[key] = run(g, s, "multisection", [], a.profile)
                CACHE.write_text(json.dumps(cache))
            b.append(cache[key])
            n.append(run(g, s, "flow-bisection", extra, a.profile))
            print(f"  {g} seed {s}: base {b[-1]['qap']} ({b[-1]['ms']:.0f} ms)  new {n[-1]['qap']} ({n[-1]['ms']:.0f} ms)", flush=True)
        m = lambda L, k: sum(x[k] for x in L) / len(L)
        rows.append((g, m(b, "qap"), m(n, "qap"), m(b, "ms"), m(n, "ms"), m(b, "init_ms"), m(n, "init_ms")))
    print(f"\n{'graph':<15}{'base QAP':>11}{'new QAP':>11}{'ratio':>7} |{'base s':>8}{'new s':>8}{'speedup':>8} |{'init base':>10}{'init new':>9}")
    print("-" * 95)
    for g, bq, nq, bt, nt, bi, ni in rows:
        print(f"{g:<15}{bq:>11.0f}{nq:>11.0f}{nq / bq:>7.3f} |{bt / 1000:>8.1f}{nt / 1000:>8.1f}{bt / nt:>7.2f}x |{bi / 1000:>10.1f}{ni / 1000:>9.1f}")
    print("-" * 95)
    gb, gn = geo([r[1] for r in rows]), geo([r[2] for r in rows])
    gbt, gnt = geo([r[3] for r in rows]), geo([r[4] for r in rows])
    print(f"{'geomean':<15}{gb:>11.0f}{gn:>11.0f}{gn / gb:>7.3f} |{gbt / 1000:>8.1f}{gnt / 1000:>8.1f}{gbt / gnt:>7.2f}x |{geo([r[5] for r in rows]) / 1000:>10.1f}{geo([r[6] for r in rows]) / 1000:>9.1f}")


main()
