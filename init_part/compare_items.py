#!/usr/bin/env python3
"""Compare flow-bisection with the kaffpa-strong multisection item by item.

Every item in data/multisection_items/<graph>/item_*/ is one kaffpa call (graph.graph, stats.json).
Flow-bisection solves the same item as a flat k-way split (InitPart --hierarchy k --distance 1: QAP = 2 * edge cut)
with the item's k and allowed imbalance. Compared: edge cut, balance (max block weight / allowed), time.
Per depth and per graph: geometric mean of the ratios new/kaffpa. Writes init_part/items_result.csv.

  python3 init_part/compare_items.py [--graphs wave,144] [--seed 0] [--extra "--flow-bisection-seeds 12"] [--worst 10]
"""
import argparse, csv, json, math, re, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "build" / "InitPart"
ITEMS = ROOT / "data" / "multisection_items"


def run(graph, k, imb, seed, extra):
    cmd = [str(BIN), "--graph", str(graph), "--hierarchy", str(k), "--distance", "1", "--imbalance", repr(imb),
           "--config", "strong", "--partitioning", "flow-bisection", "--seed", str(seed), "--threads", "1"] + extra
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    def g(key):
        m = re.search(rf"^{re.escape(key)}\s*:\s*(\S+)", out, re.M)
        return float(m.group(1)) if m else math.nan
    return {"cut": g("QAP") / 2, "maxw": g("max block w"), "lmax": g("Lmax"), "ms": g("Time (ms)"), "oload": g("#oload partitions")}


def geomean(xs):
    xs = [x for x in xs if x > 0 and not math.isnan(x)]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else math.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", default="")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--extra", default="")
    ap.add_argument("--best-of", type=int, default=1, help="run seeds seed..seed+N-1 and keep the smallest cut (like the threads of the solver)")
    ap.add_argument("--worst", type=int, default=10)
    a = ap.parse_args()
    extra = a.extra.split()
    names = [n for n in a.graphs.split(",") if n] or sorted(p.name for p in ITEMS.iterdir() if p.is_dir())

    rows = []
    for name in names:
        for d in sorted((ITEMS / name).glob("item_*")):
            st = json.loads((d / "stats.json").read_text())
            rs = [run(d / "graph.graph", st["k"], st["imbalance_allowed"], a.seed + i, extra) for i in range(a.best_of)]
            r = dict(min(rs, key=lambda x: (x["oload"] > 0, x["cut"])), ms=sum(x["ms"] for x in rs))  # time of all runs
            rows.append({"graph": name, "item": d.name, "depth": st["depth"], "n": st["n"], "k": st["k"],
                         "b_cut": st["edge_cut"], "n_cut": r["cut"],
                         "b_bal": st["max_block_weight"] / st["allowed_block_weight"], "n_bal": r["maxw"] / r["lmax"],
                         "b_ms": st["time_ms"], "n_ms": r["ms"], "n_oload": r["oload"]})
        print(f"  {name} done", flush=True)

    with open(ROOT / "init_part" / "items_result.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    def summary(label, sel):
        if not sel: return
        cr = geomean([r["n_cut"] / r["b_cut"] for r in sel if r["b_cut"] > 0])
        sp = geomean([r["b_ms"] / r["n_ms"] for r in sel])
        wins = sum(r["n_cut"] < r["b_cut"] for r in sel); ties = sum(r["n_cut"] == r["b_cut"] for r in sel)
        infeas = sum(r["n_bal"] > 1.0 for r in sel)
        bb = sum(r["b_bal"] for r in sel) / len(sel); nb = sum(r["n_bal"] for r in sel) / len(sel)
        print(f"{label:16} {len(sel):5} {cr:9.3f} {sp:8.1f}x {wins:4}/{ties:<4} {bb:8.3f} {nb:8.3f} {infeas:5}")

    hdr = f"{'':16} {'items':>5} {'cut ratio':>9} {'speedup':>9} {'win/tie':>9} {'base bal':>8} {'new bal':>8} {'infeas':>5}"
    print("\n" + hdr); print("-" * len(hdr))
    for dep in sorted({r["depth"] for r in rows}):
        summary(f"depth {dep} (k={rows[[r['depth'] for r in rows].index(dep)]['k']})", [r for r in rows if r["depth"] == dep])
    print("-" * len(hdr))
    for name in names:
        summary(name, [r for r in rows if r["graph"] == name])
    print("-" * len(hdr))
    summary("ALL", rows)
    print("cut ratio = geomean new/kaffpa (<1: new cuts less); win/tie = items where new has the smaller / equal cut;")
    print("bal = mean max block weight / allowed (<=1 feasible); infeas = items where new exceeds the allowed weight")
    print(f"\nworst {a.worst} items by cut ratio:")
    for r in sorted((r for r in rows if r["b_cut"] > 0), key=lambda r: -r["n_cut"] / r["b_cut"])[:a.worst]:
        print(f"  {r['graph']:14} {r['item']:22} n={r['n']:5} k={r['k']} cut {r['b_cut']:.0f} -> {r['n_cut']:.0f} ({r['n_cut']/r['b_cut']:.3f}), bal {r['n_bal']:.3f}")


if __name__ == "__main__":
    main()
