#!/usr/bin/env python3
"""Writes a smaller Poisson test instance for Waver-HeiProMap and the wse-spmv evaluator into data/small<side>/:

  input.graph      METIS graph of the 7-point stencil on a side^3 grid, row = i + side * (j + side * k)
  topology.graph   METIS graph of the cols x rows PE grid (PE = row * cols + col), 4-neighborhood
  base.wspvmap     any WSPVMAP1 of the matrix: its header carries the matrix fingerprint that mappings are packed with
  instance.json    sizes and paths, read by tools/refine_loop.py

The full instance of wse-spmv/run_all_mappers.py is 215 x 215 x 216 on 316 x 316 PEs (about 100 rows per PE). Here
cols = rows = round(sqrt(side^3 / 100)), e.g. side 100 -> 1M rows on 100 x 100 PEs. The files that spmv-map-eval writes
to /tmp/heipromap_*.graph are not touched.

  python3 tools/make_small_instance.py [--side 100]
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

WSE = Path(os.path.expanduser("~/Projects/wse-spmv"))
HEI = Path(__file__).resolve().parent.parent


def write_stencil_graph(path, side):
    n = side ** 3
    m = 3 * side * side * (side - 1)
    s2 = side * side
    with open(path, "w") as f:
        f.write(f"{n} {m}\n")
        for u in range(n):
            i, j, k = u % side, (u // side) % side, u // s2
            nb = []
            if k > 0: nb.append(u - s2)
            if j > 0: nb.append(u - side)
            if i > 0: nb.append(u - 1)
            if i < side - 1: nb.append(u + 1)
            if j < side - 1: nb.append(u + side)
            if k < side - 1: nb.append(u + s2)
            f.write(" ".join(str(v + 1) for v in nb) + "\n")


def write_grid_graph(path, cols, rows):
    n = cols * rows
    m = rows * (cols - 1) + cols * (rows - 1)
    with open(path, "w") as f:
        f.write(f"{n} {m}\n")
        for p in range(n):
            x, y = p % cols, p // cols
            nb = []
            if y > 0: nb.append(p - cols)
            if x > 0: nb.append(p - 1)
            if x < cols - 1: nb.append(p + 1)
            if y < rows - 1: nb.append(p + cols)
            f.write(" ".join(str(v + 1) for v in nb) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--side", type=int, default=100)
    args = ap.parse_args()
    side = args.side
    n = side ** 3
    cols = rows = round((n / 100) ** 0.5)
    out = HEI / "data" / f"small{side}"
    out.mkdir(parents=True, exist_ok=True)

    write_stencil_graph(out / "input.graph", side)
    write_grid_graph(out / "topology.graph", cols, rows)

    os.chdir(WSE)
    sys.path.insert(0, str(WSE))
    argv, sys.argv = sys.argv, [sys.argv[0]]
    import run_all_mappers as ram
    sys.argv = argv
    template = list(ram.COMMAND_TEMPLATE)
    for flag, value in (("--nx", side), ("--ny", side), ("--nz", side), ("--cols", cols), ("--rows", rows)):
        template[template.index(flag) + 1] = str(value)
    base = out / "base.wspvmap"
    subprocess.run(ram.limited(template + ["--map", "kaminpar/recbisect", "--save-mapping", str(base)]),
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    json.dump({"nx": side, "ny": side, "nz": side, "cols": cols, "rows": rows, "rows_per_pe": 100,
               "graph": str(out / "input.graph"), "topology": str(out / "topology.graph"), "base": str(base)},
              open(out / "instance.json", "w"), indent=1)
    print(f"wrote {out}: {n} rows on {cols} x {rows} PEs")


if __name__ == "__main__":
    main()
