#!/usr/bin/env python3
"""Waver on a 2D grid graph whose best mapping is known.

An L x L grid graph (5-point stencil, METIS format) on a 100 x 100 PE grid: the best mapping cuts it into tiles of
(L/100) x (L/100) vertices, one per PE, so every PE has 4 neighbors, all one hop away. This script writes the graph,
runs Waver-HeiProMap on it (extra arguments are passed on) and compares the result with the tiles:

  python3 grid_test.py [--side 1000] [--no-run] [waver args...]

Reported: the cut and the message count against the tiling's, the hops, and how many vertices sit in the cell the best
of the 8 symmetries of the tiling puts them in. Per PE: the connected pieces of its vertices and the area of their
bounding box against the tile area.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = 100  # PEs per side


def write_graph(side, path):
    n = side * side
    idx = np.arange(n).reshape(side, side)
    right = np.stack([idx[:, :-1].ravel(), idx[:, 1:].ravel()], axis=1)
    down = np.stack([idx[:-1, :].ravel(), idx[1:, :].ravel()], axis=1)
    edges = np.concatenate([right, down])
    both = np.concatenate([edges, edges[:, ::-1]])
    order = np.lexsort((both[:, 1], both[:, 0]))
    both = both[order]
    counts = np.bincount(both[:, 0], minlength=n)
    ptr = np.concatenate([[0], np.cumsum(counts)])
    with open(path, "w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u in range(n):
            f.write(" ".join(map(str, both[ptr[u]:ptr[u + 1], 1] + 1)) + "\n")
    return edges


def write_topology(path):
    """The GRID x GRID PE grid (4-neighborhood, unit weights), METIS format."""
    with open(path, "w") as f:
        f.write(f"{GRID * GRID} {2 * GRID * (GRID - 1)}\n")
        for r in range(GRID):
            for c in range(GRID):
                nb = [(r + dr) * GRID + c + dc + 1
                      for dr, dc in ((-1, 0), (0, -1), (0, 1), (1, 0))
                      if 0 <= r + dr < GRID and 0 <= c + dc < GRID]
                f.write(" ".join(map(str, nb)) + "\n")


def metrics(side, edges, part):
    n = side * side
    tile = side // GRID
    pe_r, pe_c = part // GRID, part % GRID
    r, c = np.divmod(np.arange(n), side)
    ideal_r, ideal_c = r // tile, c // tile
    best = (-1, None)
    for transpose in (False, True):
        for flip_r in (False, True):
            for flip_c in (False, True):
                a, b = (ideal_c, ideal_r) if transpose else (ideal_r, ideal_c)
                if flip_r:
                    a = GRID - 1 - a
                if flip_c:
                    b = GRID - 1 - b
                agree = int(np.sum((a == pe_r) & (b == pe_c)))
                if agree > best[0]:
                    best = (agree, (transpose, flip_r, flip_c))
    u, v = edges[:, 0], edges[:, 1]
    cut = part[u] != part[v]
    pairs = np.unique(np.sort(np.stack([part[u][cut], part[v][cut]], axis=1), axis=1), axis=0)
    hops = np.abs(pairs[:, 0] // GRID - pairs[:, 1] // GRID) + np.abs(pairs[:, 0] % GRID - pairs[:, 1] % GRID)
    ideal_cut = 2 * side * (GRID - 1)
    ideal_pairs = 2 * GRID * (GRID - 1)

    # connected pieces and bounding boxes of the vertices of every PE
    same = ~cut
    comp_count, comp = connected_components(
        coo_matrix((np.ones(int(same.sum())), (u[same], v[same])), shape=(n, n)), directed=False)
    comp_pe = np.zeros(comp_count, dtype=np.int64)
    comp_pe[comp] = part
    pieces = np.bincount(comp_pe, minlength=GRID * GRID)
    areas = np.zeros(GRID * GRID)
    for pe in range(GRID * GRID):
        members = np.nonzero(part == pe)[0]
        if len(members):
            rr, cc = r[members], c[members]
            areas[pe] = (rr.max() - rr.min() + 1) * (cc.max() - cc.min() + 1)
    print(f"vertices in their tile cell (best symmetry {best[1]}): {best[0] / n:.1%}")
    print(f"cut edges {int(cut.sum())} (tiling {ideal_cut}, {cut.sum() / ideal_cut:.2f}x)")
    print(f"messages {len(pairs)} (tiling {ideal_pairs}, {len(pairs) / ideal_pairs:.2f}x), "
          f"mean hops {hops.mean():.2f}, hops>1 {np.mean(hops > 1):.1%}, longest {hops.max()}")
    print(f"pieces per PE: mean {pieces.mean():.2f}, max {pieces.max()}, PEs with 1 piece {np.mean(pieces == 1):.1%}")
    print(f"bounding box area / tile area: mean {areas.mean() / tile ** 2:.2f}, p99 {np.percentile(areas, 99) / tile ** 2:.2f}")


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--side", type=int, default=1000)
    ap.add_argument("--no-run", action="store_true")
    known, extra = ap.parse_known_args()
    graph = "/tmp/grid2d.graph"
    topology = "/tmp/grid2d_topology.graph"
    mapping = "/tmp/grid2d_mapping.txt"
    edges = write_graph(known.side, graph)
    if not known.no_run:
        write_topology(topology)
        env = dict(os.environ, GRAPH=graph, TOPOLOGY=topology, OUT=mapping)
        subprocess.run([os.path.join(HERE, "run_waver.sh"), "--level-stats", "1", *extra], env=env, check=True)
    part = np.loadtxt(mapping, dtype=np.int64)
    metrics(known.side, edges, part)


if __name__ == "__main__":
    sys.exit(main())
