#!/usr/bin/env python3
"""Prototype: sender-group moves as negative cycles of an improvement graph.

The input is a state written by the multilevel solver with --dump-state PREFIX (graph, partition, grid, cost model,
weight limits of the cells of one level). The state before the scatter and chain refinements (pre_group) or after them
(post_group) is the start.

Improvement graph
  nodes   the blocks (PEs)
  edges   A -> C: the sender group S_AB (the vertices of A with an edge to a block B) moves to C, where C is B or a
          quotient neighbor of B. The weight is the cost change of that move alone (negative = it pays); of the groups
          to the same C only the best is kept.
  walks   a cycle A1 -> A2 -> ... -> A1 (every block gives one group and gets one) or a chain A1 -> ... -> An that ends
          in a block with room for what it gets. The weights of a walk do not add up exactly (the moves interact) and
          the weight limits are a condition on two neighboring edges, so a walk found with negative total weight is a
          candidate: its moves are applied one after the other with their exact cost change and undone if the total
          does not pay or a limit is violated.

  python3 tools/improvement_graph.py --prefix DIR/s --call 10 [--start pre_group] [--depth 4] [--branch 6]
"""
import argparse
import sys
import time
from collections import defaultdict


# ---------------------------------------------------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------------------------------------------------
def read_metis(path):
    """METIS format with vertex and edge weights (011). Returns vertex weights and the adjacency lists."""
    with open(path) as f:
        n, m, fmt = f.readline().split()
        n = int(n)
        assert fmt == "011", fmt
        vw, adj = [], []
        for _ in range(n):
            t = f.readline().split()
            vw.append(int(t[0]))
            adj.append([(int(t[i]) - 1, int(t[i + 1])) for i in range(1, len(t), 2)])
    return vw, adj


def read_meta(path):
    meta = {"cost": {}}
    with open(path) as f:
        for line in f:
            key, *values = line.split()
            if key in ("lmin", "lmax"):
                meta[key] = [int(v) for v in values]
            elif key.startswith("cost_"):
                meta["cost"][key[5:]] = float(values[0])
            elif key in ("n", "k", "W"):
                meta[key] = int(values[0])
            else:
                meta[key] = float(values[0])
    return meta


# ---------------------------------------------------------------------------------------------------------------------
# the partition with everything the moves need
# ---------------------------------------------------------------------------------------------------------------------
class State:
    def __init__(self, vw, adj, part, meta):
        self.vw, self.adj, self.part = vw, adj, list(part)
        self.k, self.W = meta["k"], meta["W"]
        self.a_end, self.b_end = meta["alpha_end"], meta["beta_end"]
        self.a_link, self.b_link = meta["alpha_link"], meta["beta_link"]
        self.lmin, self.lmax = meta["lmin"], meta["lmax"]
        self.bw = [0] * self.k
        for u, b in enumerate(self.part):
            self.bw[b] += vw[u]
        # The solver repairs the lower limit only after the refinements, so some blocks may be outside their band at
        # the start. A block may not get worse than it is: its band is widened to its start weight.
        self.lo = [min(self.lmin[b], self.bw[b]) for b in range(self.k)]
        self.hi = [max(self.lmax[b], self.bw[b]) for b in range(self.k)]
        # nbw[u]: block -> edge weight from u into that block
        self.nbw = []
        for u in range(len(vw)):
            d = defaultdict(int)
            for v, w in adj[u]:
                d[self.part[v]] += w
            self.nbw.append(d)
        # quotient rows Q[a][b]: cut edge weight between the blocks (both directions kept)
        self.Q = [defaultdict(int) for _ in range(self.k)]
        for u in range(len(vw)):
            a = self.part[u]
            for b, w in self.nbw[u].items():
                if b != a:
                    self.Q[a][b] += w
        # boundary vertices of every block
        self.bnd = [set() for _ in range(self.k)]
        for u in range(len(vw)):
            if self.is_boundary(u):
                self.bnd[self.part[u]].add(u)

    def is_boundary(self, u):
        a = self.part[u]
        return any(b != a for b in self.nbw[u])

    def dist(self, a, c):
        W = self.W
        return abs(a % W - c % W) + abs(a // W - c // W)

    def pair_cost(self, a, c, w):
        """Cost of the (directed) message a -> c with edge weight w; both directions are counted by the caller."""
        if w <= 0:
            return 0.0
        h = self.dist(a, c)
        return self.a_end + self.a_link * h + (self.b_end + self.b_link * h) * w

    def total_cost(self):
        total = 0.0
        for a in range(self.k):
            for c, w in self.Q[a].items():
                total += self.pair_cost(a, c, w)
        return total

    # -- groups ------------------------------------------------------------------------------------------------------
    def groups(self, a):
        """block b -> vertices of a with an edge to b."""
        g = defaultdict(list)
        for u in self.bnd[a]:
            for b in self.nbw[u]:
                if b != a:
                    g[b].append(u)
        return g

    def measure(self, a, group):
        """(weight, ext, own): the weight of the group, its edge weight to every other block, and to a outside it."""
        inside = set(group)
        ext, own, weight = defaultdict(int), 0, 0
        for u in group:
            weight += self.vw[u]
            for v, w in self.adj[u]:
                d = self.part[v]
                if d == a:
                    if v not in inside:
                        own += w
                else:
                    ext[d] += w
        return weight, ext, own

    def delta(self, a, c, ext, own):
        """Cost change of moving a group of a with the given ext/own to c (negative = better)."""
        Q = self.Q
        change = 0.0
        for d, w in ext.items():
            if d == c:
                continue
            # A-D loses w, C-D gains w
            change += 2.0 * (self.pair_cost(a, d, Q[a][d] - w) - self.pair_cost(a, d, Q[a][d]))
            change += 2.0 * (self.pair_cost(c, d, Q[c][d] + w) - self.pair_cost(c, d, Q[c][d]))
        w_ac = Q[a][c]
        new_ac = w_ac - ext.get(c, 0) + own
        change += 2.0 * (self.pair_cost(a, c, new_ac) - self.pair_cost(a, c, w_ac))
        return change

    # -- moves -------------------------------------------------------------------------------------------------------
    def move_vertex(self, v, a, c):
        self.part[v] = c
        self.bw[a] -= self.vw[v]
        self.bw[c] += self.vw[v]
        Q = self.Q
        touched = [v]
        for u, w in self.adj[v]:
            d = self.part[u]
            if d != a:                      # the edge was cut between a and d
                Q[a][d] -= w
                Q[d][a] -= w
                if Q[a][d] == 0:
                    del Q[a][d]
                    del Q[d][a]
            if d != c:                      # now it is cut between c and d
                Q[c][d] += w
                Q[d][c] += w
            nb = self.nbw[u]
            nb[a] -= w
            if nb[a] == 0:
                del nb[a]
            nb[c] += w
            touched.append(u)
        for x in touched:
            b = self.part[x]
            if self.is_boundary(x):
                self.bnd[b].add(x)
            else:
                self.bnd[b].discard(x)
        # v changed its block: it may sit in the old block's set
        self.bnd[a].discard(v)

    def move_group(self, group, a, c):
        for v in group:
            self.move_vertex(v, a, c)

    def check(self, blocks=None):
        """The blocks (default: all) outside their band; returns them."""
        return [b for b in (range(self.k) if blocks is None else blocks) if not (self.lo[b] <= self.bw[b] <= self.hi[b])]


# ---------------------------------------------------------------------------------------------------------------------
# the improvement graph and its search
# ---------------------------------------------------------------------------------------------------------------------
class Edge:
    __slots__ = ("a", "c", "delta", "weight", "group")

    def __init__(self, a, c, delta, weight, group):
        self.a, self.c, self.delta, self.weight, self.group = a, c, delta, weight, group


def out_edges(s, a, targets, keep_all):
    """The edges for the sender groups of block a. The targets of the group towards b are b and the `targets` quotient
    neighbors of b with the most edge weight. Of the groups to the same C only the best is kept, unless keep_all."""
    best, every = {}, []
    for b, group in s.groups(a).items():
        weight, ext, own = s.measure(a, group)
        near = sorted((x for x in s.Q[b] if x != a), key=lambda x: -s.Q[b][x])[:targets]
        for c in [b] + near:
            if c == a:
                continue
            d = s.delta(a, c, ext, own)
            e = Edge(a, c, d, weight, group)
            if keep_all:
                every.append(e)
            elif c not in best or d < best[c].delta:
                best[c] = e
    return sorted(every if keep_all else best.values(), key=lambda e: e.delta)


def search(s, edges, depth, branch, eps, min_edge):
    """All walks with a negative total weight (cycles and chains), by depth limited depth first search."""
    found = []
    lmin, lmax, bw = s.lo, s.hi, s.bw

    def extend(path, nodes, start, in_weight, total):
        cur = nodes[-1]
        taken = 0
        for e in edges[cur]:
            g = e.weight
            # limits of cur after it received the previous group and gave this one
            if not (lmin[cur] <= bw[cur] + in_weight - g <= lmax[cur]):
                continue
            taken += 1
            if taken > branch:
                break
            new_total = total + e.delta
            t = e.c
            if t == start:
                first = path[0].weight
                if lmin[start] <= bw[start] - first + g <= lmax[start] and new_total < -eps:
                    found.append((new_total, "cycle", path + [e]))
                continue
            if t in nodes:
                continue
            if bw[t] + g <= lmax[t] and new_total < -eps:
                found.append((new_total, "chain", path + [e]))
            remaining = depth - len(path) - 1
            if remaining > 0 and new_total + remaining * min_edge < -eps:
                extend(path + [e], nodes + [t], start, g, new_total)

    for a in range(s.k):
        taken = 0
        for e in edges[a]:
            if bw[a] - e.weight < lmin[a]:
                continue
            taken += 1
            if taken > branch:
                break
            t = e.c
            if bw[t] + e.weight <= lmax[t] and e.delta < -eps:
                found.append((e.delta, "chain", [e]))
            if depth > 1 and e.delta + (depth - 1) * min_edge < -eps:
                extend([e], [a, t], a, e.weight, e.delta)
    return found


def try_walk(s, walk, eps):
    """Applies the moves of a walk one after the other with their exact cost change. Keeps them if all limits hold and
    the total pays, otherwise undoes them. Returns the exact total change, or None if it was undone."""
    done, total = [], 0.0
    for e in walk:
        group = [u for u in e.group if s.part[u] == e.a]
        if not group:
            break
        weight, ext, own = s.measure(e.a, group)
        total += s.delta(e.a, e.c, ext, own)
        s.move_group(group, e.a, e.c)
        done.append((group, e.a, e.c))
    blocks = {e.a for e in walk} | {e.c for e in walk}
    ok = len(done) == len(walk) and not s.check(blocks) and total < -eps
    if ok:
        return total
    for group, a, c in reversed(done):
        s.move_group(group, c, a)
    return None


def run(s, depth, branch, targets, keep_all, eps, max_passes):
    t0 = time.time()
    edges = [out_edges(s, a, targets, keep_all) for a in range(s.k)]
    print(f"  improvement graph: {sum(len(e) for e in edges)} edges, built in {time.time() - t0:.1f}s")
    bad = set()
    kinds = defaultdict(int)
    lengths = defaultdict(int)
    cost = s.total_cost()
    for p in range(max_passes):
        t1 = time.time()
        min_edge = min((e[0].delta for e in edges if e), default=0.0)
        found = search(s, edges, depth, branch, eps, min_edge)
        found.sort(key=lambda x: x[0])
        t2 = time.time()
        dirty, kept, tried = set(), 0, 0
        for total, kind, walk in found:
            key = tuple((e.a, e.c, len(e.group)) for e in walk)
            blocks = {e.a for e in walk} | {e.c for e in walk}
            if key in bad or blocks & dirty:
                continue
            tried += 1
            exact = try_walk(s, walk, eps)
            if exact is None:
                bad.add(key)
                continue
            kept += 1
            kinds[kind] += 1
            lengths[len(walk)] += 1
            cost += exact
            for b in blocks:
                dirty.add(b)
                dirty.update(s.Q[b].keys())
        # the edges of the blocks around the changed ones are out of date
        around = set(dirty)
        for b in dirty:
            around.update(s.Q[b].keys())
        for b in around:
            edges[b] = out_edges(s, b, targets, keep_all)
        print(f"  pass {p}: {len(found)} walks found ({t2 - t1:.1f}s), {tried} tried, {kept} kept, cost {cost:.6g}, "
              f"{time.time() - t2:.1f}s to apply and update")
        if kept == 0:
            break
    print(f"  kept: {dict(kinds)}, by number of moves: {dict(sorted(lengths.items()))}")
    return cost


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prefix", required=True, help="the --dump-state prefix of the solver")
    ap.add_argument("--call", type=int, required=True, help="number of the refinement call (level) to use")
    ap.add_argument("--start", default="pre_group", choices=["pre_group", "post_group"])
    ap.add_argument("--depth", type=int, default=4, help="maximum number of moves of a walk")
    ap.add_argument("--branch", type=int, default=6, help="best outgoing edges followed per block")
    ap.add_argument("--targets", type=int, default=4, help="quotient neighbors of the receiver tried as targets")
    ap.add_argument("--keep-all", action="store_true", help="keep every group of a block towards every target as an edge")
    ap.add_argument("--passes", type=int, default=20)
    ap.add_argument("--eps", type=float, default=1e-6)
    args = ap.parse_args()

    base = f"{args.prefix}_{args.call}"
    meta = read_meta(base + ".meta")
    vw, adj = read_metis(base + ".graph")
    read_part = lambda tag: [int(x) for x in open(f"{base}_{tag}.part").read().split()]
    start = State(vw, adj, read_part(args.start), meta)
    print(f"level {args.call}: n={meta['n']} k={meta['k']} grid width {meta['W']}")
    cost0 = start.total_cost()
    print(f"cost of {args.start}: {cost0:.6g} (solver: {meta['cost'][args.start]:.6g})")
    if abs(cost0 - meta["cost"][args.start]) > 1e-6 * cost0:
        sys.exit("the cost does not match the solver's: the cost model or the input is wrong")
    for tag in ("pre_group", "post_group"):
        print(f"  solver {tag}: {meta['cost'][tag]:.6g}")
    print(f"blocks outside their band at the start: {len(start.check())} (they may not get worse)")

    final = run(start, args.depth, args.branch, args.targets, args.keep_all, args.eps, args.passes)
    print(f"final cost {final:.6g}, recomputed {start.total_cost():.6g}; balance violations {len(start.check())}")
    pre, post = meta["cost"]["pre_group"], meta["cost"]["post_group"]
    print(f"gain over {args.start}: {100 * (cost0 - final) / cost0:.2f}%   "
          f"(scatter + chain + scatter of the solver: {100 * (pre - post) / pre:.2f}% from pre_group)")


if __name__ == "__main__":
    main()
