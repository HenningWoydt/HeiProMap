#!/usr/bin/env python3
"""Serial: flow-bisection (c=8, seed 0) with --flow-bisection-prune false vs true. QAPs must be identical."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench_flow_bisection as b

a = argparse.Namespace(hierarchy="4:8:6", distance="1:10:100", imbalance="0.03", threads=1, initial_c=8,
                       memlimit="100G", no_memlimit=False, out=b.ROOT / "data" / "bench_flow_bisection")
(a.out / "mappings").mkdir(parents=True, exist_ok=True)
same, t_off, t_on, ip_off, ip_on = 0, [], [], [], []
print(f"{'config':<7}{'graph':<14}{'QAP off':>10}{'QAP on':>10}{'t off':>7}{'t on':>7}{'init ms off':>12}{'init ms on':>11}")
for c in ["fast", "eco", "strong"]:
    for g in b.DEFAULT_GRAPHS:
        res = {}
        for mode in ["false", "true"]:
            a.extra = ["--flow-bisection-prune", mode]
            res[mode] = b.run_one((g, c, "flow", 0, a))
        off, on = res["false"], res["true"]
        same += off["final_qap"] == on["final_qap"]
        t_off.append(off["time_s"]); t_on.append(on["time_s"])
        ip_off.append(off["init_part_ms"]); ip_on.append(on["init_part_ms"])
        print(f"{c:<7}{g:<14}{off['final_qap']:>10.0f}{on['final_qap']:>10.0f}{off['time_s']:>7.2f}{on['time_s']:>7.2f}"
              f"{off['init_part_ms']:>12.0f}{on['init_part_ms']:>11.0f}", flush=True)
print(f"\nidentical final QAP: {same}/15; geomean total time on/off = {b.geomean([n/o for n,o in zip(t_on,t_off)]):.3f}; "
      f"init-partitioning on/off = {b.geomean([n/o for n,o in zip(ip_on,ip_off)]):.3f}")
