#!/usr/bin/env python3
"""Which refinement of Waver-HeiProMap has the most impact?

Three experiments, all on the same graph and topology, scored by the makespan of the timed SpMV model
(Waver-Schedule) and by the numbers Waver-HeiProMap reports itself:

  leave-one-out  FULL (every refinement on) against FULL with one refinement switched off. The increase of the
                 makespan is what that refinement contributes in the presence of the others.
  add-one-in     BARE (every refinement off) against BARE with one refinement switched on. This is what the
                 refinement can do by itself. A refinement that is strong alone but weak in leave-one-out is
                 covered by another one.
  on-top         FULL with one more refinement or option switched on that is not part of FULL (the expensive or
                 experimental ones). The decrease of the makespan is what it adds to everything else.

  python3 tools/ablation.py --mode spectral|multilevel [--seeds 3] [--jobs 4] [--dose] [--out DIR]

Every run uses --config strong as its base, so the settings of a refinement that is on are those of the preset
(for example swap radius 2, annealing 5000 iterations, 20 corner rounds); the refinements are switched with the
command line options. The report is about quality; the times of the runs (their own phase timers) are in runs.json,
and in the report with --show-time. With --jobs above 1 the runs share the memory bandwidth, so the times are
somewhat higher than when a run is alone (use --jobs 1 for those).
"""
import argparse
import concurrent.futures
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def switch(option, on):
    """A 0/1 option."""
    return ([option, "0"], [option, str(on)])


# ---- the refinements of the spectral mode -------------------------------------------------------------------------
# name -> (arguments that switch it off, arguments that switch it on with the settings of --config strong)
SPECTRAL = {
    "components": {
        "message": switch("--message-refinement", 1),
        "swap": (["--swap-radius", "0"], ["--swap-radius", "2"]),
        "long-message": (["--long-message-padding", "0"], ["--long-message-padding", "1"]),
        "annealing": (["--annealing-iterations", "0"], ["--annealing-iterations", "5000"]),
        "corner": (["--corner-rounds", "0"], ["--corner-rounds", "20"]),
        "vcycle": switch("--vcycle", 1),
    },
    # the refinements of the full configuration (the V-cycle is expensive and off in every preset but super-strong)
    "full": ["message", "swap", "long-message", "annealing", "corner"],
    # options switched on top of FULL: name -> arguments
    "extras": {},
    # structural options of FULL, switched off: name -> arguments
    "structure": {},
    # parameter sweeps around the full configuration: name -> list of (label, arguments)
    "dose": {
        "message iterations": [(str(i), ["--message-refinement-iterations", str(i)]) for i in (1, 5, 20, 50)],
        "swap radius": [(str(r), ["--swap-radius", str(r)]) for r in (1, 2, 3, 4)],
        "long-message padding": [(str(p), ["--long-message-padding", str(p)]) for p in (1, 2, 3)],
        "annealing iterations": [(str(i), ["--annealing-iterations", str(i)]) for i in (1000, 5000, 20000)],
        "corner rounds": [(str(r), ["--corner-rounds", str(r)]) for r in (5, 20, 40)],
    },
}

# ---- the refinements of the multilevel mode -----------------------------------------------------------------------
MULTILEVEL = {
    "components": {
        "message": switch("--message-refinement", 1),
        "scatter": switch("--scatter-refinement", 1),
        "chain": switch("--chain-refinement", 1),
        "group-split": (["--group-split-parts", "0"],
                        ["--group-split-parts", "3", "--group-split-depth", "2", "--group-split-rounds", "2"]),
        "smooth": switch("--smooth-refinement", 1),
        "cycle": switch("--cycle-refinement", 1),
        "busy": switch("--busy-refinement", 1),
        "crowd": switch("--crowd-refinement", 1),
        "shift": switch("--shift-refinement", 1),
        "swap": (["--swap-radius", "0"], ["--swap-radius", "2"]),
        "long-message": (["--long-message-padding", "0"], ["--long-message-padding", "1"]),
        "annealing": (["--annealing-iterations", "0"], ["--annealing-iterations", "5000"]),
        "corner": (["--corner-rounds", "0"], ["--corner-rounds", "20"]),
        # not part of FULL: expensive or experimental
        "vcycle": switch("--vcycle", 1),
        "smash": (["--smash-rounds", "0"], ["--smash-rounds", "2"]),
        "window": (["--window-rounds", "0"], ["--window-rounds", "2"]),
        "relax": (["--relax-rounds", "0"], ["--relax-rounds", "2"]),
    },
    "full": ["message", "scatter", "chain", "group-split", "smooth", "cycle", "busy", "crowd", "shift", "swap",
             "long-message", "annealing", "corner"],
    "extras": {
        "terminal-propagation": ["--terminal-propagation", "1"],
        "global-split": ["--global-split", "1"],
        "diffusion-rebalance": ["--diffusion-rebalance", "1"],
        "split-refine-passes": ["--split-refine-passes", "3"],
        "rebalance-nearest": ["--rebalance-nearest", "1"],
        "coarse-slack": ["--coarse-slack", "1"],
    },
    "structure": {
        "split-kway (k-way local split)": ["--split-kway", "0"],
    },
    "dose": {},
}

MODES = {"spectral": SPECTRAL, "multilevel": MULTILEVEL}

# Waver-HeiProMap summary rows: key -> regular expression
SUMMARY = {
    "cost": r"^Total message cost\s*:\s*([0-9.e+]+)",
    "messages": r"^Message count \(\|E_pi\|\)\s*:\s*([0-9.e+]+)",
    "hop_volume": r"^Hop volume \(sum h\*S\)\s*:\s*([0-9.e+]+)",
    "max_link": r"^Max link load\s*:\s*([0-9.e+]+)",
    "max_endpoint": r"^Max PE endpoint load\s*:\s*([0-9.e+]+)",
    "time_total": r"^Total time \(s\)\s*:\s*([0-9.e+]+)",
    "ms_prepartition": r"^Prepartition\s*:\s*([0-9.e+]+)",
    "ms_refinement": r"^Refinement\s*:\s*([0-9.e+]+)",
    "ms_message": r"^  Msg Refinement\s*:\s*([0-9.e+]+)",
    "ms_swap": r"^  Swap\s*:\s*([0-9.e+]+)",
    "ms_long_message": r"^  Long Msg\s*:\s*([0-9.e+]+)",
    "ms_group": r"^  Scatter/Chain/Split\s*:\s*([0-9.e+]+)",
    "ms_busy": r"^  Busy\s*:\s*([0-9.e+]+)",
    "ms_annealing": r"^  Sim\. Annealing\s*:\s*([0-9.e+]+)",
    "ms_corner": r"^  Corner \(final\)\s*:\s*([0-9.e+]+)",
}
# the time of the refinements: the part of the run that the switches change
REFINEMENT_MS = ("ms_refinement", "ms_corner")


def merge_args(*groups):
    """Concatenates option lists; an option that occurs twice keeps the last value."""
    values = {}
    for group in groups:
        for i in range(0, len(group), 2):
            values[group[i]] = group[i + 1]
    out = []
    for key, value in values.items():
        out += [key, value]
    return out


def variants(spec, dose_wanted):
    """The list of (experiment, name, extra arguments)."""
    components, full = spec["components"], spec["full"]
    off = {c: components[c][0] for c in components}
    on = {c: components[c][1] for c in components}
    full_args = merge_args(*[on[c] for c in full])
    bare_args = merge_args(*[off[c] for c in components])
    runs = [("reference", "FULL", full_args), ("reference", "BARE", bare_args)]
    for c in full:
        runs.append(("leave-one-out", c, merge_args(full_args, off[c])))
    for name, arguments in spec["structure"].items():
        runs.append(("leave-one-out", name, merge_args(full_args, arguments)))
    for c in components:
        runs.append(("add-one-in", c, merge_args(bare_args, on[c])))
    for c in components:
        if c not in full:
            runs.append(("on-top", c, merge_args(full_args, on[c])))
    for name, arguments in spec["extras"].items():
        runs.append(("on-top", name, merge_args(full_args, arguments)))
    if dose_wanted:
        for parameter, values in spec["dose"].items():
            for label, arguments in values:
                runs.append(("dose", f"{parameter} = {label}", merge_args(full_args, arguments)))
    return runs


def check_options(binary, runs, extra):
    """Waver-HeiProMap ignores an option it does not know, which would silently make a run a copy of another one."""
    help_text = subprocess.run([binary, "--help"], capture_output=True, text=True).stdout
    known = set(re.findall(r"^\s+(--[a-z0-9-]+)", help_text, re.M))
    unknown = sorted({a for _, _, args in runs for a in args if a.startswith("--") and a not in known} |
                     {a for a in extra if a.startswith("--") and a not in known})
    if unknown:
        sys.exit("options the binary does not know (rebuild it?): " + ", ".join(unknown))


def parse_summary(text):
    result = {}
    for line in ANSI.sub("", text).splitlines():
        for key, pattern in SUMMARY.items():
            if key not in result:
                match = re.match(pattern, line)
                if match:
                    result[key] = float(match.group(1))
    return result


def wrap_memory_limit(command, limit_gb):
    """Runs the command in a systemd scope with a memory limit if that is available."""
    if limit_gb <= 0 or shutil.which("systemd-run") is None:
        return command
    return ["systemd-run", "--user", "--scope", "--quiet", "-p", f"MemoryMax={limit_gb}G"] + command


def run_one(job, args):
    experiment, name, extra, seed = job
    with tempfile.TemporaryDirectory(prefix="ablation_") as tmp:
        mapping = os.path.join(tmp, "mapping.txt")
        command = [args.binary, "--graph", args.graph, "--topology", args.topology, "--mapping", mapping,
                   "--imbalance", str(args.imbalance), "--config", "strong", "--mode", args.mode,
                   "--seed", str(seed), "--threads", str(args.threads)] + merge_args(args.extra, extra)
        started = time.time()
        solver = subprocess.run(wrap_memory_limit(command, args.memory_gb), capture_output=True, text=True)
        wall = time.time() - started
        record = {"experiment": experiment, "name": name, "seed": seed, "args": extra, "wall": wall}
        if solver.returncode != 0 or not os.path.exists(mapping):
            record["error"] = f"solver exit {solver.returncode}: {(solver.stdout + solver.stderr)[-300:]}"
            return record
        record.update(parse_summary(solver.stdout))
        schedule = [args.schedule, "--graph", args.graph, "--topology", args.topology, "--mapping", mapping,
                    "--chain", "0", "--imbalance", str(args.imbalance)]
        if args.costs:
            schedule += ["--costs", args.costs]
        scored = subprocess.run(wrap_memory_limit(schedule, args.memory_gb), capture_output=True, text=True)
        match = re.search(r"^makespan\s+(\d+)", scored.stdout, re.M)
        if scored.returncode != 0 or not match:
            record["error"] = f"schedule exit {scored.returncode}: {(scored.stdout + scored.stderr)[-300:]}"
            return record
        record["makespan"] = float(match.group(1))
        match = re.search(r"max engine load\s+(\d+)", scored.stdout)
        if match:
            record["max_engine_load"] = float(match.group(1))
        record["ms_refinements"] = sum(record.get(key, 0.0) for key in REFINEMENT_MS)
        return record


def mean_of(records, key):
    values = [r[key] for r in records if key in r]
    return statistics.mean(values) if values else float("nan")


def std_of(records, key):
    values = [r[key] for r in records if key in r]
    return statistics.stdev(values) if len(values) > 1 else 0.0


def percent(change, base):
    return 100.0 * change / base if base else float("nan")


def report(records, order, out, show_time):
    groups = {}
    for r in records:
        if "error" not in r:
            groups.setdefault((r["experiment"], r["name"]), []).append(r)
    failed = [r for r in records if "error" in r]
    lines = []

    def row(*cells):
        lines.append("| " + " | ".join(cells) + " |")

    def get(experiment, name):
        return groups.get((experiment, name), [])

    def significant(records, reference, key):
        """The difference of the means is larger than twice the larger standard deviation (needs two seeds)."""
        if min(len([r for r in records if key in r]), len([r for r in reference if key in r])) < 2:
            return False
        diff = abs(mean_of(records, key) - mean_of(reference, key))
        return diff > 2.0 * max(std_of(records, key), std_of(reference, key))

    full, bare = get("reference", "FULL"), get("reference", "BARE")
    seeds = max((len(v) for v in groups.values()), default=0)
    lines.append("# Refinement ablation\n")
    lines.append(f"Means over {seeds} seeds. Makespan is that of the timed SpMV model (Waver-Schedule), cost the message "
                 "cost of the cost model, link and engine the largest link load and the largest engine load of a PE. "
                 "A positive impact is an improvement (a lower makespan, cost or load). `*` marks a difference of the "
                 "makespan larger than twice the larger standard deviation of the two configurations: the others are "
                 "within the noise of the seeds (with one seed nothing is marked).\n")
    lines.append("## References\n")
    header = ["Configuration", "Makespan", "Cost", "Max link", "Max endpoint", "Max engine"]
    if show_time:
        header += ["Refinement time (s)", "Total time (s)"]
    row(*header)
    row(*(["---"] + ["---:"] * (len(header) - 1)))
    for name, records in (("FULL", full), ("BARE", bare)):
        if records:
            cells = [f"`{name}`", f"{mean_of(records, 'makespan'):.0f} ±{std_of(records, 'makespan'):.0f}",
                     f"{mean_of(records, 'cost'):.4g}", f"{mean_of(records, 'max_link'):.0f}",
                     f"{mean_of(records, 'max_endpoint'):.0f}", f"{mean_of(records, 'max_engine_load'):.0f}"]
            if show_time:
                cells += [f"{mean_of(records, 'ms_refinements') / 1000:.1f}", f"{mean_of(records, 'time_total'):.1f}"]
            row(*cells)
    lines.append("")

    def table(title, experiment, reference, sign, caption):
        entries = []
        for name in order.get(experiment, []):
            records = get(experiment, name)
            if records and reference:
                entries.append((sign * (mean_of(reference, "makespan") - mean_of(records, "makespan")), name, records))
        if not entries:
            return
        # the improvement of the makespan, best first: sign +1 compares FULL with the variant, -1 the variant with BARE
        entries.sort(key=lambda e: -e[0])
        lines.append(f"## {title}\n")
        lines.append(caption + "\n")
        header = ["Refinement", "Makespan", "Impact on makespan", "on cost", "on max link", "on max engine"]
        if show_time:
            header += ["Refinement time (s)"]
        row(*header)
        row(*(["---"] + ["---:"] * (len(header) - 1)))
        for gain, name, records in entries:
            def delta(key):
                d = sign * (mean_of(reference, key) - mean_of(records, key))
                return f"{percent(d, mean_of(reference, key)):+.1f}%"
            mark = " *" if significant(records, reference, "makespan") else ""
            cells = [f"`{name}`", f"{mean_of(records, 'makespan'):.0f} ±{std_of(records, 'makespan'):.0f}",
                     f"{gain:+.0f} ({percent(gain, mean_of(reference, 'makespan')):+.1f}%){mark}",
                     delta("cost"), delta("max_link"), delta("max_engine_load")]
            if show_time:
                cells += [f"{mean_of(records, 'ms_refinements') / 1000:.1f}"]
            row(*cells)
        lines.append("")

    # leave-one-out: the variant is the worse one, so the refinement's contribution is variant - FULL
    table("Leave-one-out", "leave-one-out", full, -1,
          "The refinement is switched off in `FULL`. The impact is what is lost without it: a positive number means "
          "that the refinement helps.")
    table("Add-one-in", "add-one-in", bare, +1,
          "The refinement is switched on in `BARE`, where everything else is off. The impact is what it gains by itself.")
    table("On top of FULL", "on-top", full, +1,
          "The refinement or option is switched on in `FULL` in addition. The impact is what it adds to everything else.")
    if order.get("dose"):
        lines.append("## Parameter sweeps around FULL\n")
        header = ["Setting", "Makespan", "Impact on makespan", "on cost", "on max link", "on max engine"]
        row(*header)
        row(*(["---"] + ["---:"] * (len(header) - 1)))
        for name in order["dose"]:
            records = get("dose", name)
            if not records or not full:
                continue
            def delta(key):
                return f"{percent(mean_of(full, key) - mean_of(records, key), mean_of(full, key)):+.1f}%"
            gain = mean_of(full, "makespan") - mean_of(records, "makespan")
            mark = " *" if significant(records, full, "makespan") else ""
            row(f"`{name}`", f"{mean_of(records, 'makespan'):.0f} ±{std_of(records, 'makespan'):.0f}",
                f"{gain:+.0f} ({percent(gain, mean_of(full, 'makespan')):+.1f}%){mark}", delta("cost"),
                delta("max_link"), delta("max_engine_load"))
        lines.append("")
    if failed:
        lines.append("## Failed runs\n")
        for r in failed:
            lines.append(f"- {r['experiment']} / {r['name']} / seed {r['seed']}: {r['error']}")
        lines.append("")
    text = "\n".join(lines)
    with open(out, "w") as f:
        f.write(text)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", default="spectral", choices=sorted(MODES), help="solver mode")
    parser.add_argument("--graph", default=os.path.join(ROOT, "data", "heipromap_input.graph"))
    parser.add_argument("--topology", default=os.path.join(ROOT, "data", "heipromap_topology.graph"))
    parser.add_argument("--binary", default=os.path.join(ROOT, "build", "Waver-HeiProMap"))
    parser.add_argument("--schedule", default=os.path.join(ROOT, "build", "Waver-Schedule"))
    parser.add_argument("--costs", default="", help="costs.json of the scheduler (default: its built-in values)")
    parser.add_argument("--imbalance", type=float, default=0.03)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--threads", type=int, default=1, help="threads of every Waver-HeiProMap run")
    parser.add_argument("--jobs", type=int, default=4, help="runs in parallel")
    parser.add_argument("--memory-gb", type=int, default=24, help="memory limit of every process (0: none)")
    parser.add_argument("--report-only", action="store_true", help="only write the report from the runs.json of --out")
    parser.add_argument("--show-time", action="store_true", help="add the time of the refinements to the report")
    parser.add_argument("--dose", action="store_true", help="also sweep the parameters of every refinement")
    parser.add_argument("--only", default="", help="comma separated experiments to run (default: all)")
    parser.add_argument("--out", default=os.path.join(ROOT, "ablation_results"))
    parser.add_argument("extra", nargs="*", help="further Waver-HeiProMap options, e.g. -- --spectral-rounds 2")
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    if args.report_only:
        with open(os.path.join(args.out, "runs.json")) as f:
            records = json.load(f)
        order = {}
        for r in records:
            if r["name"] not in order.setdefault(r["experiment"], []):
                order[r["experiment"]].append(r["name"])
        print(report(records, order, os.path.join(args.out, "ablation.md"), args.show_time))
        return
    for path in (args.binary, args.schedule, args.graph, args.topology):
        if not os.path.exists(path):
            sys.exit(f"missing: {path}")

    runs = variants(MODES[args.mode], args.dose)
    if args.only:
        wanted = set(args.only.split(",")) | {"reference"}
        runs = [r for r in runs if r[0] in wanted]
    check_options(args.binary, runs, args.extra)
    order = {}
    for experiment, name, _ in runs:
        order.setdefault(experiment, []).append(name)
    jobs = [(e, n, a, seed) for (e, n, a) in runs for seed in range(args.seeds)]
    print(f"{len(runs)} configurations x {args.seeds} seeds = {len(jobs)} runs, {args.jobs} at a time", flush=True)

    records = []
    raw = os.path.join(args.out, "runs.json")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, job, args) for job in jobs]
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            record = future.result()
            records.append(record)
            status = record.get("error", f"makespan {record.get('makespan', 0):.0f}, cost {record.get('cost', 0):.4g}")
            print(f"[{done}/{len(jobs)}] {record['experiment']} / {record['name']} / seed {record['seed']}: {status}",
                  flush=True)
            with open(raw, "w") as f:
                json.dump(records, f, indent=1)

    text = report(records, order, os.path.join(args.out, "ablation.md"), args.show_time)
    print("\n" + text)
    print(f"raw results: {raw}\nreport: {os.path.join(args.out, 'ablation.md')}")


if __name__ == "__main__":
    main()
