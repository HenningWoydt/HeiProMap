import os
import json
import subprocess
import sys

def build_heipa(project_root):
    build_script = os.path.join(project_root, "build.sh")
    
    # Build Bench version
    print("Building HeiPa_bench (Release, no profiler)...")
    try:
        subprocess.run([build_script], check=True, capture_output=True)
        bin_path = os.path.join(project_root, "build", "HeiPa")
        bench_path = os.path.join(project_root, "HeiPa_bench")
        if os.path.exists(bench_path): os.remove(bench_path)
        os.rename(bin_path, bench_path)
        print("HeiPa_bench built.")
    except Exception as e:
        print(f"Build failed: {e}")
        sys.exit(1)

    # Build Profile version
    print("Building HeiPa_profile (Release, with profiler)...")
    try:
        subprocess.run([build_script, "-p"], check=True, capture_output=True)
        bin_path = os.path.join(project_root, "build", "HeiPa")
        profile_path = os.path.join(project_root, "HeiPa_profile")
        if os.path.exists(profile_path): os.remove(profile_path)
        os.rename(bin_path, profile_path)
        print("HeiPa_profile built.")
    except Exception as e:
        print(f"Build failed: {e}")
        sys.exit(1)

def aggregate_profile(master, current):
    master["total_ms"] += current["total"]["total_ms"]
    for g_name, g_data in current["groups"].items():
        if g_name not in master["groups"]:
            master["groups"][g_name] = {"total_ms": 0.0, "functions": {}}
        master["groups"][g_name]["total_ms"] += g_data["total_ms"]
        for f_name, f_data in g_data["functions"].items():
            if f_name not in master["groups"][g_name]["functions"]:
                master["groups"][g_name]["functions"][f_name] = {"total_ms": 0.0, "kernels": {}}
            master["groups"][g_name]["functions"][f_name]["total_ms"] += f_data["total_ms"]
            for k_name, k_data in f_data["kernels"].items():
                if k_name not in master["groups"][g_name]["functions"][f_name]["kernels"]:
                    master["groups"][g_name]["functions"][f_name]["kernels"][k_name] = {"calls": 0, "total_ms": 0.0}
                master["groups"][g_name]["functions"][f_name]["kernels"][k_name]["calls"] += k_data["calls"]
                master["groups"][g_name]["functions"][f_name]["kernels"][k_name]["total_ms"] += k_data["total_ms"]

def print_profile_table(profile):
    total_ms = profile["total_ms"]
    if total_ms == 0:
        return

    name_w = 48
    calls_w = 10
    total_w = 12
    avg_w = 10
    pct_w = 7

    def pad(text, w):
        return (text[:w-1] + "…") if len(text) > w else text.ljust(w)

    rule = "-" * (name_w + calls_w + total_w + avg_w + pct_w + 12)
    print("\n" + rule)
    print(f"{pad('Scope', name_w)}   {pad('Calls', calls_w)}   {pad('Total ms', total_w)}   {pad('Avg ms', avg_w)}   {pad('%Tot', pct_w)}")
    print(rule)

    print(f"{pad('TOTAL', name_w)}   {pad('-', calls_w)}   {pad(f'{total_ms:.3f}', total_w)}   {pad('-', avg_w)}   {pad('100.0', pct_w)}")

    sorted_groups = sorted(profile["groups"].items(), key=lambda x: x[1]["total_ms"], reverse=True)

    for g_name, g_data in sorted_groups:
        g_pct = (g_data["total_ms"] / total_ms) * 100
        g_total_str = f"{g_data['total_ms']:.3f}"
        print(f"{pad('+-- [G] ' + g_name, name_w)}   {pad('-', calls_w)}   {pad(g_total_str, total_w)}   {pad('-', avg_w)}   {pad(f'{g_pct:.1f}', pct_w)}")
        
        sorted_functions = sorted(g_data["functions"].items(), key=lambda x: x[1]["total_ms"], reverse=True)
        for f_name, f_data in sorted_functions:
            f_pct = (f_data["total_ms"] / total_ms) * 100
            f_total_str = f"{f_data['total_ms']:.3f}"
            print(f"{pad('|   +-- [F] ' + f_name, name_w)}   {pad('-', calls_w)}   {pad(f_total_str, total_w)}   {pad('-', avg_w)}   {pad(f'{f_pct:.1f}', pct_w)}")
            
            sorted_kernels = sorted(f_data["kernels"].items(), key=lambda x: x[1]["total_ms"], reverse=True)
            for k_name, k_data in sorted_kernels:
                k_pct = (k_data["total_ms"] / total_ms) * 100
                k_avg = k_data["total_ms"] / k_data["calls"] if k_data["calls"] > 0 else 0
                k_total_str = f"{k_data['total_ms']:.3f}"
                k_avg_str = f"{k_avg:.3f}"
                print(f"{pad('|   |   +-- [K] ' + k_name, name_w)}   {pad(str(k_data['calls']), calls_w)}   {pad(k_total_str, total_w)}   {pad(k_avg_str, avg_w)}   {pad(f'{k_pct:.1f}', pct_w)}")
    
    print(rule + "\n")

def main():
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    # Ensure binaries are built correctly
    build_heipa(project_root)

    data_dir = os.path.join(project_root, "data")
    results_dir = os.path.join(data_dir, "results")
    graphs_dir = os.path.join(data_dir, "graphs")
    heipa_bench = os.path.join(project_root, "HeiPa_bench")
    heipa_profile = os.path.join(project_root, "HeiPa_profile")

    if not os.path.exists(results_dir):
        print(f"Error: Results directory {results_dir} not found.")
        sys.exit(1)

    json_files = [f for f in os.listdir(results_dir) if f.endswith(".json")]
    
    if not json_files:
        print("No collected data found in data/results/")
        sys.exit(0)

    results = []
    master_profile = {"total_ms": 0.0, "groups": {}}

    print(f"Running benchmarks on {len(json_files)} instances...")

    for f in json_files:
        with open(os.path.join(results_dir, f), 'r') as jfile:
            data = json.load(jfile)
        
        g_hash = data["graph_hash"]
        k = data["k"]
        imb = data["imbalance"]
        seed = data["seed"]
        dataset_cut = data["edge_cut"]
        dataset_time = data["time_ms"]

        graph_path = os.path.join(graphs_dir, f"{g_hash}.graph")
        if not os.path.exists(graph_path):
            continue

        # 1. Run Bench version
        cmd_bench = [
            heipa_bench,
            "--graph", graph_path,
            "--k", str(k),
            "--imbalance", str(imb),
            "--seed", str(seed),
            "--config", "fast",
            "--json-output", "1"
        ]

        try:
            result = subprocess.run(cmd_bench, capture_output=True, text=True, check=True)
            output = result.stdout
            json_start = output.find('{')
            json_end = output.rfind('}') + 1
            heipa_res = json.loads(output[json_start:json_end])
            
            heipa_cut = heipa_res["edge_cut"]
            heipa_time = heipa_res["time_ms"]

            # Balance check
            if "max_block_weight" in data and "graph_weight" in data:
                lmax_dataset = (1.0 + imb) * (data["graph_weight"] / k)
                dataset_balanced = data["max_block_weight"] <= lmax_dataset + 1
            else:
                dataset_balanced = None
            
            if "max_block_weight" in heipa_res and "lmax" in heipa_res:
                heipa_balanced = heipa_res["max_block_weight"] <= heipa_res["lmax"] + 1
            else:
                heipa_balanced = None

            results.append({
                "dataset_cut": dataset_cut,
                "heipa_cut": heipa_cut,
                "dataset_time": dataset_time,
                "heipa_time": heipa_time,
                "dataset_balanced": dataset_balanced,
                "heipa_balanced": heipa_balanced
            })

            # 2. Run Profile version
            cmd_profile = [
                heipa_profile,
                "--graph", graph_path,
                "--k", str(k),
                "--imbalance", str(imb),
                "--seed", str(seed),
                "--config", "fast",
                "--json-output", "1"
            ]
            
            result_prof = subprocess.run(cmd_profile, capture_output=True, text=True, check=True)
            output_prof = result_prof.stdout
            prof_start = output_prof.find("PROFILER_JSON_START")
            prof_end = output_prof.find("PROFILER_JSON_END")
            if prof_start != -1 and prof_end != -1:
                prof_json_str = output_prof[prof_start + len("PROFILER_JSON_START"):prof_end].strip()
                prof_json = json.loads(prof_json_str)
                aggregate_profile(master_profile, prof_json)

        except subprocess.CalledProcessError:
            pass
        except json.JSONDecodeError:
            pass

    if results:
        n = len(results)
        better_cut = sum(1 for r in results if r["heipa_cut"] < r["dataset_cut"])
        equal_cut = sum(1 for r in results if r["heipa_cut"] == r["dataset_cut"])
        worse_cut = n - better_cut - equal_cut
        faster = sum(1 for r in results if r["heipa_time"] < r["dataset_time"])
        slower = n - faster
        dominated = sum(1 for r in results if r["heipa_cut"] <= r["dataset_cut"] and r["heipa_time"] < r["dataset_time"])
        balanced_dataset = sum(1 for r in results if r["dataset_balanced"] is True)
        balanced_heipa = sum(1 for r in results if r["heipa_balanced"] is True)
        total_balance_check = sum(1 for r in results if r["dataset_balanced"] is not None)
        speedups = [r["dataset_time"] / r["heipa_time"] for r in results if r["heipa_time"] > 0]
        improvements = [(r["dataset_cut"] - r["heipa_cut"]) / r["dataset_cut"] for r in results if r["dataset_cut"] > 0]
        
        import math
        def geo_mean(iterable):
            a = [x for x in iterable if x > 0]
            if not a: return 0
            return math.exp(sum(math.log(x) for x in a) / len(a))

        g_speedup = geo_mean(speedups)

        print("\n" + "="*30 + " BENCHMARK SUMMARY " + "="*30)
        if total_balance_check > 0:
            print(f"  Dataset Stability: {balanced_dataset/total_balance_check:>6.1%} ({balanced_dataset}/{total_balance_check} valid)")
        print(f"  HeiPa Stability:   {balanced_heipa/n:>6.1%} ({balanced_heipa}/{n} valid)")

        print(f"\n[Quality Metrics (Edge Cut)]")
        print(f"  Better:  {better_cut:>3} ({better_cut/n:>5.1%})")
        print(f"  Equal:   {equal_cut:>3} ({equal_cut/n:>5.1%})")
        print(f"  Worse:   {worse_cut:>3} ({worse_cut/n:>5.1%})")
        if improvements:
            print(f"  Avg Improvement: {sum(improvements)/len(improvements):>6.2%}")
        
        print(f"\n[Performance Metrics (Time)]")
        print(f"  Faster:  {faster:>3} ({faster/n:>5.1%})")
        print(f"  Slower:  {slower:>3} ({slower/n:>5.1%})")
        print(f"  Strict Dominance: {dominated}/{n} (HeiPa is strictly better than Dataset)")
        print(f"  Geom. Mean Speedup:   {g_speedup:.2f}x")

        print("\n" + "="*30 + " AGGREGATED PROFILING " + "="*30)
        print_profile_table(master_profile)
        print("="*79)

if __name__ == "__main__":
    main()
