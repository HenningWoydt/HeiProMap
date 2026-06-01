import os
import json
import subprocess
import sys

def build_heipa(project_root):
    print("Building HeiPa (Release, no asserts, no profiler)...")
    build_script = os.path.join(project_root, "build.sh")
    try:
        # Standard build is Release, no asserts, no profiler
        subprocess.run([build_script], check=True, capture_output=True)
        print("Build successful.")
    except subprocess.CalledProcessError as e:
        print(f"Build failed: {e}")
        print(e.stderr)
        sys.exit(1)

def main():
    project_root = os.path.dirname(os.path.abspath(__file__))
    
    # Ensure binary is built correctly
    build_heipa(project_root)

    data_dir = os.path.join(project_root, "data")
    results_dir = os.path.join(data_dir, "results")
    graphs_dir = os.path.join(data_dir, "graphs")
    heipa_bin = os.path.join(project_root, "build", "HeiPa")

    if not os.path.exists(results_dir):
        print(f"Error: Results directory {results_dir} not found.")
        sys.exit(1)

    json_files = [f for f in os.listdir(results_dir) if f.endswith(".json")]
    
    if not json_files:
        print("No collected data found in data/results/")
        sys.exit(0)

    results = []
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

        # Run HeiPa in fast mode
        cmd = [
            heipa_bin,
            "--graph", graph_path,
            "--k", str(k),
            "--imbalance", str(imb),
            "--seed", str(seed),
            "--config", "fast",
            "--json-output", "1"
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Find the JSON block in the output (in case there's other text)
            output = result.stdout
            json_start = output.find('{')
            json_end = output.rfind('}') + 1
            heipa_res = json.loads(output[json_start:json_end])
            
            heipa_cut = heipa_res["edge_cut"]
            heipa_time = heipa_res["time_ms"]

            # Balance check
            # For Dataset
            if "max_block_weight" in data and "graph_weight" in data:
                lmax_dataset = (1.0 + imb) * (data["graph_weight"] / k)
                dataset_balanced = data["max_block_weight"] <= lmax_dataset + 1 # small epsilon
            else:
                dataset_balanced = None
            
            # For HeiPa
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

        except subprocess.CalledProcessError as e:
            pass
        except json.JSONDecodeError:
            pass

    if results:
        n = len(results)
        
        # Basic counts
        better_cut = sum(1 for r in results if r["heipa_cut"] < r["dataset_cut"])
        equal_cut = sum(1 for r in results if r["heipa_cut"] == r["dataset_cut"])
        worse_cut = n - better_cut - equal_cut
        
        faster = sum(1 for r in results if r["heipa_time"] < r["dataset_time"])
        slower = n - faster
        
        # Strict Dominance (Faster AND Better/Equal Cut)
        dominated = sum(1 for r in results if r["heipa_cut"] <= r["dataset_cut"] and r["heipa_time"] < r["dataset_time"])
        
        # Balance Stats
        balanced_dataset = sum(1 for r in results if r["dataset_balanced"] is True)
        balanced_heipa = sum(1 for r in results if r["heipa_balanced"] is True)
        total_balance_check = sum(1 for r in results if r["dataset_balanced"] is not None)

        # Ratios
        speedups = [r["dataset_time"] / r["heipa_time"] for r in results if r["heipa_time"] > 0]
        improvements = [(r["dataset_cut"] - r["heipa_cut"]) / r["dataset_cut"] for r in results if r["dataset_cut"] > 0]
        
        # Hall of Fame
        best_speedup = max(speedups) if speedups else 0
        best_improvement = max(improvements) if improvements else 0
        
        # Geometric Mean for Speedup (more robust for ratios)
        import math
        def geo_mean(iterable):
            a = [x for x in iterable if x > 0]
            if not a: return 0
            return math.exp(sum(math.log(x) for x in a) / len(a))

        g_speedup = geo_mean(speedups)

        print("\n" + "="*30 + " BENCHMARK SUMMARY " + "="*30)
        
        print(f"\n[Balance Scorecard]")
        if total_balance_check > 0:
            print(f"  Dataset Stability: {balanced_dataset/total_balance_check:>6.1%} ({balanced_dataset}/{total_balance_check} valid)")
        print(f"  HeiPa Stability:   {balanced_heipa/n:>6.1%} ({balanced_heipa}/{n} valid)")

        print(f"\n[Quality Metrics (Edge Cut)]")
        print(f"  Better:  {better_cut:>3} ({better_cut/n:>5.1%})")
        print(f"  Equal:   {equal_cut:>3} ({equal_cut/n:>5.1%})")
        print(f"  Worse:   {worse_cut:>3} ({worse_cut/n:>5.1%})")
        if improvements:
            print(f"  Avg Improvement: {sum(improvements)/len(improvements):>6.2%}")
            print(f"  Best Reduction:  {best_improvement:>6.2%}")
        
        worse_improvements = [i for i in improvements if i < 0]
        if worse_improvements:
            print(f"  Mean Regression (Worse instances): {sum(worse_improvements)/len(worse_improvements):>6.2%}")

        print(f"\n[Performance Metrics (Time)]")
        print(f"  Faster:  {faster:>3} ({faster/n:>5.1%})")
        print(f"  Slower:  {slower:>3} ({slower/n:>5.1%})")
        print(f"  Strict Dominance: {dominated}/{n} (HeiPa is strictly better than Dataset)")
        print(f"  Arith. Mean Speedup: {sum(speedups)/max(1,len(speedups)):.2f}x")
        print(f"  Geom. Mean Speedup:   {g_speedup:.2f}x")
        
        slowdowns = [1.0/s for s in speedups if s < 1.0]
        if slowdowns:
            print(f"  Mean Slowdown (Worse instances):   {sum(slowdowns)/len(slowdowns):.2f}x")
        
        print(f"\n[Distribution (Speedup)]")
        ranges = [
            ("Super Fast (>3x)", lambda x: x > 3),
            ("Noticeable (1.5x-3x)", lambda x: 1.5 < x <= 3),
            ("Similar (0.8x-1.5x)", lambda x: 0.8 <= x <= 1.5),
            ("Slower (<0.8x)", lambda x: x < 0.8),
        ]
        for label, func in ranges:
            count = sum(1 for s in speedups if func(s))
            bar = "█" * int((count/n)*20) if n > 0 else ""
            print(f"  {label:<20} | {count:>3} | {bar}")

        print("\n" + "="*79)

if __name__ == "__main__":
    main()
