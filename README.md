# HeiProMap: Heidelberg Process Mapping

HeiProMap is a high-performance C++ framework for graph partitioning and process mapping. It implements a multilevel approach to solve the mapping problem, where a guest graph (representing communication patterns) is mapped onto a host graph (representing a parallel architecture) to minimize communication costs.

## Features

- **Multilevel Framework**: Supports coarsening, initial partitioning, and refinement phases.
- **Coarsening Algorithms**:
  - Global Path Algorithm
  - Size-Constrained Label Propagation
  - Heavy Edge Matching
- **Partitioning Strategies**:
  - Global Multisection
  - Recursive Bisection
  - Greedy Partitioner
  - Integration with KaHIP and HeiPa
- **Refinement Techniques**:
  - Label Propagation Refinement
  - Quotient Graph Refinement
  - Flow-Based Refinement
- **Parallelism**: Leverages OpenMP and Intel TBB for multi-core performance.
- **Configurable Profiles**: Predefined settings for `fast`, `eco`, `strong`, and `super-strong` optimization levels.

## Dependencies

- **C++17** compatible compiler (GCC, Clang)
- **CMake** (version 3.16 or higher)
- **OpenMP**
- **Intel TBB**
- **KaHIP** (Heidelberg Graph Partitioning)

## Building

To build the project, use the provided `build.sh` script or standard CMake commands:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

## Usage

HeiProMap provides several executables tailored for different graph-related optimization tasks:

- **`HeiProMap`**: The primary **process mapping** solver.
  - **Use Case**: Mapping a guest graph (e.g., communication patterns of a parallel application) onto a hierarchical host graph (e.g., a supercomputer's network of nodes, sockets, and cores).
  - **Goal**: Minimize the total communication cost (Quadratic Assignment Problem) by placing frequently communicating processes closer together in the architecture hierarchy.

- **`HeiPa`**: A high-quality **graph partitioner**.
  - **Use Case**: Traditional graph partitioning where an input graph needs to be divided into $k$ roughly equal-sized blocks.
  - **Goal**: Minimize the number of edges cut between different blocks (edge-cut objective). It serves as a standalone, high-performance alternative to tools like METIS or KaHIP.

- **`Dyn-HeiProMap`**: A solver for **dynamic graph mapping**.
  - **Use Case**: Scenarios where the guest graph changes over time (e.g., dynamic mesh refinement, evolving social networks, or time-varying communication patterns).
  - **Goal**: Efficiently update the mapping after graph changes without recomputing everything from scratch. It features an interactive shell and batch processing for incremental updates and refinement.

### Basic Example (Process Mapping)

```bash
./HeiProMap --graph <path_to_graph> --hierarchy 4:8:6 --distance 1:10:100 --config eco --threads 8
```

### Command Line Options

- `--graph`, `-g`: Path to the input graph file (METIS format).
- `--hierarchy`, `-h`: Target architecture hierarchy (e.g., `nodes:sockets:cores`).
- `--distance`, `-d`: Communication distances between hierarchy levels.
- `--config`, `-c`: Optimization profile (`fast`, `eco`, `strong`, `super-strong`).
- `--threads`, `-t`: Number of threads to use.
- `--imbalance`, `-e`: Allowed imbalance (default: `0.03`).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

**Henning Woydt**  
Email: [henning.woydt@informatik.uni-heidelberg.de](mailto:henning.woydt@informatik.uni-heidelberg.de)  
GitHub: [https://github.com/HenningWoydt/HeiProMap](https://github.com/HenningWoydt/HeiProMap)
