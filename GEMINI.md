# HeiProMap Agent Instructions

This file provides foundational guidance for AI agents (like Gemini CLI) to effectively build, use, and test the HeiProMap framework.

## Project Architecture
- **Header-Only Library**: The core logic is implemented in `src/` and exposed via the `HeiProMapLib` interface library.
- **Dependencies**: The project depends on **KaHIP** and **TBB**. The `build.sh` script is responsible for downloading and building these locally in `extern/local/`.
- **Parallelism**: Uses OpenMP and Intel TBB.

## Build Procedures
The primary way to build the project is through the `build.sh` script.

- **Standard Build**:
  ```bash
  ./build.sh
  ```
  *Note: Building is silent by default. Use `-v` or `--verbose` to see detailed output.*

- **Debug Build**:
  ```bash
  ./build.sh -d
  ```

- **Enable Profiler/Asserts**:
  ```bash
  ./build.sh -p  # Enable Profiler
  ./build.sh -a  # Enable Asserts
  ```

## Usage and Use Cases
Refer to the `README.md` for detailed command-line options. The main executables are:
- **`HeiProMap`**: Process mapping (guest graph to architecture hierarchy).
- **`HeiPa`**: High-quality graph partitioning (edge-cut minimization).
- **`Dyn-HeiProMap`**: Incremental mapping for dynamic graphs.

## Development Mandates
- **Exception Handling**: The project is compiled with **C++ exceptions disabled** by default (`-fno-exceptions`). They can be re-enabled for specific builds with the `--exceptions` flag.
- **Header Changes**: Since the library is header-only, changes to `src/*.h` files will affect all executables.
- **Safety**: Do not modify the `extern/` directory manually; use the `build.sh` script to manage dependencies.
