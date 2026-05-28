#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
ENABLE_PROFILER="OFF"
ENABLE_ASSERTS="OFF"
ENABLE_EXCEPTIONS="OFF"
BUILD_TYPE="Release"

BUILD_TESTING="OFF"
RUN_TESTS="OFF"
VERBOSE="OFF"

show_help() {
  echo "Usage: $0 [options]"
  echo "Options:"
  echo "  -v, --verbose     Show detailed build output (default: silent)"
  echo "  -p, --profiler    Enable the profiler (ENABLE_PROFILER=ON)"
  echo "  -a, --asserts     Enable assertions (ENABLE_ASSERTS=ON)"
  echo "  --exceptions      Enable C++ exceptions (default: disabled)"
  echo "  -d, --debug       Build in Debug mode"
  echo "  -t, --test        Build tests (BUILD_TESTING=ON, ENABLE_ASSERTS=ON)"
  echo "  --run-tests       Run tests using gtest-parallel (requires -t)"
  echo "  -h, --help        Show this help message"
  exit 0
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -v|--verbose)
      VERBOSE="ON"
      shift
      ;;
    -p|--profiler)
      ENABLE_PROFILER="ON"
      shift
      ;;
    -a|--asserts)
      ENABLE_ASSERTS="ON"
      shift
      ;;
    --exceptions)
      ENABLE_EXCEPTIONS="ON"
      shift
      ;;
    -d|--debug)
      BUILD_TYPE="Debug"
      shift
      ;;
    -t|--test)
      BUILD_TESTING="ON"
      ENABLE_ASSERTS="ON"
      shift
      ;;
    --run-tests)
      RUN_TESTS="ON"
      shift
      ;;
    -h|--help)
      show_help
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

calc_jobs() {
  local cores
  cores=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
  local j=$((cores - 2))
  if [ "$j" -lt 1 ]; then j=1; fi
  echo "$j"
}
JOBS="${MAX_THREADS:-$(calc_jobs)}"

execute() {
  if [ "${VERBOSE}" == "ON" ]; then
    "$@"
  else
    local log_file
    log_file=$(mktemp)
    if ! "$@" > "$log_file" 2>&1; then
      echo "Error executing: $*" >&2
      cat "$log_file"
      rm -f "$log_file"
      exit 1
    fi
    rm -f "$log_file"
  fi
}

echo "Root: $ROOT"
echo "Jobs: $JOBS"

KAHIP_RELEASE_SO="${ROOT}/extern/local/kahip-release/lib/libkahip.so"
KAHIP_DEBUG_SO="${ROOT}/extern/local/kahip-debug/lib/libkahip.so"

if [ -f "$KAHIP_RELEASE_SO" ] && [ -f "$KAHIP_DEBUG_SO" ]; then
  echo "KaHIP libraries found; skipping download and build."
else
  echo "KaHIP libraries not found; building from source..."

  # -----------------------------
  # Fetch dependencies
  # -----------------------------
  rm -rf extern
  mkdir -p extern

  echo "Downloading KaHIP 3.19..."
  (
    cd extern
    wget -q -O v3.19.tar.gz https://github.com/KaHIP/KaHIP/archive/refs/tags/v3.19.tar.gz
    tar -xzf v3.19.tar.gz
    mv KaHIP-3.19 KaHIP
    rm -f v3.19.tar.gz
  )

  # -----------------------------
  # Build KaHIP Release
  # -----------------------------
  echo "Building KaHIP Release..."
  rm -rf "${ROOT}/extern/KaHIP/build-release"
  mkdir -p "${ROOT}/extern/KaHIP/build-release"

  (
    cd "${ROOT}/extern/KaHIP/build-release"
    execute cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${ROOT}/extern/local/kahip-release" \
      -DNOMPI=ON \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DCMAKE_C_FLAGS_RELEASE="-Ofast -DNDEBUG -march=native" \
      -DCMAKE_CXX_FLAGS_RELEASE="-Ofast -DNDEBUG -march=native"
    execute cmake --build . --target install --parallel "$JOBS"
  )

  # -----------------------------
  # Build KaHIP Debug (Valgrind-friendly)
  # -----------------------------
  echo "Building KaHIP Debug..."
  rm -rf "${ROOT}/extern/KaHIP/build-debug" "${ROOT}/extern/local/kahip-debug"
  mkdir -p "${ROOT}/extern/KaHIP/build-debug"

  (
    cd "${ROOT}/extern/KaHIP/build-debug"
    execute cmake .. \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX="${ROOT}/extern/local/kahip-debug" \
      -DNOMPI=ON \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DCMAKE_C_FLAGS="-O0 -g -march=x86-64 -mtune=generic -fno-omit-frame-pointer -mno-avx512f -mno-avx512vl -mno-avx512dq -mno-avx512bw -mno-avx512cd -mno-avx2 -mno-avx -mno-fma" \
      -DCMAKE_CXX_FLAGS="-O0 -g -march=x86-64 -mtune=generic -fno-omit-frame-pointer -mno-avx512f -mno-avx512vl -mno-avx512dq -mno-avx512bw -mno-avx512cd -mno-avx2 -mno-avx -mno-fma" \
      -DCMAKE_C_FLAGS_DEBUG="" \
      -DCMAKE_CXX_FLAGS_DEBUG=""
    execute cmake --build . --target install --parallel 1 --verbose
  )
fi

# -----------------------------
# TBB
# -----------------------------
TBB_LOCAL="${ROOT}/extern/local/tbb"

if [ -f "${TBB_LOCAL}/lib/libtbb.so" ]; then
  echo "Local TBB found; skipping build."
else
  echo "TBB not found locally; building from source..."
  TBB_VERSION="v2021.13.0"
  (
    cd "${ROOT}/extern"
    wget -q -O tbb.tar.gz "https://github.com/oneapi-src/oneTBB/archive/refs/tags/${TBB_VERSION}.tar.gz"
    tar -xzf tbb.tar.gz
    rm -f tbb.tar.gz
    mv oneTBB-* oneTBB
    execute cmake -S oneTBB -B oneTBB/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${TBB_LOCAL}" \
      -DTBB_TEST=OFF
    execute cmake --build oneTBB/build --target install --parallel "$JOBS"
  )
fi

# -----------------------------
# GoogleTest
# -----------------------------
GTEST_LOCAL="${ROOT}/extern/local/googletest"
if [ -d "${GTEST_LOCAL}" ]; then
  echo "Local GoogleTest found; skipping build."
else
  echo "GoogleTest not found locally; building from source..."
  GTEST_VERSION="v1.14.0"
  (
    cd "${ROOT}/extern"
    wget -q -O googletest.tar.gz "https://github.com/google/googletest/archive/refs/tags/${GTEST_VERSION}.tar.gz"
    tar -xzf googletest.tar.gz
    rm -f googletest.tar.gz
    mv googletest-* googletest-src
    execute cmake -S googletest-src -B googletest-src/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${GTEST_LOCAL}" \
      -DBUILD_GMOCK=OFF \
      -DINSTALL_GTEST=ON
    execute cmake --build googletest-src/build --target install --parallel "$JOBS"
  )
fi

# -----------------------------
# gtest-parallel
# -----------------------------
GP_LOCAL="${ROOT}/extern/gtest-parallel"
if [ -d "${GP_LOCAL}" ]; then
  echo "gtest-parallel found; skipping download."
else
  echo "gtest-parallel not found locally; downloading..."
  (
    cd "${ROOT}/extern"
    execute git clone -q https://github.com/google/gtest-parallel.git
  )
fi

# -----------------------------
# Build HeiProMap
# -----------------------------
echo "Building HeiProMap (${BUILD_TYPE}, Profiler=${ENABLE_PROFILER}, Asserts=${ENABLE_ASSERTS})..."
rm -rf "${ROOT}/build"
mkdir "${ROOT}/build"

CMAKE_EXTRA_ARGS="-DCMAKE_PREFIX_PATH=${TBB_LOCAL}\;${GTEST_LOCAL} -DENABLE_PROFILER=${ENABLE_PROFILER} -DENABLE_ASSERTS=${ENABLE_ASSERTS} -DENABLE_EXCEPTIONS=${ENABLE_EXCEPTIONS} -DBUILD_TESTING=${BUILD_TESTING}"
  execute cmake -S "${ROOT}" -B "${ROOT}/build" -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" ${CMAKE_EXTRA_ARGS}
  execute cmake --build "${ROOT}/build" --parallel "$JOBS" --target HeiProMap
  execute cmake --build "${ROOT}/build" --parallel "$JOBS" --target Dyn-HeiProMap
  execute cmake --build "${ROOT}/build" --parallel "$JOBS" --target HeiPa

  if [ "${BUILD_TESTING}" == "ON" ]; then
    execute cmake --build "${ROOT}/build" --parallel "$JOBS" --target HeiProMapTests
  fi

  if [ "${RUN_TESTS}" == "ON" ]; then
    if [ "${BUILD_TESTING}" == "OFF" ]; then
      echo "Error: Cannot run tests without building them first. Use -t or --test." >&2
      exit 1
    fi
    echo "Running tests with gtest-parallel..."
    python3 "${GP_LOCAL}/gtest-parallel" "${ROOT}/build/tests/HeiProMapTests"
  fi

