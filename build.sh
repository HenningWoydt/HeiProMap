#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
for arg in "$@"; do
  echo "Unknown argument: $arg" >&2
  exit 1
done

calc_jobs() {
  local cores
  cores=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
  local j=$((cores - 2))
  if [ "$j" -lt 1 ]; then j=1; fi
  echo "$j"
}
JOBS="${MAX_THREADS:-$(calc_jobs)}"

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
    cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${ROOT}/extern/local/kahip-release" \
      -DNOMPI=ON \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DCMAKE_C_FLAGS_RELEASE="-Ofast -DNDEBUG -march=native" \
      -DCMAKE_CXX_FLAGS_RELEASE="-Ofast -DNDEBUG -march=native"
    cmake --build . --target install --parallel "$JOBS"
  )

  # -----------------------------
  # Build KaHIP Debug (Valgrind-friendly)
  # -----------------------------
  echo "Building KaHIP Debug..."
  rm -rf "${ROOT}/extern/KaHIP/build-debug" "${ROOT}/extern/local/kahip-debug"
  mkdir -p "${ROOT}/extern/KaHIP/build-debug"

  (
    cd "${ROOT}/extern/KaHIP/build-debug"
    cmake .. \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX="${ROOT}/extern/local/kahip-debug" \
      -DNOMPI=ON \
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
      -DCMAKE_C_FLAGS="-O0 -g -march=x86-64 -mtune=generic -fno-omit-frame-pointer -mno-avx512f -mno-avx512vl -mno-avx512dq -mno-avx512bw -mno-avx512cd -mno-avx2 -mno-avx -mno-fma" \
      -DCMAKE_CXX_FLAGS="-O0 -g -march=x86-64 -mtune=generic -fno-omit-frame-pointer -mno-avx512f -mno-avx512vl -mno-avx512dq -mno-avx512bw -mno-avx512cd -mno-avx2 -mno-avx -mno-fma" \
      -DCMAKE_C_FLAGS_DEBUG="" \
      -DCMAKE_CXX_FLAGS_DEBUG=""
    cmake --build . --target install --parallel 1 --verbose
  )
fi

# -----------------------------
# TBB
# -----------------------------
TBB_LOCAL="${ROOT}/extern/local/tbb"

if ldconfig -p 2>/dev/null | grep -q libtbb.so; then
  echo "System TBB found; skipping local build."
elif [ -f "${TBB_LOCAL}/lib/libtbb.so" ]; then
  echo "Local TBB found; skipping build."
else
  echo "TBB not found; building from source..."
  TBB_VERSION="v2021.13.0"
  (
    cd "${ROOT}/extern"
    wget -q -O tbb.tar.gz "https://github.com/oneapi-src/oneTBB/archive/refs/tags/${TBB_VERSION}.tar.gz"
    tar -xzf tbb.tar.gz
    rm -f tbb.tar.gz
    mv oneTBB-* oneTBB
    cmake -S oneTBB -B oneTBB/build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX="${TBB_LOCAL}" \
      -DTBB_TEST=OFF
    cmake --build oneTBB/build --target install --parallel "$JOBS"
  )
fi

# -----------------------------
# Build HeiProMap Release
# -----------------------------
echo "Building HeiProMap Release..."
rm -rf "${ROOT}/build"
mkdir "${ROOT}/build"

CMAKE_EXTRA_ARGS=""
if [ -f "${TBB_LOCAL}/lib/libtbb.so" ]; then
  CMAKE_EXTRA_ARGS="-DCMAKE_PREFIX_PATH=${TBB_LOCAL}"
fi

cmake -S "${ROOT}" -B "${ROOT}/build" -DCMAKE_BUILD_TYPE=Release ${CMAKE_EXTRA_ARGS}
cmake --build "${ROOT}/build" --parallel "$JOBS" --target HeiProMap
