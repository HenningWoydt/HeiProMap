#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
FAST=0

for arg in "$@"; do
  case "$arg" in
    --fast)
      FAST=1
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      echo "Usage: $0 [--fast]" >&2
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

echo "Root: $ROOT"
echo "Jobs: $JOBS"
echo "Fast: $FAST"

if [ "$FAST" -eq 0 ]; then
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
else
  echo "Fast mode enabled; skipping KaHIP download and build."
fi

# -----------------------------
# Build HeiProMap Release
# -----------------------------
echo "Building HeiProMap Release..."
rm -rf "${ROOT}/build"
mkdir "${ROOT}/build"
cmake -S "${ROOT}" -B "${ROOT}/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "${ROOT}/build" --parallel "$JOBS" --target HeiProMap
