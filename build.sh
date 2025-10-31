#!/bin/bash

ROOT=${PWD}

# ----- pick a reasonable parallelism (leave 2 cores free) -----
calc_jobs() {
  local cores
  cores=$( nproc 2>/dev/null \
        || getconf _NPROCESSORS_ONLN 2>/dev/null \
        || sysctl -n hw.ncpu 2>/dev/null \
        || echo 4 )
  local j=$(( cores - 2 ))
  if [ "$j" -lt 1 ]; then j=1; fi
  echo "$j"
}
JOBS="${MAX_THREADS:-$(calc_jobs)}"
echo "Building with $JOBS parallel jobs (override with MAX_THREADS)."

echo "Root          : ${ROOT}"
echo "Using C compiler: ${GCC:-<system default>}"

# make local folder for all includes
mkdir -p extern

# --- Download KaHIP 3.19 ---
echo "Downloading KaHIP 3.19..."
if (
  cd extern \
  && rm -f v3.19.tar.gz \
  && rm -rf KaHIP \
  && wget -q https://github.com/KaHIP/KaHIP/archive/refs/tags/v3.19.tar.gz \
  && tar -xzf v3.19.tar.gz \
  && mv KaHIP-3.19 KaHIP \
  && rm -f v3.19.tar.gz
); then
  echo "KaHIP 3.19 downloaded and extracted successfully."
else
  echo "Failed to download KaHIP!" >&2
  exit 1
fi

# --- Download maxflow 3.04 ---
echo "Downloading maxflow 3.04..."
if (
  cd extern \
  && rm -f maxflow-v3.04.src.zip \
  && rm -rf maxflow-v3.04.src \
  && wget -q https://pub.ista.ac.at/~vnk/software/maxflow-v3.04.src.zip \
  && unzip -q maxflow-v3.04.src.zip \
  && mv maxflow-v3.04.src maxflow \
  && rm -f maxflow-v3.04.src.zip
); then
  echo "maxflow 3.04 downloaded and extracted successfully."
else
  echo "Failed to download maxflow!" >&2
  exit 1
fi

# --- build KaHIP ---
echo "Building KaHIP 3.19..."
if (
  cd "${ROOT}/extern/KaHIP" \
  && mkdir -p build && cd build \
  && cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${ROOT}/extern/local/kahip" \
    -DNOMPI=ON \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    > /dev/null 2>&1 \
  && make install -j "$JOBS" > /dev/null 2>&1
); then
  echo "KaHIP 3.19 build completed successfully."
else
  echo "KaHIP 3.19 build failed!" >&2
  exit 1
fi
cd "${ROOT}"

# make directory
mkdir build
cd build

# build
cmake .. -DCMAKE_BUILD_TYPE=Release && cd ${ROOT}
cmake --build build --parallel "$(get_num_cores)" --target HeiProMap
cmake --build build --parallel "$(get_num_cores)" --target heipromap  # the library
