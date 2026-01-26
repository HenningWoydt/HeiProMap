#!/bin/bash

ROOT=${PWD}
GCC=$(which gcc || true)

echo "Root          : ${ROOT}"
echo "Using C compiler: ${GCC:-<system default>}"

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

# make local folder for all includes
rm -rf extern
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
  && rm -rf maxflow \
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

# --- Download GKlib (latest release) ---
echo "Downloading GKlib..."
if (
  cd extern \
  && rm -f gklib.tar.gz \
  && rm -rf GKlib \
  && wget -q https://github.com/KarypisLab/GKlib/archive/refs/heads/master.tar.gz -O gklib.tar.gz \
  && tar -xzf gklib.tar.gz \
  && mv GKlib-master GKlib \
  && rm -f gklib.tar.gz
); then
  echo "GKlib downloaded and extracted successfully."
else
  echo "Failed to download GKlib!" >&2
  exit 1
fi

# --- Download METIS 5.2.1 ---
echo "Downloading METIS 5.2.1..."
if (
  cd extern \
  && rm -f metis-5.2.1.tar.gz \
  && rm -rf METIS \
  && wget -q https://github.com/KarypisLab/METIS/archive/refs/tags/v5.2.1.tar.gz -O metis-5.2.1.tar.gz \
  && tar -xzf metis-5.2.1.tar.gz \
  && mv METIS-5.2.1 METIS \
  && rm -f metis-5.2.1.tar.gz
); then
  echo "METIS 5.2.1 downloaded and extracted successfully."
else
  echo "Failed to download METIS v5.2.1!" >&2
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


# --- Download Mt-KaHyPar 1.5.3 ---
echo "Downloading Mt-KaHyPar 1.5.3..."
if (
  cd extern \
  && rm -f v1.5.3.tar.gz \
  && rm -rf MtKaHyPar \
  && wget -q https://github.com/kahypar/mt-kahypar/archive/refs/tags/v1.5.3.tar.gz \
  && tar -xzf v1.5.3.tar.gz \
  && mv mt-kahypar-1.5.3 MtKaHyPar \
  && rm -f v1.5.3.tar.gz
); then
  echo "Mt-KaHyPar 1.5.3 downloaded and extracted successfully."
else
  echo "Failed to download Mt-KaHyPar 1.5.3!" >&2
  exit 1
fi
cd "${ROOT}"

# install GKLIB into a local folder
export CFLAGS="-Wall -Wno-error=pedantic -Wno-error -D_GNU_SOURCE -DHAVE_STRDUP=1"
export CPPFLAGS="-Wall -Wno-error=pedantic -Wno-error -D_GNU_SOURCE -DHAVE_STRDUP=1"

echo "Building GKlib..."
if cd "${ROOT}/extern/GKlib" && rm -rf build \
  && make config prefix="${ROOT}/extern/local/gklib" cc="${GCC}" \
  && make install; then
  echo "GKlib build completed successfully."
else
  echo "GKlib build failed!" >&2
  exit 1
fi
cd "${ROOT}"

echo "Building METIS..."
if cd "${ROOT}/extern/METIS" \
  && rm -rf build \
  && make config prefix="${ROOT}/extern/local/metis" gklib_path="${ROOT}/extern/local/gklib" cc="${GCC}" > /dev/null 2>&1 \
  && make install > /dev/null 2>&1; then
  echo "METIS build completed successfully."
else
  echo "METIS build failed!" >&2
  exit 1
fi
cd "${ROOT}"

# -------------------------
# Build / install Mt-KaHyPar
# -------------------------
echo "Building Mt-KaHyPar..."
(
  cd "${ROOT}/extern/MtKaHyPar"
  rm -rf build
  mkdir -p build
  cd build

  cmake .. \
    -DCMAKE_INSTALL_PREFIX="${ROOT}/extern/local/mt-kahypar" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DKAHYPAR_DOWNLOAD_TBB=ON \
    -DKAHYPAR_DOWNLOAD_BOOST=ON \
    -DKAHYPAR_USE_64_BIT_IDS=ON \
    -DKAHYPAR_ENABLE_THREAD_PINNING=OFF \
    -DKAHYPAR_DISABLE_ASSERTIONS=ON \
    -DCMAKE_C_COMPILER="${GCC:-gcc}" \
    -DCMAKE_CXX_COMPILER="${GXX:-g++}"
cmake --build . --parallel "$JOBS" --target install-mtkahypar
)
echo "Mt-KaHyPar build completed successfully."

# --- bundle Mt-KaHyPar downloaded runtime deps into the install prefix ---
echo "Bundling Mt-KaHyPar downloaded Boost/TBB into prefix..."

MTK_BUILD="${ROOT}/extern/MtKaHyPar/build"
MTK_PREFIX="${ROOT}/extern/local/mt-kahypar"

# Mt-KaHyPar may install into lib64 on some systems (common on HPC)
if [ -d "${MTK_PREFIX}/lib64" ]; then
  MTK_LIBDIR="${MTK_PREFIX}/lib64"
else
  MTK_LIBDIR="${MTK_PREFIX}/lib"
fi
mkdir -p "${MTK_LIBDIR}"

# Locate downloaded shared libs in the build tree
TBB_SO=$(find "${MTK_BUILD}" -type f -name "libtbb.so*" | head -n 1)
TBBMALLOC_SO=$(find "${MTK_BUILD}" -type f -name "libtbbmalloc.so*" | head -n 1)

# --- Boost libs (downloaded by KAHYPAR_DOWNLOAD_BOOST=ON) ---
BOOST_PO_SO=$(find "${MTK_BUILD}" -type f -name "libboost_program_options.so*" | head -n 1)
BOOST_SYS_SO=$(find "${MTK_BUILD}" -type f -name "libboost_system.so*" | head -n 1)
BOOST_FS_SO=$(find "${MTK_BUILD}" -type f -name "libboost_filesystem.so*" | head -n 1)
BOOST_THR_SO=$(find "${MTK_BUILD}" -type f -name "libboost_thread.so*" | head -n 1)
BOOST_CONT_SO=$(find "${MTK_BUILD}" -type f -name "libboost_container.so*" | head -n 1)

if [ -z "${TBB_SO}" ] || [ -z "${TBBMALLOC_SO}" ] || [ -z "${BOOST_PO_SO}" ]; then
  echo "ERROR: Could not locate required downloaded libs in ${MTK_BUILD}" >&2
  echo "  TBB_SO=${TBB_SO}" >&2
  echo "  TBBMALLOC_SO=${TBBMALLOC_SO}" >&2
  echo "  BOOST_PO_SO=${BOOST_PO_SO}" >&2
  echo "Hint: find ${MTK_BUILD} -name 'libboost_*.so*' -o -name 'libtbb*.so*'" >&2
  exit 1
fi

# Copy TBB + tbbmalloc (with symlinks)
cp -a "$(dirname "${TBB_SO}")"/libtbb.so* "${MTK_LIBDIR}/"
cp -a "$(dirname "${TBBMALLOC_SO}")"/libtbbmalloc.so* "${MTK_LIBDIR}/"

# Copy Boost program_options (+ common deps if present)
cp -a "$(dirname "${BOOST_PO_SO}")"/libboost_program_options.so* "${MTK_LIBDIR}/"
[ -n "${BOOST_SYS_SO}" ] && cp -a "$(dirname "${BOOST_SYS_SO}")"/libboost_system.so* "${MTK_LIBDIR}/"
[ -n "${BOOST_FS_SO}" ]  && cp -a "$(dirname "${BOOST_FS_SO}")"/libboost_filesystem.so* "${MTK_LIBDIR}/"
[ -n "${BOOST_THR_SO}" ] && cp -a "$(dirname "${BOOST_THR_SO}")"/libboost_thread.so* "${MTK_LIBDIR}/"
[ -n "${BOOST_CONT_SO}" ] && cp -a "$(dirname "${BOOST_CONT_SO}")"/libboost_container.so* "${MTK_LIBDIR}/"


# Optional: ensure mt-kahypar itself searches next to itself at runtime
if command -v patchelf >/dev/null 2>&1; then
  if [ -f "${MTK_LIBDIR}/libmtkahypar.so.1.5.3" ]; then
    echo "Setting RPATH on libmtkahypar to \$ORIGIN"
    patchelf --set-rpath '$ORIGIN' "${MTK_LIBDIR}/libmtkahypar.so.1.5.3" || true
  fi
fi

echo "Bundled libs in ${MTK_LIBDIR}:"
ls -1 "${MTK_LIBDIR}" | egrep 'mtkahypar|libboost_|libtbb' || true

# make directory
rm -rf build
mkdir build
cd build

# build
cmake .. -DCMAKE_BUILD_TYPE=Release && cd ${ROOT}
cmake --build build --parallel "$JOBS" --target HeiProMap
cmake --build build --parallel "$JOBS" --target heipromap  # the library
