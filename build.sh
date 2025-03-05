#!/bin/bash

ROOT=${PWD}

function get_num_cores {
  grep -c ^processor /proc/cpuinfo;
}

# update all submodules
git submodule update --init --recursive

# make directory
mkdir build
cd build

# build
cmake .. -DCMAKE_BUILD_TYPE=Release && cd ${ROOT}
cmake --build build --parallel "$(get_num_cores)" --target HeiProMap
cmake --build build --parallel "$(get_num_cores)" --target MtHeiProMap
cmake --build build --parallel "$(get_num_cores)" --target MtHeiProMapRefinement
# cmake --build build --parallel "$(get_num_cores)" --target spm_gtest
# cd ${ROOT}

# run tests
# cd build
# ./spm_gtest
