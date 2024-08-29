#ifndef HEIDELBERGPROCESSMAPPING_PARALLEL_UTILS_H
#define HEIDELBERGPROCESSMAPPING_PARALLEL_UTILS_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <omp.h>

#include "../../definitions.h"

namespace HeiProMap {

void parallel_read_partition(const std::string &mapping_in,
                             std::vector<partition_t> &partition,
                             u64 n_threads);

    void parallel_write_partition(std::vector<partition_t> &partition,
                                  const std::string &mapping_out,
                                  u64 n_threads);

}

#endif //HEIDELBERGPROCESSMAPPING_PARALLEL_UTILS_H
