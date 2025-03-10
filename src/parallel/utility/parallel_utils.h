/*******************************************************************************
 * MIT License
 *
 * This file is part of HeiProMap.
 *
 * Copyright (C) 2025 Henning Woydt <henning.woydt@informatik.uni-heidelberg.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef HEIPROMAP_PARALLEL_UTILS_H
#define HEIPROMAP_PARALLEL_UTILS_H

#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "../../definitions.h"
#include "../../macros.h"
#include "../../commons/utils.h"

namespace HeiProMap {

void parallel_read_partition(const std::string &mapping_in,
                             std::vector<partition_t> &partition,
                             u64 n_threads);

    void parallel_write_partition(std::vector<partition_t> &partition,
                                  const std::string &mapping_out,
                                  u64 n_threads);

}

#endif //HEIPROMAP_PARALLEL_UTILS_H
