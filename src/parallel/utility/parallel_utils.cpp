#include "parallel_utils.h"
#include "../../macros.h"
#include "../../serial/utility/utils.h"

namespace HeiProMap {

    void parallel_read_partition(const std::string &mapping_in,
                                 std::vector<partition_t> &partition,
                                 u64 n_threads) {
        char *file_arr = nullptr;
        size_t file_size = 0;
        int fd;

        // Open the file
        fd = open(mapping_in.c_str(), O_RDONLY);
        if (fd == -1) {
            std::cerr << "File " << mapping_in << " does not exist!" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Get the file size
        struct stat fileInfo;
        if (fstat(fd, &fileInfo) == -1) {
            std::cerr << "File " << mapping_in << " Could not get file size!" << std::endl;
            close(fd);
            exit(EXIT_FAILURE);
        }
        file_size = fileInfo.st_size;

        // Memory-map the file
        file_arr = static_cast<char *>(mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (file_arr == MAP_FAILED) {
            std::cerr << "File " << mapping_in << " Could not map the file!" << std::endl;
            close(fd);
            exit(EXIT_FAILURE);
        }

        size_t n = partition.size();
#pragma omp parallel default(none) firstprivate(n, file_arr, file_size, n_threads) shared(partition) num_threads(n_threads)
        {
            size_t thread_id = omp_get_thread_num();
            vertex_t base_range = floor((f64) n / (f64) n_threads);
            vertex_t rem = n % n_threads;

            vertex_t start_u;
            vertex_t end_u;
            if (thread_id < rem) {
                start_u = thread_id * (base_range + 1);
                end_u = start_u + base_range + 1;
            } else {
                start_u = rem * (base_range + 1) + (thread_id - rem) * base_range;
                end_u = start_u + base_range;
            }

            size_t i = 0;
            vertex_t u = 0;
            while(true){

                if (u < start_u) {
                    // this line should not be read by this thread
                    move_while_not(file_arr, i, '\n', file_size);
                    ++i;
                    u += 1;
                    continue;
                }

                if (u >= end_u) {
                    // this thread has read everything
                    break;
                }

                while (i < file_size && file_arr[i] != '\n') {
                    // read in the partition
                    partition_t p = 0;
                    while (i < file_size && file_arr[i] != '\n') { p = p * 10 + (file_arr[i++] - '0'); }
                    partition[u] = p;
                }

                u += 1;
                i += 1;

                if(u == partition.size()){
                    break;
                }
            }
        }

        // Clean up
        munmap(file_arr, file_size);
        close(fd);
    }

    void parallel_write_partition(std::vector<partition_t> &partition, const std::string &mapping_out, u64 n_threads) {
        std::vector<char*> arrs(n_threads);
        std::vector<size_t> sizes(n_threads);
        size_t n = partition.size();

#pragma omp parallel default(none) firstprivate(n, n_threads) shared(partition, arrs, sizes) num_threads(n_threads)
        {
            size_t thread_id = omp_get_thread_num();
            vertex_t base_range = floor((f64) n / (f64) n_threads);
            vertex_t rem = n % n_threads;

            vertex_t start_u;
            vertex_t end_u;
            if (thread_id < rem) {
                start_u = thread_id * (base_range + 1);
                end_u = start_u + base_range + 1;
            } else {
                start_u = rem * (base_range + 1) + (thread_id - rem) * base_range;
                end_u = start_u + base_range;
            }

            // each thread populates its own memory
            size_t max_size =(end_u - start_u)* 16;
            size_t size = 0;
            char* arr = (char *) malloc(max_size * sizeof(char));
            char buffer[16];

            for(vertex_t u = start_u; u < end_u; ++u){
                sprintf(buffer, "%d", u);
                size_t i = 0;
                while(buffer[i] != '\0'){
                    arr[size++] = buffer[i++];
                }
                arr[size++] = '\n';
            }

            arrs[thread_id] = arr;
            sizes[thread_id] = size;
        }

        std::ofstream ofs;
        ofs.open(mapping_out);
        for(size_t i = 0; i < n_threads; ++i) {
            ofs.write(arrs[i], sizes[i]);
        }
        ofs.close();

        for(char* p : arrs){
            free(p);
        }
    }
}
