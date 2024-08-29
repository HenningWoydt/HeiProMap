#ifndef HEIDELBERGPROCESSMAPPING_ARENA_ALLOCATOR_H
#define HEIDELBERGPROCESSMAPPING_ARENA_ALLOCATOR_H

#include <cstddef>
#include <cstdlib>
#include <vector>
#include <cstring>

namespace HeiProMap {

    template<typename T>
    class ArenaAllocator {
    private:
        std::vector<T *> previous_base_addr;
        T *base_addr = nullptr;
        size_t capacity = 0;
        size_t used_size = 0;


    public:

        ~ArenaAllocator() {
            for (T *addr: previous_base_addr) { free(addr); }
            free(base_addr);
        }

        // initialize for at least n elements
        void initialize(size_t n) {
            capacity = ((n + 64 - 1) / 64) * 64; // round to nearest multiple of 64
            base_addr = (T *) malloc(capacity * sizeof(T));
        }

        T* get_base_addr(){
            return base_addr;
        }


        inline T* request(size_t n) {
            T* addr = nullptr;

            if (used_size + n <= capacity) {
                addr = base_addr + used_size;
                used_size += n;
            } else {
                previous_base_addr.push_back(base_addr);
                base_addr = (T *) malloc(capacity * sizeof(T));
                used_size = 0;

                addr = base_addr + used_size;
                used_size += n;
            }

            return addr;
        }

        T operator [](size_t i) const {return base_addr[i];}
        T & operator [](size_t i) {return base_addr[i];}
    };

}

#endif //HEIDELBERGPROCESSMAPPING_ARENA_ALLOCATOR_H
