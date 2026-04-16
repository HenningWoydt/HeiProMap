#ifndef HEIPROMAP_MEMORY_STACK_H
#define HEIPROMAP_MEMORY_STACK_H

#include <cstdint>
#include <cstdlib>

namespace HeiProMap {
    inline uint64_t align64(uint64_t n) { return (n + 63) & ~uint64_t(63); }

    class MemoryStack {
        static constexpr size_t ALIGNMENT = 64;
        char *buffer = nullptr;
        size_t capacity = 0;
        size_t offset = 0;

    public:
        ~MemoryStack() { std::free(buffer); }

        MemoryStack() = default;

        MemoryStack(const MemoryStack &) = delete;

        MemoryStack &operator=(const MemoryStack &) = delete;

        MemoryStack(MemoryStack &&other) noexcept: buffer(other.buffer), capacity(other.capacity), offset(other.offset) {
            other.buffer = nullptr;
            other.capacity = 0;
            other.offset = 0;
        }

        MemoryStack &operator=(MemoryStack &&other) noexcept {
            if (this != &other) {
                std::free(buffer);
                buffer = other.buffer;
                capacity = other.capacity;
                offset = other.offset;
                other.buffer = nullptr;
                other.capacity = 0;
                other.offset = 0;
            }
            return *this;
        }

        void ensure(size_t total_bytes) {
            total_bytes = align64(total_bytes);
            if (total_bytes > capacity) {
                std::free(buffer);
                buffer = static_cast<char *>(std::aligned_alloc(ALIGNMENT, total_bytes));
                capacity = total_bytes;
            }
        }

        void *get_memory(size_t size) {
            size = align64(size);
            void *ptr = buffer + offset;
            offset += size;
            return ptr;
        }

        void clear() { offset = 0; }
    };
}

#endif // HEIPROMAP_MEMORY_STACK_H
