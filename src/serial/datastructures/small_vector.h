#ifndef HEIDELBERGPROCESSMAPPING_SMALL_VECTOR_H
#define HEIDELBERGPROCESSMAPPING_SMALL_VECTOR_H

#include <cstddef>
#include <cstdlib>
#include <vector>
#include <cstring>

#include "arena_allocator.h"

namespace HeiProMap {

    template<typename T>
    class SmallVector {
    private:
        T *m_start, *m_end;
        size_t m_capacity;

    public:
        SmallVector(){
            m_start = nullptr;
            m_end = nullptr;
            m_capacity = 0;
        };

        SmallVector(T* t_start, T* t_end) {
            m_start = t_start;
            m_end = t_end;
            m_capacity = std::distance(t_start, t_end);
        }

        T* begin(){ return m_start; }
        T* end(){ return m_end; }
        const T* begin() const { return m_start; }
        const T* end() const { return m_end; }
        T operator [](size_t i) const {return m_start[i];}
        T & operator [](size_t i) {return m_start[i];}

        size_t size() const { return std::distance(m_start, m_end); }

        void push_back(const T &t, ArenaAllocator<T> &t_alloc) {
            if (size() == m_capacity) {
                size_t new_capacity = std::max((size_t) 1, m_capacity * 2);
                T *new_addr = t_alloc.request(new_capacity);

                memcpy(new_addr, m_start, m_capacity * sizeof(T));
                m_start = new_addr;
                m_end = new_addr + m_capacity;
                m_capacity = new_capacity;
            }
            *m_end = t;
            m_end += 1;
        }

        void insert(size_t pos, const T& value, ArenaAllocator<T> &t_alloc) {
            if (size() == m_capacity) {
                size_t new_capacity = std::max((size_t) 1, m_capacity * 2);
                T *new_addr = t_alloc.request(new_capacity);

                memcpy(new_addr, m_start, pos * sizeof(T));
                new_addr[pos] = value;
                memcpy(new_addr + pos+1, m_start + pos, (m_capacity - pos) * sizeof(T));

                // Update pointers and capacity
                m_start = new_addr;
                m_end = new_addr + m_capacity + 1;
                m_capacity = new_capacity;
            } else {
                for (size_t i = size(); i > pos; --i) { m_start[i] = m_start[i - 1]; }
                m_start[pos] = value;
                m_end++;
            }
        }

        void erase(size_t pos) {
            for (size_t i = pos; i < size() - 1; ++i) { m_start[i] = m_start[i + 1]; }
            m_end--;
        }
    };

}

#endif //HEIDELBERGPROCESSMAPPING_SMALL_VECTOR_H
