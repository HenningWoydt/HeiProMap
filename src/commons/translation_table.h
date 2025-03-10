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

#ifndef HEIPROMAP_TRANSLATION_TABLE_H
#define HEIPROMAP_TRANSLATION_TABLE_H

#include <unordered_map>

#include "../definitions.h"
#include "../macros.h"

namespace HeiProMap {
    template<typename T>
    class TranslationTable {
    private:
        T* m_translation_o_to_n = nullptr;
        T* m_translation_n_to_o = nullptr;

    public:
        /**
         * Default constructor.
         */
        TranslationTable() = default;

        ~TranslationTable() {
            free(m_translation_o_to_n);
            free(m_translation_n_to_o);
        }

        /**
         * Adds a translation from o to n and from n to o.
         *
         * @param o Old value.
         * @param n New value.
         */
        void add(const T o, const T n) {
            m_translation_o_to_n[o] = n;

            m_translation_n_to_o[n] = o;
        }

        void reserve(size_t n_space, size_t o_space){
            free(m_translation_n_to_o);
            free(m_translation_o_to_n);

            size_t n_space_64 = round_up_64(n_space + 1);
            size_t o_space_64 = round_up_64(o_space + 1);
            m_translation_n_to_o = (T*) aligned_alloc(64, n_space_64 * sizeof(T));
            m_translation_o_to_n = (T*) aligned_alloc(64, o_space_64 * sizeof(T));
        }

        /**
         * Get the new value for o.
         *
         * @param o Old value
         */
        T get_n(const T o) const {
            return m_translation_o_to_n[o];
        }

        /**
         * Get the old value for n.
         *
         * @param n New value.
         */
        T get_o(const T n) const {
            return m_translation_n_to_o[n];
        }
    };
}

#endif //HEIPROMAP_TRANSLATION_TABLE_H
