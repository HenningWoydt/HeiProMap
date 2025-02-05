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

#include "../../definitions.h"
#include "../../macros.h"

namespace HeiProMap {
    class TranslationTable {
        std::unordered_map<u64, u64> m_translation_o_to_n;
        std::unordered_map<u64, u64> m_translation_n_to_o;

    public:
        /**
         * Default constructor.
         */
        TranslationTable() = default;

        /**
         * Initializes the translation table with the Identity mapping.
         */
        explicit TranslationTable(u64 n) {
            for (u64 u = 0; u < n; ++u) {
                add(u, u);
            }
        }

        /**
         * Adds a translation from o to n and from n to o.
         *
         * @param o Old value.
         * @param n New value.
         */
        void add(u64 o, u64 n) {
            ASSERT(m_translation_o_to_n.find(o) == m_translation_o_to_n.end());
            ASSERT(m_translation_n_to_o.find(n) == m_translation_n_to_o.end());

            m_translation_o_to_n[o] = n;
            m_translation_n_to_o[n] = o;
        }

        /**
         * Get the new value for o.
         *
         * @param o Old value
         */
        u64 get_n(u64 o) const {
            ASSERT(m_translation_o_to_n.find(o) != m_translation_o_to_n.end());

            return m_translation_o_to_n.at(o);
        }

        /**
         * Get the old value for n.
         *
         * @param n New value.
         */
        u64 get_o(u64 n) const {
            ASSERT(m_translation_n_to_o.find(n) != m_translation_n_to_o.end());

            return m_translation_n_to_o.at(n);
        }

        /**
         * Removes all entries.
         */
        void clear() {
            m_translation_o_to_n.clear();
            m_translation_n_to_o.clear();
        }

        void merge(TranslationTable& tt) {
            m_translation_n_to_o.merge(tt.m_translation_n_to_o);
            m_translation_o_to_n.merge(tt.m_translation_o_to_n);
        }

        void print() {
            std::cout << "Old to New" << std::endl;
            for (const auto& [key, value] : m_translation_o_to_n) {
                std::cout << key << " : " << value << std::endl;
            }
            std::cout << "New to old" << std::endl;
            for (const auto& [key, value] : m_translation_n_to_o) {
                std::cout << key << " : " << value << std::endl;
            }
        }
    };
}

#endif //HEIPROMAP_TRANSLATION_TABLE_H
