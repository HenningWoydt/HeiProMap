#include <gtest/gtest.h>
#include "src/utility/translation_table.h"

namespace HeiProMap {

TEST(TranslationTableTest, BasicFunctionality) {
    TranslationTable<int> tt;
    tt.reserve(10, 10);
    
    tt.add(0, 5);
    tt.add(1, 6);
    
    EXPECT_EQ(tt.get_n(0), 5);
    EXPECT_EQ(tt.get_n(1), 6);
    EXPECT_EQ(tt.get_o(5), 0);
    EXPECT_EQ(tt.get_o(6), 1);
}

} // namespace HeiProMap
