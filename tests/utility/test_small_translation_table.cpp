#include <gtest/gtest.h>
#include "src/utility/small_translation_table.h"

namespace HeiProMap {

TEST(SmallTranslationTableTest, BasicTranslation) {
    SmallTranslationTable<int> tt;
    tt.add(10, 100);
    tt.add(20, 200);

    EXPECT_EQ(tt.get_n(10), 100);
    EXPECT_EQ(tt.get_o(200), 20);
}

TEST(SmallTranslationTableTest, Clear) {
    SmallTranslationTable<int> tt;
    tt.add(10, 100);
    
    tt.clear();
    
    // Using .at() would throw, which we can't test without exceptions.
    // Instead, we just ensure it can be reused.
    tt.add(30, 300);
    EXPECT_EQ(tt.get_n(30), 300);
}

} // namespace HeiProMap
