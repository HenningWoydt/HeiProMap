#include <gtest/gtest.h>
#include "src/utility/flat_translation_table.h"

namespace HeiProMap {

TEST(FlatTranslationTableTest, BasicTranslation) {
    FlatTranslationTable<int> tt(2);
    tt.add(0, 10);
    tt.add(1, 20);

    EXPECT_EQ(tt.get_n(0), 10);
    EXPECT_EQ(tt.get_o(10), 0);
}

TEST(FlatTranslationTableTest, Resize) {
    FlatTranslationTable<int> tt(2);
    tt.add(0, 10);

    tt.resize(4);
    tt.add(2, 30);
    
    EXPECT_EQ(tt.get_n(2), 30);
}

} // namespace HeiProMap
