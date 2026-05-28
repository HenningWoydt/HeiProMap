#include <gtest/gtest.h>
#include "src/utility/utils.h"

namespace HeiProMap {

TEST(UtilsTest, Split) {
    std::string s = "a b c";
    std::vector<std::string> expected = {"a", "b", "c"};
    EXPECT_EQ(split(s, ' '), expected);
}

TEST(UtilsTest, ConvertTo) {
    EXPECT_EQ(convert_to<int>("123"), 123);
    EXPECT_FLOAT_EQ(convert_to<float>("3.14"), 3.14f);
}

TEST(UtilsTest, VectorMath) {
    std::vector<int> v = {1, 2, 3, 4};
    EXPECT_EQ(prod<long long>(v), 24);
    EXPECT_EQ(sum<int>(v), 10);
    EXPECT_EQ(max(v), 4);
    EXPECT_EQ(min(v), 1);
    EXPECT_FLOAT_EQ(avg(v), 2.5);
    EXPECT_EQ(argmin(v), 0);
    EXPECT_EQ(argmax(v), 3);
}

TEST(UtilsTest, Duplicates) {
    std::vector<int> v1 = {1, 2, 3, 4};
    std::vector<int> v2 = {1, 2, 2, 4};
    EXPECT_TRUE(no_duplicates(v1));
    EXPECT_FALSE(no_duplicates(v2));
    EXPECT_TRUE(no_duplicates_sorted(v1));
    EXPECT_FALSE(no_duplicates_sorted(v2));
}

TEST(UtilsTest, Trim) {
    EXPECT_EQ(trim("  hello  "), "hello");
    EXPECT_EQ(trim("hello  "), "hello");
    EXPECT_EQ(trim("  hello"), "hello");
    EXPECT_EQ(trim("hello"), "hello");
}

} // namespace HeiProMap
