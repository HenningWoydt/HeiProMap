#include <gtest/gtest.h>
#include "src/utility/JSON_utils.h"

namespace HeiProMap {

TEST(JSONUtilsTest, TestIntegerTypes) {
    EXPECT_EQ(to_JSON_value((u8)5), "5");
    EXPECT_EQ(to_JSON_value((s32)-10), "-10");
}

TEST(JSONUtilsTest, TestFloatingPointTypes) {
    EXPECT_EQ(to_JSON_value(3.14f), "3.140000");
    EXPECT_EQ(to_JSON_value(2.718), "2.718000");
}

TEST(JSONUtilsTest, TestString) {
    EXPECT_EQ(to_JSON_value(std::string("hello")), "\"hello\"");
}

TEST(JSONUtilsTest, TestJSONString) {
    JSONString js{"{\"key\":\"value\"}"};
    EXPECT_EQ(to_JSON_value(js), "{\"key\":\"value\"}");
}

TEST(JSONUtilsTest, TestVector) {
    std::vector<int> v = {1, 2, 3};
    EXPECT_EQ(to_JSON_value(v), "[1, 2, 3]");
    
    std::vector<std::string> vs = {"a", "b"};
    EXPECT_EQ(to_JSON_value(vs), "[\"a\", \"b\"]");
    
    std::vector<int> empty_v;
    EXPECT_EQ(to_JSON_value(empty_v), "[]");
}

TEST(JSONUtilsTest, TestMap) {
    std::map<std::string, int> m = {{"one", 1}, {"two", 2}};
    EXPECT_EQ(to_JSON_value(m), "{\"one\" : 1, \"two\" : 2}");
    
    std::map<int, std::string> m2 = {{1, "one"}, {2, "two"}};
    EXPECT_EQ(to_JSON_value(m2), "{1 : \"one\", 2 : \"two\"}");

    std::map<int, int> empty_m;
    EXPECT_EQ(to_JSON_value(empty_m), "{}");
}

} // namespace HeiProMap
