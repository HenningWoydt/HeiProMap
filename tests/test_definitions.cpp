#include <gtest/gtest.h>
#include "src/definitions.h"

namespace HeiProMap {

TEST(DefinitionsTest, MoveStruct) {
    Move m(1, 2, 3);
    EXPECT_EQ(m.u, 1);
    EXPECT_EQ(m.u_id, 2);
    EXPECT_EQ(m.to_move_id, 3);
}

} // namespace HeiProMap
