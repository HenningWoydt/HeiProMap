#ifndef HEIDELBERGPROCESSMAPPING_MACROS_H
#define HEIDELBERGPROCESSMAPPING_MACROS_H

#include <iostream>
#include <string>

namespace HeiProMap {
#ifndef ASSERT_ENABLED
#define ASSERT_ENABLED true
#endif

#ifndef HEAVYASSERT_ENABLED
#define HEAVYASSERT_ENABLED false
#endif

#if (ASSERT_ENABLED)
    // Use ASSERT for quick operations like O(1) operations, for other Asserts use HEAVYASSERT
#define ASSERT(condition) if(!(condition)) {std::cerr << "Error in file " << __FILE__ << " in function " << __FUNCTION__ << " at line " << __LINE__ << "!" << std::endl; abort(); } ((void)0)
#else
#define ASSERT(condition) ((void)0)
#endif


#if (HEAVYASSERT_ENABLED)
    // Use HEAVYASSERT for expensive operations like O(n), O(n^2) operations, for faster Asserts use ASSERT
#define HEAVYASSERT(condition) if(!(condition)) {std::cerr << "Error in file " << __FILE__ << " in function " << __FUNCTION__ << " at line " << __LINE__ << "!" << std::endl; abort(); } ((void)0)
#else
#define HEAVYASSERT(condition) ((void)0)
#endif
}

#endif //HEIDELBERGPROCESSMAPPING_MACROS_H
