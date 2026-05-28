file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS "*.cpp")

add_executable(HeiProMapTests
    ${TEST_SOURCES}
)

# Exclude from default 'all' target if BUILD_TESTING is OFF
if(NOT BUILD_TESTING)
    set_target_properties(HeiProMapTests PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()

# GTest might be provided via FetchContent (targets gtest, etc.)
# or via find_package (target GTest::gtest or GTest::GTest)
if(TARGET GTest::gtest)
    set(GTEST_LIBS GTest::gtest)
elseif(TARGET GTest::GTest)
    set(GTEST_LIBS GTest::GTest)
else()
    set(GTEST_LIBS gtest)
endif()

target_link_libraries(HeiProMapTests
    PRIVATE
        HeiProMapLib
        ${GTEST_LIBS}
)

include(GoogleTest)
gtest_discover_tests(HeiProMapTests)
