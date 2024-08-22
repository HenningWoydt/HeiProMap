cmake_minimum_required(VERSION 3.22)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(SerialProcessMapping)
set(PROJECT_VENDOR "Henning Woydt")
set(PROJECT_CONTACT "henning.woydt@informatik.uni-heidelberg.de")
set(PROJECT_URL "https://github.com/HenningWoydt/HeiProMap")
set(PROJECT_DESCRIPTION "HeiProMap: Heidelberg Process Mapping")
set(PROJECT_VERSION "1.0.0")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# KaHIP library
set(NOMPI ON CACHE BOOL "For KaHIP: disable all targets that depend on MPI (kaffpaE, ParHIP)")
include_directories(extern/KaHIP)
link_directories(extern/KaHIP)
add_subdirectory(extern/KaHIP)

# Compiler Flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g3 -DNDEBUG -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall -Wextra -pedantic")

# Find all source and header files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "src/*.h")

# Find all test files
file(GLOB_RECURSE Test_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")
file(GLOB_RECURSE Test_HEADERS CONFIGURE_DEPENDS "tests/*.h")

add_executable(HeiProMap application/HeiProMap.cpp ${HEADERS} ${SOURCES})
add_executable(MtHeiProMap application/MtHeiProMap.cpp ${HEADERS} ${SOURCES})
add_executable(MtHeiProMapRefinement application/MtHeiProMapRefinement.cpp ${HEADERS} ${SOURCES})

# Link OpenMP library
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    target_link_libraries(MtHeiProMap OpenMP::OpenMP_CXX)
    target_link_libraries(MtHeiProMapRefinement OpenMP::OpenMP_CXX)
endif ()

# Link KaHIP library
target_link_libraries(HeiProMap kahip)
target_link_libraries(MtHeiProMap kahip)
target_link_libraries(MtHeiProMapRefinement kahip)

add_subdirectory(external_tools/googletest)
include(gmock)
enable_testing()

# Test executable
add_executable(HeiProMap_gtest ${HEADERS} ${SOURCES} ${Test_HEADERS} ${Test_SOURCES})

# Link to Google Test
target_link_libraries(HeiProMap_gtest gtest gtest_main)
