cmake_minimum_required(VERSION 3.22)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(HeiProMap)
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
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -fno-exceptions")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g3 -DNDEBUG -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g3 -DNDEBUG -fno-inline -fno-omit-frame-pointer -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -fno-exceptions")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -fno-exceptions")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -pedantic -fno-exceptions")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# Find all interfaces source and header files
file(GLOB_RECURSE SOURCES_INTERFACE CONFIGURE_DEPENDS "src/interfaces/*.cpp")
file(GLOB_RECURSE HEADERS_INTERFACE CONFIGURE_DEPENDS "src/interfaces/*.h")

# Find all serial source and header files
file(GLOB_RECURSE SOURCES_SERIAL CONFIGURE_DEPENDS "src/serial/*.cpp")
file(GLOB_RECURSE HEADERS_SERIAL CONFIGURE_DEPENDS "src/serial/*.h")

# Find all parallel source and header files
file(GLOB_RECURSE SOURCES_PARALLEL CONFIGURE_DEPENDS "src/parallel/*.cpp")
file(GLOB_RECURSE HEADERS_PARALLEL CONFIGURE_DEPENDS "src/parallel/*.h")

# Find all serial test files
file(GLOB_RECURSE Test_SOURCES_SERIAL CONFIGURE_DEPENDS "tests/serial/*.cpp")
file(GLOB_RECURSE Test_HEADERS_SERIAL CONFIGURE_DEPENDS "tests/serial/*.h")

# Find all serial test files
file(GLOB_RECURSE Test_SOURCES_PARALLEL CONFIGURE_DEPENDS "tests/parallel/*.cpp")
file(GLOB_RECURSE Test_HEADERS_PARALLEL CONFIGURE_DEPENDS "tests/parallel/*.h")

add_executable(HeiProMap application/HeiProMap.cpp ${HEADERS_INTERFACE} ${SOURCES_INTERFACE} ${HEADERS_SERIAL} ${SOURCES_SERIAL})
# add_executable(MtHeiProMap application/MtHeiProMap.cpp ${HEADERS_INTERFACE} ${SOURCES_INTERFACE} ${HEADERS_SERIAL} ${SOURCES_SERIAL} ${HEADERS_PARALLEL} ${SOURCES_PARALLEL})
# add_executable(MtHeiProMapRefinement application/MtHeiProMapRefinement.cpp ${HEADERS_INTERFACE} ${SOURCES_INTERFACE} ${HEADERS_SERIAL} ${SOURCES_SERIAL} ${HEADERS_PARALLEL} ${SOURCES_PARALLEL})

# Link OpenMP library
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    # target_link_libraries(MtHeiProMap OpenMP::OpenMP_CXX)
    # target_link_libraries(MtHeiProMapRefinement OpenMP::OpenMP_CXX)
endif ()

# Link KaHIP library
target_link_libraries(HeiProMap kahip_static)
# target_link_libraries(MtHeiProMap kahip)
# target_link_libraries(MtHeiProMapRefinement kahip)

add_subdirectory(external_tools/googletest)
include(gmock)
enable_testing()

# Test executable
add_executable(HeiProMap_gtest ${HEADERS_INTERFACE} ${SOURCES_INTERFACE} ${HEADERS_SERIAL} ${SOURCES_SERIAL} ${Test_HEADERS_SERIAL} ${Test_SOURCES_SERIAL})
# add_executable(MtHeiProMap_gtest ${HEADERS_INTERFACE} ${SOURCES_INTERFACE} ${HEADERS_SERIAL} ${SOURCES_SERIAL} ${Test_HEADERS_SERIAL} ${Test_SOURCES_SERIAL} ${Test_HEADERS_PARALLEL} ${Test_SOURCES_PARALLEL})

# Link OpenMP library
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    # target_link_libraries(MtHeiProMap_gtest OpenMP::OpenMP_CXX)
endif ()

# Link to Google Test
target_link_libraries(HeiProMap_gtest gtest gtest_main)
# target_link_libraries(MtHeiProMap_gtest gtest gtest_main)
