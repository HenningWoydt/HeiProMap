cmake_minimum_required(VERSION 3.16)

project(HeiProMap LANGUAGES CXX)

set(PROJECT_VENDOR "Henning Woydt")
set(PROJECT_CONTACT "henning.woydt@informatik.uni-heidelberg.de")
set(PROJECT_URL "https://github.com/HenningWoydt/HeiProMap")
set(PROJECT_DESCRIPTION "HeiProMap: Heidelberg Process Mapping")
set(PROJECT_VERSION "1.0.0")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------
# KaMinPar (shared-memory partitioner, built from extern/KaMinPar)
# Added before HeiProMap's global compiler flags (-fno-exceptions, -fwhole-program, ...) so it builds with its own.
# ---------------------------
set(KAMINPAR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/KaMinPar")
if(NOT EXISTS "${KAMINPAR_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "KaMinPar not found under ${KAMINPAR_SOURCE_DIR}; run build.sh to fetch it")
endif()
set(KAMINPAR_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_APPS OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_WITH_SPARSEHASH OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_WITH_CCACHE OFF CACHE BOOL "" FORCE)
set(KAMINPAR_BUILD_WITH_MTUNE_NATIVE OFF CACHE BOOL "" FORCE)
set(KAMINPAR_ENABLE_TIMERS OFF CACHE BOOL "" FORCE)
set(KAMINPAR_64BIT_EDGE_IDS ON CACHE BOOL "" FORCE)
set(INSTALL_KAMINPAR OFF CACHE BOOL "" FORCE)
add_subdirectory("${KAMINPAR_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/KaMinPar" EXCLUDE_FROM_ALL)

# ---------------------------
# Options
# ---------------------------
option(ENABLE_ASSERTS "Enable assertions" OFF)

add_definitions(-DHEIPROMAP_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

if(ENABLE_ASSERTS)
    add_definitions(-DENABLE_ASSERTS=1)
    add_definitions(-DHEAVYASSERT_ENABLED=1)
endif()

# ---------------------------
# Language / compiler setup
# ---------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(ENABLE_EXCEPTIONS "Enable C++ exceptions" OFF)

if(NOT ENABLE_EXCEPTIONS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
endif()

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3 -DNDEBUG -fno-omit-frame-pointer -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -rdynamic -march=x86-64 -mtune=generic -fno-omit-frame-pointer -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")

# Enable IPO only where it makes sense
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL TRUE)

# ---------------------------
# OpenMP
# ---------------------------
find_package(OpenMP REQUIRED)

# ---------------------------
# TBB
# ---------------------------
find_package(TBB REQUIRED)

# ---------------------------
# KaHIP selection by config
# ---------------------------
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    # set(KAHIP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/kahip-debug")
    set(KAHIP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/kahip-release")
else()
    set(KAHIP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/kahip-release")
endif()

find_library(KAHIP_LIBRARY
        NAMES kahip libkahip
        PATHS "${KAHIP_ROOT}/lib"
        NO_DEFAULT_PATH
)

find_path(KAHIP_INCLUDE_DIR
        NAMES kaHIP_interface.h
        PATHS "${KAHIP_ROOT}/include"
        NO_DEFAULT_PATH
)

message(STATUS "Build type          : ${CMAKE_BUILD_TYPE}")
message(STATUS "ENABLE_ASSERTS      : ${ENABLE_ASSERTS}")
message(STATUS "Selected KaHIP root : ${KAHIP_ROOT}")
message(STATUS "Found KaHIP library : ${KAHIP_LIBRARY}")
message(STATUS "Found KaHIP include : ${KAHIP_INCLUDE_DIR}")
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP        : TRUE")
    message(STATUS "OpenMP version      : ${OpenMP_CXX_VERSION}")
else()
    message(STATUS "Found OpenMP        : FALSE")
endif()
if(TBB_FOUND)
    message(STATUS "Found TBB           : TRUE")
else()
    message(STATUS "Found TBB           : FALSE")
endif()
message(STATUS "TBB version         : ${TBB_VERSION}")
if(TARGET TBB::tbb)
    get_target_property(_TBB_INCLUDE TBB::tbb INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "TBB include         : ${_TBB_INCLUDE}")
endif()

if(NOT KAHIP_LIBRARY OR NOT KAHIP_INCLUDE_DIR)
    message(FATAL_ERROR "KaHIP not found under ${KAHIP_ROOT}")
endif()

# ---------------------------
# Sources
# ---------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "src/*.h")

# ---------------------------
# Common include dirs
# ---------------------------
set(HEIPROMAP_PRIVATE_INCLUDES
        ${KAHIP_INCLUDE_DIR})

# ---------------------------
# Library
# ---------------------------
add_library(HeiProMapLib INTERFACE)
target_include_directories(HeiProMapLib INTERFACE ${PROJECT_SOURCE_DIR})
target_include_directories(HeiProMapLib INTERFACE ${HEIPROMAP_PRIVATE_INCLUDES})
target_link_libraries(HeiProMapLib INTERFACE
        ${KAHIP_LIBRARY}
        OpenMP::OpenMP_CXX
        TBB::tbb
        KaMinPar::KaMinPar
        m
)

# ---------------------------
# Executables
# ---------------------------
add_executable(HeiProMap app/HeiProMap.cpp)
target_link_libraries(HeiProMap PRIVATE HeiProMapLib)
set_target_properties(HeiProMap PROPERTIES BUILD_RPATH "${KAHIP_ROOT}/lib")
target_link_options(HeiProMap PRIVATE "-Wl,--disable-new-dtags")

add_executable(Dyn-HeiProMap app/Dyn-HeiProMap.cpp)
target_link_libraries(Dyn-HeiProMap PRIVATE HeiProMapLib)
set_target_properties(Dyn-HeiProMap PROPERTIES BUILD_RPATH "${KAHIP_ROOT}/lib")
target_link_options(Dyn-HeiProMap PRIVATE "-Wl,--disable-new-dtags")

add_executable(HeiPa app/HeiPa.cpp)
target_link_libraries(HeiPa PRIVATE HeiProMapLib)
set_target_properties(HeiPa PROPERTIES BUILD_RPATH "${KAHIP_ROOT}/lib")
target_link_options(HeiPa PRIVATE "-Wl,--disable-new-dtags")

add_executable(Deep-HeiProMap app/Deep-HeiProMap.cpp)
target_link_libraries(Deep-HeiProMap PRIVATE HeiProMapLib)
set_target_properties(Deep-HeiProMap PROPERTIES BUILD_RPATH "${KAHIP_ROOT}/lib")
target_link_options(Deep-HeiProMap PRIVATE "-Wl,--disable-new-dtags")

add_executable(Waver-HeiProMap app/Waver-HeiProMap.cpp)
target_link_libraries(Waver-HeiProMap PRIVATE HeiProMapLib)
set_target_properties(Waver-HeiProMap PROPERTIES BUILD_RPATH "${KAHIP_ROOT}/lib")
target_link_options(Waver-HeiProMap PRIVATE "-Wl,--disable-new-dtags")
