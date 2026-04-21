cmake_minimum_required(VERSION 3.16)

project(HeiProMap LANGUAGES CXX)

set(PROJECT_VENDOR "Henning Woydt")
set(PROJECT_CONTACT "henning.woydt@informatik.uni-heidelberg.de")
set(PROJECT_URL "https://github.com/HenningWoydt/HeiProMap")
set(PROJECT_DESCRIPTION "HeiProMap: Heidelberg Process Mapping")
set(PROJECT_VERSION "1.0.0")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------
# Language / compiler setup
# ---------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3 -DNDEBUG -fno-omit-frame-pointer -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -march=x86-64 -mtune=generic -fno-omit-frame-pointer -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")

# Enable IPO only where it makes sense
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL TRUE)

set(ENABLE_PROFILER ON CACHE BOOL "Enable profiler")
set(ENABLE_ASSERTS OFF CACHE BOOL "Enable assertions")

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
message(STATUS "ENABLE_PROFILER     : ${ENABLE_PROFILER}")
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
# Executable
# ---------------------------
add_executable(HeiProMap
        main.cpp
        ${HEADERS}
        ${SOURCES}
)

target_include_directories(HeiProMap PRIVATE ${HEIPROMAP_PRIVATE_INCLUDES})

target_link_libraries(HeiProMap PUBLIC
        ${KAHIP_LIBRARY}
        OpenMP::OpenMP_CXX
        TBB::tbb
        m
)

set_target_properties(HeiProMap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib"
)

target_link_options(HeiProMap PRIVATE
        "-Wl,--disable-new-dtags"
)

if(ENABLE_PROFILER)
    target_compile_definitions(HeiProMap PRIVATE ENABLE_PROFILER=1)
else()
    target_compile_definitions(HeiProMap PRIVATE ENABLE_PROFILER=0)
endif()

if(ENABLE_ASSERTS)
    target_compile_definitions(HeiProMap PRIVATE ASSERT_ENABLED=true)
endif()