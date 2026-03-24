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

set(ENABLE_PROFILER OFF CACHE BOOL "Enable profiler")
message(STATUS "ENABLE_PROFILER    : ${ENABLE_PROFILER}")

# ---------------------------
# OpenMP
# ---------------------------
find_package(OpenMP REQUIRED)

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
message(STATUS "Selected KaHIP root : ${KAHIP_ROOT}")
message(STATUS "Found KaHIP library : ${KAHIP_LIBRARY}")
message(STATUS "Found KaHIP include : ${KAHIP_INCLUDE_DIR}")
message(STATUS "Found OpenMP        : ${OpenMP_CXX_FOUND}")

if(NOT KAHIP_LIBRARY OR NOT KAHIP_INCLUDE_DIR)
    message(FATAL_ERROR "KaHIP not found under ${KAHIP_ROOT}")
endif()

# ---------------------------
# Sources
# ---------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "src/*.h")

file(GLOB_RECURSE MAXFLOW_ALG_bk_HEADERS CONFIGURE_DEPENDS "extern/maxflow_algorithms/nbk/*.h")

# ---------------------------
# Common include dirs
# ---------------------------
set(HEIPROMAP_PRIVATE_INCLUDES
        ${KAHIP_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/extern/maxflow_algorithms/nbk
)

# ---------------------------
# Executable
# ---------------------------
add_executable(HeiProMap
        main.cpp
        ${HEADERS}
        ${SOURCES}
        ${MAXFLOW_ALG_bk_HEADERS}
)

target_include_directories(HeiProMap PRIVATE ${HEIPROMAP_PRIVATE_INCLUDES})

target_link_libraries(HeiProMap PUBLIC
        ${KAHIP_LIBRARY}
        OpenMP::OpenMP_CXX
        m
)

set_target_properties(HeiProMap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib"
)

target_link_options(HeiProMap PRIVATE
        "-Wl,--disable-new-dtags"
)

# ---------------------------
# Object library
# ---------------------------
add_library(libheipromap OBJECT
        ${HEADERS}
        ${SOURCES}
        ${MAXFLOW_ALG_bk_HEADERS}
)

target_include_directories(libheipromap PRIVATE ${HEIPROMAP_PRIVATE_INCLUDES})

target_compile_definitions(libheipromap PRIVATE
        ENABLE_PROFILER=${ENABLE_PROFILER}
)

target_link_libraries(libheipromap PUBLIC
        OpenMP::OpenMP_CXX
)

# ---------------------------
# Shared library
# ---------------------------
add_library(heipromap SHARED
        lib/libheipromap.cpp
        $<TARGET_OBJECTS:libheipromap>
)

target_include_directories(heipromap
        PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/interface
        PRIVATE ${KAHIP_INCLUDE_DIR}
)

target_compile_definitions(HeiProMap PRIVATE
        ENABLE_PROFILER=${ENABLE_PROFILER}
)

target_link_libraries(heipromap PRIVATE
        ${KAHIP_LIBRARY}
        OpenMP::OpenMP_CXX
)

set_target_properties(heipromap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib"
)

target_link_options(heipromap PRIVATE
        "-Wl,--disable-new-dtags"
)