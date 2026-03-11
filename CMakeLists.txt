cmake_minimum_required(VERSION 3.16)

project(HeiProMap LANGUAGES CXX)

set(PROJECT_VENDOR "Henning Woydt")
set(PROJECT_CONTACT "henning.woydt@informatik.uni-heidelberg.de")
set(PROJECT_URL "https://github.com/HenningWoydt/HeiProMap")
set(PROJECT_DESCRIPTION "HeiProMap: Heidelberg Process Mapping")
set(PROJECT_VERSION "1.0.0")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler Flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3 -DNDEBUG -fno-omit-frame-pointer -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -pedantic")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# ---------------------------
# KaHIP
# ---------------------------
set(KAHIP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/kahip")

find_library(KAHIP_LIBRARY
        NAMES kahip libkahip
        PATHS ${KAHIP_ROOT}/lib
        NO_DEFAULT_PATH
)
find_path(KAHIP_INCLUDE_DIR
        NAMES kaHIP_interface.h
        PATHS ${KAHIP_ROOT}/include
        NO_DEFAULT_PATH
)

message(STATUS "Found KaHIP library: ${KAHIP_LIBRARY}")
message(STATUS "Found KaHIP include: ${KAHIP_INCLUDE_DIR}")

# ---------------------------
# METIS
# ---------------------------
set(METIS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/metis")

find_library(METIS_LIBRARY
        NAMES metis libmetis
        PATHS ${METIS_ROOT}/lib
        NO_DEFAULT_PATH
)
find_path(METIS_INCLUDE_DIR
        NAMES metis.h
        PATHS ${METIS_ROOT}/include
        NO_DEFAULT_PATH
)

message(STATUS "Found METIS library: ${METIS_LIBRARY}")
message(STATUS "Found METIS include: ${METIS_INCLUDE_DIR}")

# ---------------------------
# GKlib (needed when METIS is static libmetis.a)
# ---------------------------
set(GKLIB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/gklib")

find_library(GKLIB_LIBRARY
        NAMES GKlib gklib libGKlib libgklib
        PATHS ${GKLIB_ROOT}/lib
        NO_DEFAULT_PATH
)
find_path(GKLIB_INCLUDE_DIR
        NAMES gk_lib.h GKlib.h
        PATHS ${GKLIB_ROOT}/include
        NO_DEFAULT_PATH
)

message(STATUS "Found GKlib library: ${GKLIB_LIBRARY}")
message(STATUS "Found GKlib include: ${GKLIB_INCLUDE_DIR}")

# ---------------------------
# Mt-KaHyPar (CMake package)
# ---------------------------
set(MTK_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/mt-kahypar")
list(APPEND CMAKE_PREFIX_PATH "${MTK_ROOT}")

find_package(MtKaHyPar CONFIG REQUIRED)

# optional: help runtime loader find the .so when using the local install
if (EXISTS "${MTK_ROOT}/lib64")
    set(MTK_LIBDIR "${MTK_ROOT}/lib64")
else()
    set(MTK_LIBDIR "${MTK_ROOT}/lib")
endif()

message(STATUS "Found Mt-KaHyPar via CMake package: MtKaHyPar::mtkahypar")
message(STATUS "Mt-KaHyPar libdir: ${MTK_LIBDIR}")

# ---------------------------
# Fail early if something is missing
# ---------------------------
if (NOT KAHIP_LIBRARY OR NOT KAHIP_INCLUDE_DIR)
    message(FATAL_ERROR "KaHIP not found under ${KAHIP_ROOT}")
endif()
if (NOT METIS_LIBRARY OR NOT METIS_INCLUDE_DIR)
    message(FATAL_ERROR "METIS not found under ${METIS_ROOT}")
endif()
if (NOT GKLIB_LIBRARY OR NOT GKLIB_INCLUDE_DIR)
    message(FATAL_ERROR "GKlib not found under ${GKLIB_ROOT}")
endif()

# ---------------------------
# Sources
# ---------------------------
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "src/*.h")

file(GLOB_RECURSE MAXFLOW_ALG_bk_HEADERS CONFIGURE_DEPENDS "extern/maxflow_algorithms/bk/*.h")
file(GLOB_RECURSE MAXFLOW_ALG_hi_pr_SOURCES CONFIGURE_DEPENDS "extern/maxflow_algorithms/hi_pr/*.cpp")
list(FILTER MAXFLOW_ALG_hi_pr_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

# ---------------------------
# Executable
# ---------------------------
add_executable(HeiProMap
        main.cpp
        ${HEADERS}
        ${SOURCES}
        ${MAXFLOW_ALG_hi_pr_SOURCES}
        ${MAXFLOW_ALG_bk_HEADERS}
)

target_include_directories(HeiProMap PRIVATE
        ${KAHIP_INCLUDE_DIR}
        ${METIS_INCLUDE_DIR}
        ${GKLIB_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/extern/maxflow_algorithms/bk
        ${CMAKE_CURRENT_SOURCE_DIR}/extern/maxflow_algorithms/hi_pr
)

# Link order matters for static libs: METIS then GKlib
target_link_libraries(HeiProMap PUBLIC
        ${KAHIP_LIBRARY}
        ${METIS_LIBRARY}
        ${GKLIB_LIBRARY}
        MtKaHyPar::mtkahypar
        m
)

set_target_properties(HeiProMap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib;${METIS_ROOT}/lib;${GKLIB_ROOT}/lib;${MTK_LIBDIR}"
)

# ---------------------------
# Shared library
# ---------------------------
add_library(libheipromap OBJECT
        ${HEADERS}
        ${SOURCES}
        ${MAXFLOW_ALG_bk_HEADERS}
        ${MAXFLOW_ALG_hi_pr_SOURCES}
)

target_include_directories(libheipromap PRIVATE
        ${KAHIP_INCLUDE_DIR}
        ${METIS_INCLUDE_DIR}
        ${GKLIB_INCLUDE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/extern/maxflow_algorithms/bk
        ${CMAKE_CURRENT_SOURCE_DIR}/extern/maxflow_algorithms/hi_pr
)

add_library(heipromap SHARED
        lib/libheipromap.cpp
        $<TARGET_OBJECTS:libheipromap>
)

target_include_directories(heipromap
        PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/interface
        PRIVATE ${KAHIP_INCLUDE_DIR} ${METIS_INCLUDE_DIR} ${GKLIB_INCLUDE_DIR}
)

target_link_libraries(heipromap PRIVATE
        ${KAHIP_LIBRARY}
        ${METIS_LIBRARY}
        ${GKLIB_LIBRARY}
        MtKaHyPar::mtkahypar
        m
)

set_target_properties(heipromap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib;${METIS_ROOT}/lib;${GKLIB_ROOT}/lib;${MTK_LIBDIR}"
)


foreach(tgt HeiProMap heipromap)
    target_link_options(${tgt} PRIVATE
            "-Wl,--disable-new-dtags"
            "-Wl,-rpath-link,${MTK_LIBDIR}"
    )
    set_target_properties(${tgt} PROPERTIES
            BUILD_RPATH "${KaHIP_ROOT}/lib;${MTK_LIBDIR}"
    )
endforeach()
