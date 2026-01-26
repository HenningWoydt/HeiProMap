cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

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

# Fail early if something is missing
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

file(GLOB_RECURSE MAXFLOW_SOURCES CONFIGURE_DEPENDS "extern/maxflow/*.cpp")
file(GLOB_RECURSE MAXFLOW_HEADERS CONFIGURE_DEPENDS "extern/maxflow/*.h")

foreach(maxflow_src ${MAXFLOW_SOURCES})
    set_source_files_properties(${maxflow_src} PROPERTIES COMPILE_FLAGS "-w")
endforeach()

# ---------------------------
# Executable
# ---------------------------
add_executable(HeiProMap
        main.cpp
        ${HEADERS}
        ${SOURCES}
        ${MAXFLOW_SOURCES}
        ${MAXFLOW_HEADERS}
)

target_include_directories(HeiProMap PRIVATE
        ${KAHIP_INCLUDE_DIR}
        ${METIS_INCLUDE_DIR}
        ${GKLIB_INCLUDE_DIR}
)

# IMPORTANT: static link order -> METIS then GKlib
target_link_libraries(HeiProMap PRIVATE
        ${KAHIP_LIBRARY}
        ${METIS_LIBRARY}
        ${GKLIB_LIBRARY}
        m
)

# If you rely on shared libs from those local folders at runtime
set_target_properties(HeiProMap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib;${METIS_ROOT}/lib;${GKLIB_ROOT}/lib"
)

# ---------------------------
# Shared library
# ---------------------------
add_library(libheipromap OBJECT
        ${HEADERS}
        ${SOURCES}
        ${MAXFLOW_SOURCES}
        ${MAXFLOW_HEADERS}
)

target_include_directories(libheipromap PRIVATE
        ${KAHIP_INCLUDE_DIR}
        ${METIS_INCLUDE_DIR}
        ${GKLIB_INCLUDE_DIR}
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
        m
)

set_target_properties(heipromap PROPERTIES
        BUILD_RPATH "${KAHIP_ROOT}/lib;${METIS_ROOT}/lib;${GKLIB_ROOT}/lib"
)

install(TARGETS heipromap DESTINATION lib)
