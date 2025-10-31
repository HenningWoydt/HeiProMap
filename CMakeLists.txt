cmake_minimum_required(VERSION 3.16)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(HeiProMap)
set(PROJECT_VENDOR "Henning Woydt")
set(PROJECT_CONTACT "henning.woydt@informatik.uni-heidelberg.de")
set(PROJECT_URL "https://github.com/HenningWoydt/HeiProMap")
set(PROJECT_DESCRIPTION "HeiProMap: Heidelberg Process Mapping")
set(PROJECT_VERSION "1.0.0")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_language(CXX)
set(CXX ON)

# --- KaHIP library ---
set(KaHIP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/extern/local/kahip")
find_library(KAHIP_LIBRARY
        NAMES kahip libkahip
        PATHS ${KaHIP_ROOT}/lib
        NO_DEFAULT_PATH
)
find_path(KAHIP_INCLUDE_DIR
        NAMES kaHIP_interface.h
        PATHS ${KaHIP_ROOT}/include
        NO_DEFAULT_PATH
)

message(STATUS "Found KaHIP: ${KAHIP_LIBRARY}")
include_directories(${KAHIP_INCLUDE_DIR})

# Compiler Flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3 -DNDEBUG -fno-omit-frame-pointer -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -pedantic")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# Find all serial source and header files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS "src/*.h")

# Find all maxflow files
file(GLOB_RECURSE MAXFLOW_SOURCES CONFIGURE_DEPENDS "extern/maxflow/*.cpp")
file(GLOB_RECURSE MAXFLOW_HEADERS CONFIGURE_DEPENDS "extern/maxflow/*.h")
foreach(maxflow_src ${MAXFLOW_SOURCES})
    set_source_files_properties(${maxflow_src} PROPERTIES
            COMPILE_FLAGS "-w"              # GCC/Clang: disable all warnings
            # or: COMPILE_OPTIONS -w        # equivalent, more modern
    )
endforeach()

add_executable(HeiProMap main.cpp ${HEADERS} ${SOURCES} ${MAXFLOW_SOURCES} ${MAXFLOW_HEADERS})

# Link KaHIP library
target_link_libraries(HeiProMap ${KAHIP_LIBRARY})

# Add the library subdirectory
add_library(libheipromap OBJECT ${HEADERS} ${SOURCES} ${MAXFLOW_SOURCES} ${MAXFLOW_HEADERS})
add_library(heipromap SHARED lib/libheipromap.cpp $<TARGET_OBJECTS:libheipromap>)
target_include_directories(heipromap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/interface)
target_link_libraries(heipromap ${KAHIP_LIBRARY}) # Link KaHIP library
install(TARGETS heipromap DESTINATION lib)

