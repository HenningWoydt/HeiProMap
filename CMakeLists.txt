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

# KaHIP library
set(NOMPI ON CACHE BOOL "For KaHIP: disable all targets that depend on MPI (kaffpaE, ParHIP)")
include_directories(extern/KaHIP)
link_directories(extern/KaHIP)
add_subdirectory(extern/KaHIP)

# if(NOT GUROBI_DIR)
#     message(FATAL_ERROR "GUROBI_DIR not set — please pass it using -DGUROBI_DIR=/opt/gurobi1201/linux64")
# endif()

# message(STATUS "GUROBI_DIR is ${GUROBI_DIR}")
# message(STATUS "GUROBI_HOME is $ENV{GUROBI_HOME}")

# find_package(GUROBI REQUIRED)
# include_directories(${GUROBI_INCLUDE_DIRS})

# message(STATUS "GUROBI_INCLUDE_DIRS is ${GUROBI_INCLUDE_DIRS}")
# message(STATUS "GUROBI_CXX_LIBRARY is ${GUROBI_CXX_LIBRARY}")
# message(STATUS "GUROBI_LIBRARY is ${GUROBI_LIBRARY}")

# Compiler Flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic -Wno-unused-variable -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3 -DNDEBUG -fno-omit-frame-pointer -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_MINSIZEREL "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -pedantic")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# Find all commons source and header files
file(GLOB_RECURSE SOURCES_COMMONS CONFIGURE_DEPENDS "src/commons/*.cpp")
file(GLOB_RECURSE HEADERS_COMMONS CONFIGURE_DEPENDS "src/commons/*.h")

# Find all serial source and header files
file(GLOB_RECURSE SOURCES_SERIAL CONFIGURE_DEPENDS "src/serial/*.cpp")
file(GLOB_RECURSE HEADERS_SERIAL CONFIGURE_DEPENDS "src/serial/*.h")

# Find all parallel source and header files
file(GLOB_RECURSE SOURCES_PARALLEL CONFIGURE_DEPENDS "src/parallel/*.cpp")
file(GLOB_RECURSE HEADERS_PARALLEL CONFIGURE_DEPENDS "src/parallel/*.h")

# Find all maxflow files
file(GLOB_RECURSE MAXFLOW_SOURCES CONFIGURE_DEPENDS "extern/maxflow-v3.04.src/*.cpp")
file(GLOB_RECURSE MAXFLOW_HEADERS CONFIGURE_DEPENDS "extern/maxflow-v3.04.src/*.h")

# add_executable(HeiProMap application/HeiProMap.cpp ${HEADERS_COMMONS} ${SOURCES_COMMONS} ${HEADERS_SERIAL} ${SOURCES_SERIAL} ${MAXFLOW_SOURCES} ${MAXFLOW_HEADERS})
# add_executable(MtHeiProMap application/MtHeiProMap.cpp ${HEADERS_COMMONS} ${SOURCES_COMMONS} ${HEADERS_PARALLEL} ${SOURCES_PARALLEL})
add_executable(DeepHeiProMap application/DeepHeiProMap.cpp ${HEADERS_COMMONS} ${SOURCES_COMMONS} ${HEADERS_SERIAL} ${SOURCES_SERIAL} ${MAXFLOW_SOURCES} ${MAXFLOW_HEADERS})

# Link OpenMP library
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    # target_link_libraries(MtHeiProMap OpenMP::OpenMP_CXX)
    target_link_libraries(DeepHeiProMap OpenMP::OpenMP_CXX)
endif ()

# Link KaHIP library
# target_link_libraries(HeiProMap kahip)
# target_link_libraries(MtHeiProMap kahip)
target_link_libraries(DeepHeiProMap kahip)

# target_link_libraries(HeiProMap ${GUROBI_CXX_LIBRARY})
# target_link_libraries(HeiProMap ${GUROBI_LIBRARY})
