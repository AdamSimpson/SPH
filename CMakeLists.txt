cmake_minimum_required(VERSION 2.8)
project(sph)

set(SRCS geometry.c hash.c fileio.c communication.c fluid.c)

add_executable(sph ${SRCS})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph ${MPI_LIBRARIES})

find_library(M_LIB m)
target_link_libraries(sph ${M_LIB})
