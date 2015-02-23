cmake_minimum_required(VERSION 2.8)
project(sph)

set(SRCS geometry.c hash.c fileio.c communication.c fluid.c)

find_library(M_LIB m)
target_link_libraries(sph ${M_LIB})

add_executable(sph ${SRCS})
