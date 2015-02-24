cmake_minimum_required(VERSION 2.8)
project(sph)

file(GLOB SRC_FILES src/*.c)

add_executable(sph ${SRC_FILES})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(sph ${MPI_LIBRARIES})

find_library(M_LIB m)
target_link_libraries(sph ${M_LIB})

install(TARGETS sph DESTINATION bin)

