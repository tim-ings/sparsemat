cmake_minimum_required(VERSION 3.14)
project(sparsemat C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -fopenmp")
set(CMAKE_C_STANDARD 99)
#find_package(OPENMP REQUIRED)

add_executable(sparsemat main.c matcsr.c matcsr.h ll_float.c ll_float.h matcoo.c matcoo.h matcsc.c matcsc.h list.c list.h)


message(STATUS "Checking OpenMP")
find_package(OpenMP)
IF(OPENMP_FOUND)
    message("Found OpenMP!")
    # add flags for OpenMP
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
    message("Missed OpenMP!")
ENDIF()