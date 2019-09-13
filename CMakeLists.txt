cmake_minimum_required(VERSION 3.14)
project(sparsemat C)

set(CMAKE_C_STANDARD 99)

add_executable(sparsemat main.c matcsr.c matcsr.h ll_float.c ll_float.h matcoo.c matcoo.h matcsc.c matcsc.h)