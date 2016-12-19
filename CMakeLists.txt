cmake_minimum_required(VERSION 3.6)
project(mol)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp rng.cpp rng.h main.c aux.c aux.h xspeedtest.c argparse.c argparse.h hep_original.cpp)
add_executable(mol ${SOURCE_FILES})