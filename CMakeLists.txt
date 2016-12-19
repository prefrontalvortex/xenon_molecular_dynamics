cmake_minimum_required(VERSION 3.5)
project(mol)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -lm -static-libgcc -lgfortran -std=c++11 ")
set(CMAKE_CXX_STANDARD 11)


set(SOURCE_FILES main.cpp rng.cpp rng.h main.c aux.c aux.h xspeedtest.c argparse.c argparse.h )
add_executable(mol ${SOURCE_FILES})
set(SOURCE_FILES main.c aux.c aux.h argparse.c argparse.h)
add_executable(mol ${SOURCE_FILES})
set_target_properties(mol PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin/ )

target_link_libraries(mol m Threads::Threads)
