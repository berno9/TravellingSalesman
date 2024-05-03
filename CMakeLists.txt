cmake_minimum_required(VERSION 3.26.4)
project(TravellingSalesman)

set(CMAKE_CXX_STANDARD 17)

add_executable(TravellingSalesman src/main.cpp
        src/headers/Graph.h
        src/headers/Script.h
        src/classes/Script.cpp
        src/headers/TSPSolver.h
        src/classes/TSPSolver.cpp
        src/headers/Menu.h
        src/classes/Menu.cpp)
