cmake_minimum_required(VERSION 3.26.4)
project(TravellingSalesman)

set(CMAKE_CXX_STANDARD 17)

add_executable(TravellingSalesman src/main.cpp
        src/headers/Graph.h
        src/headers/Script.h
        src/classes/Script.cpp
        src/headers/ToyEdge.h
        src/classes/ToyEdge.cpp)
