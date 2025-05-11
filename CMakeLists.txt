cmake_minimum_required(VERSION 3.31)
project(delaunay)

set(CMAKE_CXX_STANDARD 20)

add_executable(delaunay main.cpp
        frontier.h
        frontier.cpp
        geometry.h
        delaunay.h
        delaunay.cpp)

add_library(triangulate_api SHARED triangulate_api.cpp frontier.h
        frontier.cpp
        geometry.h
        delaunay.h
        delaunay.cpp)

set_target_properties(triangulate_api PROPERTIES
        OUTPUT_NAME "triangulate"
        SUFFIX ".dylib"
)
