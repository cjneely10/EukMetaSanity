cmake_minimum_required(VERSION 3.17)
project(MergeRegions)

set(CMAKE_CXX_STANDARD 17)

add_executable(MergeRegions main.cpp Parser.cpp Parser.h Gene.cpp Gene.h Record.h exceptions.h Region.h)