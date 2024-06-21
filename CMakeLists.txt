cmake_minimum_required(VERSION 3.10)
project(segrec)

set(CMAKE_CXX_STANDARD 20)



SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -O3 -lz -lm -Wall -march=native")

set(CMAKE_BUILD_TYPE Release)



add_executable(segrec main.cpp define.h genespeciestreeutil.h hashtable.h newicklex.h node.h SegmentalReconciler.h treeinfo.h treeiterator.h util.h genespeciestreeutil.cpp newicklex.cpp node.cpp SegmentalReconciler.cpp treeinfo.cpp treeiterator.cpp)

