cmake_minimum_required(VERSION 3.28)
project(C__)

set(CMAKE_CXX_STANDARD 14)

add_executable(C__
        EKF_Test.cpp
        INCLUDE/Matrix.h
        SRC/Matrix.cpp
        EKF_GEOBS.cpp
        INCLUDE/Mjday.h
        INCLUDE/Mjday.h
        SRC/Mjday.cpp
        SRC/Mjday.cpp
)
