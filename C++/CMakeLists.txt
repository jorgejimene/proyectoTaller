cmake_minimum_required(VERSION 3.28)
project(C__)

set(CMAKE_CXX_STANDARD 14)

add_executable(C__
        EKF_Test.cpp
        INCLUDE/Matrix.h
        SRC/Matrix.cpp
        EKF_GEOBS.cpp
        INCLUDE/Mjday.h
        SRC/Mjday.cpp
        INCLUDE/R_x.h
        SRC/R_x.cpp
        INCLUDE/R_y.h
        SRC/R_y.cpp
        INCLUDE/R_z.h
        SRC/R_z.cpp
        INCLUDE/TimeUpdate.h
        SRC/TimeUpdate.cpp
        INCLUDE/Legendre.h
        SRC/Legendre.cpp
        INCLUDE/NutAngles.h
        SRC/NutAngles.cpp
        SRC/SAT_Const.cpp
        INCLUDE/SAT_Const.h
        INCLUDE/Position.h
        SRC/Position.cpp
        SRC/Unit.cpp
        INCLUDE/Unit.h
        SRC/Unit.cpp
)
