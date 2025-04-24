cmake_minimum_required(VERSION 3.28)
project(C__)

set(CMAKE_CXX_STANDARD 14)

# For MinGW/GCC
if(MINGW OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_definitions(-DUSE_MATH_DEFINES)
endif()

include_directories(INCLUDE)

add_executable(C__
        EKF_Test.cpp
        EKF_GEOBS.cpp
        SRC/Matrix.cpp
        SRC/Mjday.cpp
        SRC/R_x.cpp
        SRC/R_y.cpp
        SRC/R_z.cpp
        SRC/TimeUpdate.cpp
        SRC/Legendre.cpp
        SRC/NutAngles.cpp
        SRC/SAT_Const.cpp
        SRC/Position.cpp
        SRC/Unit.cpp
        SRC/AzElPa.cpp
        SRC/Frac.cpp
        SRC/Sign_.cpp
)