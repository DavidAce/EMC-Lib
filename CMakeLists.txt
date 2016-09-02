cmake_minimum_required(VERSION 2.8.9)
project(EMC-CLion)
 
#Bring the headers, such as Student.h into the project
include_directories(source)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -std=c++11 -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -std=c++11 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

file(GLOB SOURCES "source/*.cpp")
 
add_executable(EMC ${SOURCES} test.cpp EMC.h source/EMC.h source/objective_function.cpp source/objective_function.hpp source/constants.cpp)

target_compile_options(EMC PUBLIC -std=c++11 -fopenmp -Wall -Wfloat-conversion -O3 -march=native -mno-avx -Wno-deprecated-declarations -DEIGEN_NO_DEBUG -DNDEBUG)

