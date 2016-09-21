cmake_minimum_required(VERSION 2.8.9)
project(EMC-Lib)
 
include_directories(source)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -std=c++14 -Wno-deprecated-declarations")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -g -std=c++14 -march=native -Wno-deprecated-declarations -DNDEBUG")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
# Require MPI for this project:
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})



##################################
# Check for c++14 compiler support
##################################
include(CheckCXXCompilerFlag)
set(ENABLE_CXXFLAGS_TO_CHECK
        -std=c++14
        -std=c++1y
        -std=gnu++1y
        -std=gnu++14
        -std=c++11
        -std=gnu++11
        -std=gnu++0x
        -std=c++0x)

foreach(flag ${ENABLE_CXXFLAGS_TO_CHECK})
    string(REPLACE "-std=" "_" flag_var ${flag})
    string(REPLACE "+" "x" flag_var ${flag_var})
    check_cxx_compiler_flag("${flag}" COMPILER_HAS_CXX_FLAG${flag_var})
    if(COMPILER_HAS_CXX_FLAG${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        break()
    endif()
endforeach()


###################################
#Change flags depending on host
###################################
site_name(this_host)
if ("${this_host}" MATCHES "triolith")
    message("Host: ${this_host} (Intel) ")
    set(CMAKE_CXX_COMPILER mpiicpc)
    set(CMAKE_CXX_FLAGS_DEBUG "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g -O3 -ip -xavx -funroll-loops -DNDEBUG -cxxlib=/software/apps/gcc/5.3.0/build01/ -gxx-name=g++ -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -cxxlib=/software/apps/gcc/5.3.0/build01/ -gxx-name=g++ -ip -xavx -funroll-loops -DNDEBUG")
else()
    message("Host: ${this_host} (GNU)")
    set(CMAKE_CXX_COMPILER mpicxx)
    set(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wno-deprecated-declarations")
    set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -g -O3 -funroll-loops -march=native -frename-registers -Wno-deprecated-declarations -DNDEBUG")
endif()
message("Compiling with: ${CMAKE_CXX_COMPILER_ID}")


file(GLOB_RECURSE EMC_SOURCES "source/*.cpp")
file(GLOB_RECURSE EMC_HEADERS "source/*.hpp")

set (EMC_INCLUDE_DIRS "")
foreach (_headerFile ${EMC_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND EMC_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES EMC_INCLUDE_DIRS)

include_directories(${EMC_INCLUDE_DIRS})
add_executable(EMC example_main_field.cpp ${EMC_SOURCES} example_main_field.cpp)

#target_compile_options(EMC PUBLIC -std=c++14 -fopenmp -Wall -Wfloat-conversion -O3 -march=native -mno-avx -Wno-deprecated-declarations -DEIGEN_NO_DEBUG -DNDEBUG)

