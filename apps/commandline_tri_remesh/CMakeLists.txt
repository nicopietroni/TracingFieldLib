#cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 3.16)
project(tri_remesh)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()
# if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
# endif()
# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
set(TRACEFIELD_DIR = "../../include" CACHE INTERNAL "")
set(EIGENPATH  = "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")
set(IGLLIB  = "../../external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${EIGENPATH})
include_directories(${IGLLIB})

add_executable( ${PROJECT_NAME} main.cpp)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if(ipo_supported)
    set_property(TARGET tri_remesh PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(WARNING "IPO not supported: ${ipo_error}")
endif()