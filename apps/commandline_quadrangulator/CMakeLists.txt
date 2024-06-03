#cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 3.16)
project(quad_commandline)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -finline-functions -finline-limit=1000")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -finline-functions -finline-limit=1000")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-generate")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")


set(TRACEFIELD_DIR = "../../include" CACHE INTERNAL "")
set(EIGENPATH  = "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")

option(USE_DIRECTIONAL "Use directional feature" OFF)
if(USE_DIRECTIONAL)
set(DIRECTIONALPATH  = "../../external_lib/Directional/include/" CACHE INTERNAL "")
endif()
set(IGLLIB  = "../../external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${EIGENPATH})
option(USE_DIRECTIONAL "Use directional feature" OFF)
if(USE_DIRECTIONAL)
include_directories(${DIRECTIONALPATH})
endif()
include_directories(${IGLLIB})
add_executable( ${PROJECT_NAME} main.cpp)
