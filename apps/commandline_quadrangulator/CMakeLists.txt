#cmake_minimum_required(VERSION 3.1)
cmake_minimum_required(VERSION 3.16)
project(quad_commandline)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # macOS specific flags
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
else ()
    # Flags for other operating systems, such as Windows
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif ()
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(LUT_PATH "../../../include/Lut.txt")
add_definitions(-DLUT_PATH=\"${LUT_PATH}\")

# Enable AddressSanitizer and UndefinedBehaviorSanitizer in debug mode
# set(SANITIZE_FLAGS "-fsanitize=address -fsanitize=undefined")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SANITIZE_FLAGS} -g -O1")
# set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${SANITIZE_FLAGS}")


list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-generate")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-use")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -fstrict-aliasing -march=native")


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

#if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
#endif()