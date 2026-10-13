cmake_minimum_required(VERSION 3.5)
project(quad_commandline LANGUAGES CXX)

cmake_policy(SET CMP0069 NEW)

# Use C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release build if not specified
#if(NOT CMAKE_BUILD_TYPE)
# set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
#endif()

# Options: Release, Debug, Sanity, RelDebug
set(MY_BUILD_TYPE "Release" CACHE STRING "Choose build type: Release, Debug, Sanity, RelDebug")
set_property(CACHE MY_BUILD_TYPE PROPERTY STRINGS Release Debug Sanity RelDebug)

if(MY_BUILD_TYPE STREQUAL "Sanity")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
    add_compile_definitions(SANITY_CHECKS)

    add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address,undefined)

elseif(MY_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)

elseif(MY_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

elseif(MY_BUILD_TYPE STREQUAL "RelDebug")
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
else()
    message(FATAL_ERROR "Unknown MY_BUILD_TYPE: ${MY_BUILD_TYPE}")
endif()

message(STATUS "MY_BUILD_TYPE = ${MY_BUILD_TYPE}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Define LUT path as a compile-time macro
set(LUT_PATH "../../../include/Lut.txt")
add_compile_definitions(LUT_PATH="${LUT_PATH}")

# Include custom CMake modules
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Include directories
set(TRACEFIELD_DIR "../../include" CACHE INTERNAL "")
set(EIGENPATH "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")
set(IGLLIB "../../external_lib/libigl/include/" CACHE INTERNAL "")

# Create executable
add_executable(${PROJECT_NAME} main.cpp)

# Add include directories to target
target_include_directories(${PROJECT_NAME} PRIVATE
  ${TRACEFIELD_DIR}
  ${EIGENPATH}
  ${IGLLIB}
)

# Enable threading support
find_package(Threads REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Threads::Threads)

# Enable Interprocedural Optimization (LTO) if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported AND NOT MY_BUILD_TYPE STREQUAL "RelDebug")
  set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
  message(WARNING "IPO/LTO not supported: ${ipo_error}")
endif()

# Compiler-specific optimization flags
if (MSVC)

  # Enable full optimization and link-time optimization in Release mode
  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Release>:/O2 /DNDEBUG /GL>
    $<$<CONFIG:RelWithDebInfo>:/O2 /DNDEBUG>
  )

  # Enable link-time code generation
  target_link_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Release>:/LTCG>
  )

  # Allow larger object files
  target_compile_options(${PROJECT_NAME} PRIVATE /bigobj)

else()

  # Enable high-level optimization in Release mode
  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    $<$<CONFIG:RelWithDebInfo>:-O3 -DNDEBUG -g>
  )

  # Optimize specifically for the local CPU architecture
  # Improves vectorization and instruction scheduling
  # (also in RelWithDebInfo, so the profiled code matches Release)
  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>:-march=native>
  )

  # Optional: enable loop unrolling (can help in heavy numeric code)
  target_compile_options(${PROJECT_NAME} PRIVATE $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>:-funroll-loops>)

  # Optional (use carefully): faster floating point math
  # May break strict IEEE compliance
  # target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:Release>:-ffast-math>)

endif()