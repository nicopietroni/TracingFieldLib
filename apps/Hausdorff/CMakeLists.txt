cmake_minimum_required(VERSION 3.16)
project(hausdorff)

set(CMAKE_CXX_STANDARD 17)

# ==============================
# ===== USER SWITCHES (EDIT ONLY THESE) =====
# ==============================
set(BUILD_MODE "Release")      # Debug | Release
set(ENABLE_SANITIZER OFF)      # ON | OFF

# ==============================
set(CMAKE_BUILD_TYPE ${BUILD_MODE} CACHE STRING "Build type" FORCE)
message(STATUS "Build mode: ${BUILD_MODE}")
message(STATUS "Sanitizer:  ${ENABLE_SANITIZER}")

set(TRACEFIELD_DIR "../../include" CACHE INTERNAL "")
set(EIGENPATH "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")
set(IGLLIB "../../external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${EIGENPATH})
include_directories(${IGLLIB})

add_executable(${PROJECT_NAME} main.cpp)

# Build flags
if(BUILD_MODE STREQUAL "Debug")
    target_compile_options(${PROJECT_NAME} PRIVATE -g -O0)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -O3 -DNDEBUG)
endif()

# Sanitizer flags
if(ENABLE_SANITIZER)
    message(STATUS "ASan + UBSan enabled")
    target_compile_options(${PROJECT_NAME} PRIVATE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
    )
    target_link_options(${PROJECT_NAME} PRIVATE
        -fsanitize=address,undefined
    )
endif()

# IPO solo in Release
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported AND BUILD_MODE STREQUAL "Release")
    set_property(TARGET hausdorff PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()