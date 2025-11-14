cmake_minimum_required(VERSION 3.5)
project(icp)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(IMGUI_DIR "/Users/134750/Desktop/imgui/" CACHE INTERNAL "")
set(TRACEFIELD_DIR = "../../include" CACHE INTERNAL "")
set(GLDRAW_DIR = "../include/GLFunctions/" CACHE INTERNAL "")
set(EIGENPATH  = "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")

set(IGLLIB  = "../../external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${IMGUI_DIR})
include_directories(${IMGUI_DIR}/backends)
include_directories(/opt/homebrew/include)
include_directories(${EIGENPATH})


include_directories(${IGLLIB})

#include_directories(${BROWS_DIR})
# Find GLEW (installed via Homebrew)
find_library(GLEW_LIB NAMES GLEW glew PATHS /opt/homebrew/lib /usr/local/lib)
find_path(GLEW_INCLUDE_DIR GL/glew.h PATHS /opt/homebrew/include /usr/local/include)
if(NOT GLEW_LIB OR NOT GLEW_INCLUDE_DIR)
  message(FATAL_ERROR "GLEW not found")
endif()

add_executable(${PROJECT_NAME} main.cpp
               ${IMGUI_DIR}/imgui.cpp
               ${IMGUI_DIR}/imgui_demo.cpp
               ${IMGUI_DIR}/imgui_draw.cpp
               ${IMGUI_DIR}/imgui_tables.cpp
               ${IMGUI_DIR}/imgui_widgets.cpp
               ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
               ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
               ../../external_lib/FileBrowser/ImGuiFileBrowser.cpp
               #${GLEWPATH}/src/glew.c
               )

target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIR})


target_link_libraries(${PROJECT_NAME}
    ${GLEW_LIB} 
    "-framework OpenGL"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
    "-L/opt/homebrew/lib"
    "-lglfw"
)
