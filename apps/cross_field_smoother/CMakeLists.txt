cmake_minimum_required(VERSION 3.1)
project(cross_field_smoother)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(IMGUI_DIR "/Users/134750/Desktop/imgui/" CACHE INTERNAL "")
set(GLEWPATH "/Users/134750/Desktop/code/lib/glew" CACHE INTERNAL "")
set(TRACEFIELD_DIR = "../../include" CACHE INTERNAL "")
set(GLDRAW_DIR = "../include/GLFunctions/" CACHE INTERNAL "")
set(EIGENPATH  = "../../external_lib/eigen-3.4.0/" CACHE INTERNAL "")
set(DIRECTIONALPATH  = "../../external_lib/Directional/include/" CACHE INTERNAL "")
set(IGLLIB  = "../../external_lib/libigl/include/" CACHE INTERNAL "")
#set(BROWS_DIR  = "../../external_lib/FileBrowser" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${IMGUI_DIR})
include_directories(${IMGUI_DIR}/backends)
include_directories(/opt/homebrew/include)
include_directories(${GLEWPATH}/include)
include_directories(${EIGENPATH})
include_directories(${DIRECTIONALPATH})
include_directories(${IGLLIB})
#include_directories(${BROWS_DIR})


add_executable(${PROJECT_NAME} main.cpp
               ${IMGUI_DIR}/imgui.cpp
               ${IMGUI_DIR}/imgui_demo.cpp
               ${IMGUI_DIR}/imgui_draw.cpp
               ${IMGUI_DIR}/imgui_tables.cpp
               ${IMGUI_DIR}/imgui_widgets.cpp
               ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
               ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
               #../../external_lib/FileBrowser/ImGuiFileBrowser.cpp
               ${GLEWPATH}/src/glew.c)


target_link_libraries(${PROJECT_NAME}
    "-framework OpenGL"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
    "-L/opt/homebrew/lib"
    "-lglfw"
)
