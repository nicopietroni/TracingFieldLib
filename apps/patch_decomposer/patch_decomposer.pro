TEMPLATE = app
CONFIG += app_bundle c++11
CONFIG -= qt

TARGET = patch_decomposer

TRACEFIELD_DIR = "../../include"
IMGUI_DIR  = "/Users/134750/Desktop/imgui/"
GLEWPATH   = "/Users/134750/Desktop/code/lib/glew"
IGLLIB    = "../../external_lib/libigl/include/"
#DIRECTIONALPATH = "../../external_lib/Directional/include/"
#BROWS_DIR  = "../../external_lib/FileBrowser"
EIGENPATH  = "../../external_lib/eigen-3.4.0/"

INCLUDEPATH += $$IMGUI_DIR
INCLUDEPATH += $$TRACEFIELD_DIR
INCLUDEPATH += $$IMGUI_DIR/backends
INCLUDEPATH += /opt/homebrew/include
INCLUDEPATH += $$GLEWPATH/include
INCLUDEPATH += $$EIGENPATH
INCLUDEPATH += $$IGLLIB
INCLUDEPATH += $$DIRECTIONALPATH
INCLUDEPATH += $$BROWS_DIR

SOURCES = main.cpp
SOURCES += $$IMGUI_DIR/imgui.cpp $$IMGUI_DIR/imgui_demo.cpp $$IMGUI_DIR/imgui_draw.cpp $$IMGUI_DIR/imgui_tables.cpp $$IMGUI_DIR/imgui_widgets.cpp
SOURCES += $$IMGUI_DIR/backends/imgui_impl_glfw.cpp $$IMGUI_DIR/backends/imgui_impl_opengl2.cpp
SOURCES += $$GLEWPATH/src/glew.c
#SOURCES += $$BROWS_DIR/ImGuiFileBrowser.cpp
HEADERS += ../io_mesh/io_obj.h
HEADERS += $$TRACEFIELD_DIR/UI/orbitball.h
HEADERS += $$GLDRAW_DIR/GLDrawBase.h
HEADERS += $$TRACEFIELD_DIR/math_utils.h


HEADERS += CMakeLists.txt

LIBS += -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo
LIBS += -L/opt/homebrew/lib -lglfw
