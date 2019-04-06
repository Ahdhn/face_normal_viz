cmake_minimum_required(VERSION 3.1)
project(face_normal_viz_libigl)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

set( INPUT_PATH "input/" )
get_filename_component(INPUT_PATH ${INPUT_PATH} ABSOLUTE)
message( "INPUT_PATH= " ${INPUT_PATH})
add_definitions( -DINPUT_DIR=${INPUT_PATH})

# Add your project files
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw)
