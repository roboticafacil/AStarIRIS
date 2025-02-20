# CMakeList.txt: archivo del proyecto de CMake de nivel superior, establezca la configuración global
# e incluya los subproyectos aquí.
#
cmake_minimum_required (VERSION 3.8)

# Project global configuration
project(AStarIRIS)

# Use C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Additional options for GCC/Clang
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
    add_compile_options(-Wall -Wextra -Wno-error=unused-parameter)
    add_compile_options(-Wall -Wextra -Wno-error)
    add_compile_options(-Wall -Wextra -Werror -Wno-sign-compare -Wno-unused-parameter)
endif()

# Windows setup
if (WIN32)
    message(STATUS "Setting up for Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
    add_compile_definitions(_ITERATOR_DEBUG_LEVEL=0)
endif()

# Linux/Unix setup
if (UNIX AND NOT WIN32)
    message(STATUS "Setting up for Linux/Unix")
    add_compile_options(-fPIC) # Common options for shared libraries in Linux
endif()

# Incluye los subproyectos
add_subdirectory(AStarIRIS)
