# CMakeList.txt: archivo del proyecto de CMake de nivel superior, establezca la configuración global
# e incluya los subproyectos aquí.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()
#add_compile_definitions($<$<CONFIG:Debug>:_ITERATOR_DEBUG_LEVEL=0>)
#SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
#SET( CMAKE_MODULE_PATH "$ENV{EIGEN3_DIR}" )

#set(CMAKE_PREFIX_PATH "C:/Users/34696/source/eigen-3.4.0/cmake")
#set(Eigen3_DIR "C:/Users/34696/source/eigen-3.4.0/cmake")
#MESSAGE(${CMAKE_PREFIX_PATH})
#MESSAGE(${Eigen3_DIR})

project ("AStarIRIS")

# Incluya los subproyectos.
add_subdirectory ("AStarIRIS")
