# CMakeList.txt: proyecto de CMake para AStarIRIS, incluya el origen y defina
# la lógica específica del proyecto aquí.
#

# Agregue un origen al ejecutable de este proyecto.
#add_executable (closestPoint_test "closestPoint_test.cpp" "closestPoint_test.h")
#set (CMAKE_CXX_STANDARD 11)
set(GCC_COMPILE_FLAGS "${GCC_COMPILE_FLAGS} /std:c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
#add_compile_definitions($<$<CONFIG:Debug>:_ITERATOR_DEBUG_LEVEL=0>)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include(FetchContent)
FetchContent_Declare(matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
        GIT_TAG origin/master)

        FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

LIST(APPEND CMAKE_PROGRAM_PATH  "C:/Program Files/gnuplot/bin")
add_compile_definitions(WITH_MATPLOT)

INCLUDE_DIRECTORIES("include")
INCLUDE_DIRECTORIES("C:/Program Files (x86)/Eigen3/include/eigen3")
INCLUDE_DIRECTORIES("C:/Program Files/Mosek/10.2/tools/platform/win64x86/h")
INCLUDE_DIRECTORIES("C:/Users/34696/source/qhull-2020.2/src")
#INCLUDE_DIRECTORIES("C:/Users/34696/source/repos/AStarIRIS/out/build/x64-releaseDeb/_deps/matplotplusplus-src/source")
#INCLUDE_DIRECTORIES("C:/Users/34696/source/matplotplusplus-1.2.1/source")

# TODO: Agregue pruebas y destinos de instalación si es necesario.
#add_subdirectory ("src")
#add_subdirectory ("tests")
SET(ConicLib_SRC "src/ConicSet.cpp" "src/EigenNdArray.cpp" "src/Range.cpp" "src/Point.cpp" "src/Circle.cpp" "src/Ellipsoid.cpp" "src/Polyhedron.cpp" "src/PolyhedronV.cpp" "src/PolyhedronObstacleCircularRobot.cpp")
SET(GraphLib_SRC "src/Graph.cpp" "src/Node.cpp" "src/PointNode.cpp" "src/CircleNode.cpp" "src/EllipsoidNode.cpp" "src/PolyhedronNode.cpp" "src/PolyhedronVNode.cpp" "src/PolyhedronObstacleCircularRobotNode.cpp" "src/PolyhedronTerminalNode.cpp")
SET(GCSLib_SRC "src/GCS.cpp" "src/NavGraph.cpp" "src/CObs.cpp" "src/CObsConic.cpp" "src/IRISConic.cpp"  "src/AStarIRISConic.cpp" "src/Solver.cpp" "src/MinDistanceSolver.cpp" "src/ConvexRelaxationMinDistanceSolver.cpp" "src/ConvexRelaxationSolver.cpp" "src/ExpandableIRISConic.cpp" "src/MCIRISConic.cpp")
SET(MosekLib "C:/Program Files/Mosek/10.2/tools/platform/win64x86/bin/fusion64_10_2.lib" "C:/Program Files/Mosek/10.2/tools/platform/win64x86/bin/mosek64_10_2.lib")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(QHullLib "C:/Users/34696/source/qhull-2020.2/lib/qhullstatic_d.lib" "C:/Users/34696/source/qhull-2020.2/lib/qhullcpp_d.lib")
endif()
     
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(QHullLib "C:/Users/34696/source/qhull-2020.2/lib/qhullstatic_r.lib" "C:/Users/34696/source/qhull-2020.2/lib/qhullcpp.lib")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    SET(QHullLib "C:/Users/34696/source/qhull-2020.2/lib/qhullstatic_r.lib" "C:/Users/34696/source/qhull-2020.2/lib/qhullcpp.lib")
endif()

ADD_LIBRARY(ConicLib STATIC ${ConicLib_SRC})
ADD_LIBRARY(GraphLib STATIC ${GraphLib_SRC} )
ADD_LIBRARY(GCSLib STATIC ${GCSLib_SRC})
add_executable (closestPoint_test "tests/closestPoint_test.cpp")
add_executable (removeConstraints_test "tests/removeConstraints_test.cpp")
add_executable (closestPointExpandingEllipsoid_test "tests/closestPointExpandingEllipsoid_test.cpp")
add_executable (conicSetsPerformance_test "tests/conicSetsPerformance.cpp")
add_executable (graph_test "tests/graph_test.cpp")
add_executable (IRIS_test "tests/IRIS_test.cpp")
add_executable (Expandable_IRIS_test "tests/Expandable_IRIS_test.cpp")
add_executable (AstarIRIS_phase1_test "tests/AStar_IRIS_test.cpp")
add_executable (MC_IRIS_test "tests/MC_IRIS_test.cpp")
add_executable (matplot_test "tests/matplot_test.cpp")

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET closestPoint_test PROPERTY CXX_STANDARD 20)
#endif()

TARGET_LINK_LIBRARIES(GraphLib ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(GCSLib ConicLib GraphLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(closestPoint_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(removeConstraints_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(closestPointExpandingEllipsoid_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(conicSetsPerformance_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(graph_test ConicLib GraphLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(Expandable_IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(AstarIRIS_phase1_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(MC_IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(matplot_test PUBLIC matplot ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})



