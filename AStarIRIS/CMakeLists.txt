# Subdirectory project setup
file(GLOB_RECURSE SOURCES src/*.cpp)
file(GLOB_RECURSE HEADERS include/*.h)

# Eigen3 library is required
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(include)

# Windows setup
if (WIN32)
    message(STATUS "Setting up Windows dependencies")
    include_directories(
        "C:/Program Files (x86)/Eigen3/include/eigen3"
        "C:/Program Files/Mosek/10.2/tools/platform/win64x86/h"
        "$ENV{HOME}/source/qhull-2020.2/src"
    )

    # Setup Windows libraries
    set(MosekLib
        "C:/Program Files/Mosek/10.2/tools/platform/win64x86/bin/fusion64_10_2.lib"
        "C:/Program Files/Mosek/10.2/tools/platform/win64x86/bin/mosek64_10_2.lib"
    )

    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(QHullLib
            "$ENV{HOME}/source/qhull-2020.2/lib/qhullstatic_d.lib"
            "$ENV{HOME}/source/qhull-2020.2/lib/qhullcpp_d.lib"
        )
    elseif (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(QHullLib
            "$ENV{HOME}/source/qhull-2020.2/lib/qhullstatic_r.lib"
            "$ENV{HOME}/source/qhull-2020.2/lib/qhullcpp.lib"
        )
    endif()
endif()

# Unix setup
if (UNIX AND NOT WIN32)
    message(STATUS "Setting up Linux/Unix dependencies")

    include_directories("$ENV{HOME}/mosek_lib/mosek/10.2/tools/platform/linux64x86/h"
                        "$ENV{HOME}/project_ws/lib/eigen"
                        "/usr/include/libqhullcpp"
                        "/usr/include"
                        "$ENV{HOME}/AStarIRIS/AStarIRIS/include")

    
    set(QHullLib qhullcpp)

    set(MosekLib 
                    "$ENV{HOME}/mosek_lib/mosek/10.2/tools/platform/linux64x86/bin/libfusion64.so.10.2"
                    "$ENV{HOME}/mosek_lib/mosek/10.2/tools/platform/linux64x86/bin/libmosek64.so.10.2")
endif()

SET(ConicLib_SRC "src/ConicSet.cpp" "src/EigenUtils.cpp" "src/Range.cpp" "src/Point.cpp" "src/Sphere.cpp" "src/Ellipsoid.cpp" "src/Polyhedron.cpp" "src/PolyhedronV.cpp" "src/PolyhedronObstacleCircularRobot.cpp")
SET(GraphLib_SRC "src/Graph.cpp" "src/Node.cpp" "src/Edge.cpp" "src/PointNode.cpp" "src/SphereNode.cpp" "src/EllipsoidNode.cpp" "src/PolyhedronNode.cpp" "src/PolyhedronVNode.cpp" "src/PolyhedronObstacleCircularRobotNode.cpp" "src/PolyhedronTerminalNode.cpp")
SET(GCSLib_SRC   "src/GCS.cpp" "src/NavGraph.cpp" "src/CObs.cpp" "src/CObsConic.cpp" "src/IRISConic.cpp"  "src/AStarIRISConic.cpp" "src/SPP_GCS.cpp" "src/MIPMinDistanceSPP_GCS.cpp" "src/ConvexRelaxationMinDistanceSPP_GCS.cpp" "src/PerspectiveSPP_GCS.cpp" "src/ExpandableIRISConic.cpp" "src/MCIRISConic.cpp" "src/ConvexRelaxationBezierCurveSPP_GCS.cpp")

ADD_LIBRARY(ConicLib STATIC ${ConicLib_SRC})
ADD_LIBRARY(GraphLib STATIC ${GraphLib_SRC} )
ADD_LIBRARY(GCSLib STATIC ${GCSLib_SRC})

add_executable (closestPoint_test "tests/closestPoint_test.cpp")
add_executable (isInsideSeparatingHyperplanes_test "tests/isInsideSeparatingHyperplanes_test.cpp")
add_executable (removeConstraints_test "tests/removeConstraints_test.cpp")
add_executable (closestPointExpandingEllipsoid_test "tests/closestPointExpandingEllipsoid_test.cpp")
add_executable (conicSetsPerformance_test "tests/conicSetsPerformance.cpp")
add_executable (graph_test "tests/graph_test.cpp")
add_executable (IRIS_test "tests/IRIS_test.cpp")
add_executable (Expandable_IRIS_test "tests/Expandable_IRIS_test.cpp")
add_executable (graph_cycle_analysis "tests/graph_cycle_analysis.cpp")
add_executable (AstarIRIS_test "tests/AStar_IRIS_test.cpp")
add_executable (MC_IRIS_test "tests/MC_IRIS_test.cpp")
#add_executable (AstarIRIS_maze "tests/Astar_IRIS_maze.cpp")
add_executable (Multiagent_IRIS_test "tests/Multiagent_IRIS_test.cpp")

#if (CMAKE_VERSION VERSION_GREATER 3.12)
#  set_property(TARGET closestPoint_test PROPERTY CXX_STANDARD 20)
#endif()

TARGET_LINK_LIBRARIES(GraphLib ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(GCSLib ConicLib GraphLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(closestPoint_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(isInsideSeparatingHyperplanes_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(removeConstraints_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(closestPointExpandingEllipsoid_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(conicSetsPerformance_test ConicLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(graph_test ConicLib GraphLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(Expandable_IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(graph_cycle_analysis ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(AstarIRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(MC_IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
#TARGET_LINK_LIBRARIES(AstarIRIS_maze ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})
TARGET_LINK_LIBRARIES(Multiagent_IRIS_test ConicLib GraphLib GCSLib ${MosekLib} ${QHullLib} ${EIGEN3_LIBRARIES})




