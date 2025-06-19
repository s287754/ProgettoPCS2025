cmake_minimum_required(VERSION 3.9)

project(GeodesicPolyhedra LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND GeodesicPolyhedra_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND GeodesicPolyhedra_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
find_package(GTest REQUIRED)

list(APPEND GeodesicPolyhedra_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND GeodesicPolyhedra_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND GeodesicPolyhedra_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND GeodesicPolyhedra_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

# Insert Sources
################################################################################
add_subdirectory(src)

list(APPEND GeodesicPolyhedra_SOURCES ${GeodesicPolyhedra_sources})
list(APPEND GeodesicPolyhedra_HEADERS ${GeodesicPolyhedra_headers})
list(APPEND GeodesicPolyhedra_INCLUDE ${GeodesicPolyhedra_includes})

add_subdirectory(src_test)

list(APPEND GeodesicPolyhedra_TEST_SOURCES ${src_test_sources})
list(APPEND GeodesicPolyhedra_TEST_HEADERS ${src_test_headers})
list(APPEND GeodesicPolyhedra_TEST_INCLUDE ${src_test_includes})

add_subdirectory(ExportParaview)

list(APPEND GeodesicPolyhedra_SOURCES ${src_paraview_sources})
list(APPEND GeodesicPolyhedra_HEADERS ${src_paraview_headers})
list(APPEND GeodesicPolyhedra_INCLUDE ${src_paraview_includes})

# Create executable
################################################################################
add_executable(${PROJECT_NAME}
      main.cpp
      ${GeodesicPolyhedra_SOURCES}
      ${GeodesicPolyhedra_HEADERS})


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Polygon/Cell0Ds.csv ${CMAKE_CURRENT_SOURCE_DIR}/Polygon/Cell1Ds.csv ${CMAKE_CURRENT_SOURCE_DIR}/Polygon/Cell2Ds.csv
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(${PROJECT_NAME} ${GeodesicPolyhedra_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GeodesicPolyhedra_INCLUDE})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

add_executable(${PROJECT_NAME}_TEST
    main_test.cpp
    ${GeodesicPolyhedra_SOURCES}
    ${GeodesicPolyhedra_HEADERS}
    ${GeodesicPolyhedra_TEST_SOURCES}
    ${GeodesicPolyhedra_TEST_HEADERS})

target_link_libraries(${PROJECT_NAME}_TEST ${GeodesicPolyhedra_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME}_TEST PRIVATE ${GeodesicPolyhedra_INCLUDE} ${GeodesicPolyhedra_TEST_INCLUDE})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)