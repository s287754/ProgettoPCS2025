list(APPEND GeodesicPolyhedra_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND GeodesicPolyhedra_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)

list(APPEND GeodesicPolyhedra_headers ${CMAKE_CURRENT_SOURCE_DIR}/Geometry.hpp)
list(APPEND GeodesicPolyhedra_sources ${CMAKE_CURRENT_SOURCE_DIR}/Geometry.cpp)

list(APPEND GeodesicPolyhedra_headers ${CMAKE_CURRENT_SOURCE_DIR}/DataGenerator.hpp)
list(APPEND GeodesicPolyhedra_sources ${CMAKE_CURRENT_SOURCE_DIR}/DataGenerator.cpp)
list(APPEND GeodesicPolyhedra_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(GeodesicPolyhedra_sources ${GeodesicPolyhedra_sources} PARENT_SCOPE)
set(GeodesicPolyhedra_headers ${GeodesicPolyhedra_headers} PARENT_SCOPE)
set(GeodesicPolyhedra_includes ${GeodesicPolyhedra_includes} PARENT_SCOPE)