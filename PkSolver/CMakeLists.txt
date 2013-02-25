cmake_minimum_required(VERSION 2.8.7)

set(LIBRARY_NAME PkSolver)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(${LIBRARY_NAME} SHARED ${LIBRARY_NAME}.cxx ${LIBRARY_NAME}.h)
target_link_libraries(${LIBRARY_NAME} ${ITK_LIBRARIES})

