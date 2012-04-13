#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.4)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
#-----------------------------------------------------------------------------
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

# project (PkSolver) 
SET(LIBRARY_NAME PkSolver)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_LIBRARY(${LIBRARY_NAME} SHARED ${LIBRARY_NAME}.cxx ${LIBRARY_NAME}.h)
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${ITK_LIBRARIES})




