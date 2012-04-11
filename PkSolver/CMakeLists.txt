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
#SET(${LIBRARY_NAME}_SOURCE D:/Codes/Slicer4/Extensions/PkModeling/PkSolver/${LIBRARY_NAME}.cxx 
#			   ${CMAKE_CURRENT_SOURCE_DIR}/${LIBRARY_NAME}.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#ADD_LIBRARY(${LIBRARY_NAME}Lib SHARED ${LIBRARY_NAME}_SOURCE)
ADD_LIBRARY(${LIBRARY_NAME}Lib ${LIBRARY_NAME}.cxx ${LIBRARY_NAME}.h)
TARGET_LINK_LIBRARIES(${LIBRARY_NAME}Lib ${ITK_LIBRARIES})




