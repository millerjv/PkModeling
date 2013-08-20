cmake_minimum_required(VERSION 2.8.7)

set(LIBRARY_NAME PkSolver)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(${LIBRARY_NAME} STATIC ${LIBRARY_NAME}.cxx ${LIBRARY_NAME}.h)
target_link_libraries(${LIBRARY_NAME} ${ITK_LIBRARIES})
if (CMAKE_SYSTEM MATCHES "Linux")
  set_target_properties(${LIBRARY_NAME} PROPERTIES COMPILE_FLAGS "-fPIC")
endif ()

# No need to install the library as it is statically linked to other libraries and executables in this extension
#install(TARGETS ${LIBRARY_NAME}
#     RUNTIME DESTINATION ${PkModeling_INSTALL_SUPPORT_RUNTIME_DESTINATION} COMPONENT RuntimeLibraries
#     LIBRARY DESTINATION ${PkModeling_INSTALL_SUPPORT_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
#     ARCHIVE DESTINATION ${PkModeling_INSTALL_SUPPORT_ARCHIVE_DESTINATION} COMPONENT Development
#     )

