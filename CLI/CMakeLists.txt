
#-----------------------------------------------------------------------------
set(MODULE_NAME PkModeling)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  ${${MODULE_NAME}_SOURCE_DIR}/PkSolver
  )

if(ITK_INSTALL_PREFIX)
  message("ITK_INSTALL_PREFIX is " ${ITK_INSTALL_PREFIX})
  list(APPEND MODULE_INCLUDE_DIRECTORIES
    ${ITK_INSTALL_PREFIX}/include/vxl/core
    ${ITK_INSTALL_PREFIX}/include/vxl/vcl
    )
endif()

#-----------------------------------------------------------------------------
set(MODULE_SRCS
  ${MODULE_NAME}.cxx
  itkSignalIntensityToS0ImageFilter.h
  itkSignalIntensityToS0ImageFilter.hxx
  itkSignalIntensityToConcentrationImageFilter.h
  itkSignalIntensityToConcentrationImageFilter.hxx
  itkConcentrationToQuantitativeImageFilter.h
  itkConcentrationToQuantitativeImageFilter.hxx
  )

#-----------------------------------------------------------------------------
set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES} PkSolver
  )

#
# ITK
#
set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # See Libs/ITKFactoryRegistration/CMakeLists.txt
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  RUNTIME_OUTPUT_DIRECTORY ${PkModeling_CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  LIBRARY_OUTPUT_DIRECTORY ${PkModeling_CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  ARCHIVE_OUTPUT_DIRECTORY ${PkModeling_CMAKE_ARCHIVE_OUTPUT_DIRECTORY}  
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()
