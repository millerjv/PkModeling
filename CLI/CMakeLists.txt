
#-----------------------------------------------------------------------------
set(MODULE_NAME PkModeling)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  ${${MODULE_NAME}_SOURCE_DIR}/PkSolver
  )

#-----------------------------------------------------------------------------
set(MODULE_SRCS
  ${MODULE_NAME}.cxx
  itkSignalIntensityToS0ImageFilter.h
  itkSignalIntensityToS0ImageFilter.hxx
  itkSignalIntensityToConcentrationImageFilter.h
  itkSignalIntensityToConcentrationImageFilter.h
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
#SEMMacroBuildCLI(
#  NAME ${MODULE_NAME}
#  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
#  ADDITIONAL_SRCS ${MODULE_SRCS}
#  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
#  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

