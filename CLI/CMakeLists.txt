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

#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${MODULE_NAME}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}  
  ADDITIONAL_SRCS ${MODULE_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}  
  )

