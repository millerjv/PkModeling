# Configure XML description
set(PKM PkModeling)
configure_file(
  ${PKM}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PKM}.xml
  @ONLY
  )
  
# Additional includes
set(cli_module_include_directories
   ../PkSolver
  )

# Additional source files
set(cli_module_SRCS
	PkModeling.cxx
	itkS0CalculationFilter.h
	itkS0CalculationFilter.hxx
	itkConvertSignalIntensitiesToConcentrationValuesFilter.h
    itkConvertSignalIntensitiesToConcentrationValuesFilter.hxx
	itkCalculateQuantificationParametersFilter.h
	itkCalculateQuantificationParametersFilter.hxx
	../PkSolver/PkSolver.h
	../PkSolver/PkSolver.cxx
  )
  
# Additional Target libraries
# set(cli_module_target_libraries
#  ${ITK_LIBRARIES} PkModeling
#)


#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${PKM}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}  
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  ADDITIONAL_SRCS ${cli_module_SRCS}
  TARGET_LIBRARIES ${ITK_LIBRARIES}  
  )

#slicerMacroBuildCLI(
 # NAME FILL_THIS_IN_WITH_NAME_OF_CONCENTRATION_TO_PARAMETER_IMAGES
 # INCLUDE_DIRECTORIES ${cli_module_include_directories}
 # TARGET_LIBRARIES ${cli_module_target_libraries}
 # ADDITIONAL_SOURCES ${cli_module_SRCS}
 # )
