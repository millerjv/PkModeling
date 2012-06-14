# Configure XML description
set(PKM PkModeling)
configure_file(
  ${PKM}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PKM}.xml
  @ONLY
  )
  
# Additional includes
set(cli_module_include_directories
   ${PkModeling_SOURCE_DIR}/PkSolver
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
  )
  
# Additional Target libraries
set(cli_module_target_libraries
 ${ITK_LIBRARIES} PkSolver
)


#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${PKM}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}  
  ADDITIONAL_SRCS ${cli_module_SRCS}
  TARGET_LIBRARIES ${cli_module_target_libraries}  
  )

