# Configure XML description
set(SignalIntensitiesToConcentrationValuesCLI_NAME SignalIntensitiesToConcentrationValues)
configure_file(
  ${SignalIntensitiesToConcentrationValuesCLI_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${SignalIntensitiesToConcentrationValuesCLI_NAME}.xml
  @ONLY
  )
  
# Additional includes
set(cli_module_include_directories
   ${ITK_LIBRARIES} PkSolver
  )

# Additional source files
set(cli_module_SRCS
	SignalIntensitiesToConcentrationValues.cxx
	itkS0ForTimeSeriesInQulume.h
	itkS0ForTimeSeriesInQulume.hxx
	itkSignalIntensityToConcentration.h
        itkSignalIntensityToConcentration.hxx
  )
  
# Additional Target libraries
set(cli_module_target_libraries
  ${ITK_LIBRARIES} PkSolver
  )


#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME ${SignalIntensitiesToConcentrationValuesCLI_NAME}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}  
  ADDITIONAL_SRCS ${cli_module_SRCS}
  TARGET_LIBRARIES ${cli_module_target_libraries}  
  )

#slicerMacroBuildCLI(
 # NAME FILL_THIS_IN_WITH_NAME_OF_CONCENTRATION_TO_PARAMETER_IMAGES
 # INCLUDE_DIRECTORIES ${cli_module_include_directories}
 # TARGET_LIBRARIES ${cli_module_target_libraries}
 # ADDITIONAL_SOURCES ${cli_module_SRCS}
 # )
