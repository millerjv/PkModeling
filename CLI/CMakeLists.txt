# Additional includes
set(cli_module_include_directories
  )

# Additional source files
set(cli_module_SRCS
  )
  
# Additional Target libraries
set(cli_module_target_libraries
  ${ITK_LIBRARIES} PkSolver
  )


#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
slicerMacroBuildCLI(
  NAME FILL_THIS_IN_WITH_NAME_OF_INTENSITY_TO_CONCENTRATION
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  TARGET_LIBRARIES ${cli_module_target_libraries}
  ADDITIONAL_SOURCES ${cli_module_SRCS}
  )

slicerMacroBuildCLI(
  NAME FILL_THIS_IN_WITH_NAME_OF_CONCENTRATION_TO_PARAMETER_IMAGES
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  TARGET_LIBRARIES ${cli_module_target_libraries}
  ADDITIONAL_SOURCES ${cli_module_SRCS}
  )
