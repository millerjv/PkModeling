
set(Launcher_Command ${Slicer_LAUNCH_COMMAND})

set(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../../Data)

set(CLP ${EXTENSION_NAME})

add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

set(testname ${CLP}Test000)
add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint --help
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

set(testname ${CLP}Test)
add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
	ModuleEntryPoint
	--inputFourDNrrd  ${INPUT}/FourDImages.nrrd    
	--T1 1600
	--TR 3.984
	--FA 15
	--relaxivity 0.0049
	--outputFourDNrrd  ${INPUT}/ConcentrationQulume.nrrd   
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
