set(Launcher_Command ${Slicer_LAUNCH_COMMAND})

set(DCEMRIDATA ${CMAKE_CURRENT_SOURCE_DIR}/../Data/DCEMRIData)
set(DUKEDATA ${CMAKE_CURRENT_SOURCE_DIR}/../Data/DukeData)
set(SYNDATA ${CMAKE_CURRENT_SOURCE_DIR}/../Data/SyntheticData)

set(CLP ${EXTENSION_NAME})

add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

set(testname ${CLP}TestHelp)
add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint --help
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

set(testname ${CLP}TestPkModeling)
add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
	ModuleEntryPoint
	#--inputFourDNrrd  ${DCEMRIDATA}/case30Concentration.nrrd    
	#--inputFourDNrrd  ${SYNDATA}/SyntheticQulumeConcentration.nrrd    
	--inputFourDNrrd  ${DUKEDATA}/SyntheticDuke.nrrd    
	#--inputFourDNrrd  ${DUKEDATA}/SyntheticDukeSmallQulumeConcentration.nrrd    
	--T1Tissue 1000
    --T1Blood 1440
	--TR 5
	--FA 25
	--relaxivity 0.0049  
	--hematocrit 0.45
	--aucTimeInterval 90
	#--inputMask ${DCEMRIDATA}/maskVolume.nrrd
	#--outputKtrans ${DCEMRIDATA}/ktransVolume.nrrd  
	#--outputVe ${DCEMRIDATA}/veVolume.nrrd
	#--outputMaxSlope ${DCEMRIDATA}/maxSlopeVolume.nrrd
	#--outputAUC ${DCEMRIDATA}/AUCVolume.nrrd  
	#--inputMask ${SYNDATA}/SyntheticMask.nrrd
	#--outputKtrans ${SYNDATA}/ktransVolumeSyn.nrrd  
	#--outputVe ${SYNDATA}/veVolumeSyn.nrrd
	#--outputMaxSlope ${SYNDATA}/maxSlopeVolumeSyn.nrrd
	#--outputAUC ${SYNDATA}/AUCVolumeSyn.nrrd  
	--inputMask ${DUKEDATA}/SyntheticDukeMask.nrrd
	--outputKtrans ${DUKEDATA}/ktransVolumeSyn.nrrd  
	--outputVe ${DUKEDATA}/veVolumeSyn.nrrd
	--outputMaxSlope ${DUKEDATA}/maxSlopeVolumeSyn.nrrd
	--outputAUC ${DUKEDATA}/AUCVolumeSyn.nrrd  
	#--inputMask ${DUKEDATA}/SyntheticDukeSmallQulumeMask.nrrd
	#--outputKtrans ${DUKEDATA}/ktransVolumeSynSmallQulume.nrrd  
	#--outputVe ${DUKEDATA}/veVolumeSynSmallQulume.nrrd
	#--outputMaxSlope ${DUKEDATA}/maxSlopeVolumeSynSmallQulume.nrrd
	#--outputAUC ${DUKEDATA}/AUCVolumeSynSmallQulume.nrrd    
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
