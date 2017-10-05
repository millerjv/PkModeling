set(Launcher_Command ${Slicer_LAUNCH_COMMAND})

set(TEMP ${PkModeling_BINARY_DIR}/Testing/Temporary)

#set(DCEMRIDATA ${CMAKE_SOURCE_DIR}/Data/DCEMRIData)
#set(DUKEDATA_DIR ${CMAKE_SOURCE_DIR}/Data/DukeData)
#set(SYNDATA ${CMAKE_SOURCE_DIR}/Data/SyntheticData)
set(QINPROSTATE001 ${CMAKE_SOURCE_DIR}/Data/SampledPhantoms/QINProstate001)

set(CLP ${MODULE_NAME})

#-----------------------------------------------------------------------------
add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib ${SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES})
set_target_properties(${CLP}Test PROPERTIES COMPILE_FLAGS "-DMODULE_STATIC")
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})
if(TARGET ITKFactoryRegistration)
  target_compile_definitions(${CLP}Test PUBLIC HAS_ITK_FACTORY_REGISTRATION)
endif()

#-----------------------------------------------------------------------------
set(testName ${CLP}TestHelp)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint --help
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
if(QINPROSTATE001)
  set(testName QINProstate001)
  add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
    --compare ${QINPROSTATE001}/Baseline/phantom-ktrans.nrrd
    ${TEMP}/${testName}-ktrans.nrrd
    ModuleEntryPoint        
                --T1Tissue 1597
                --T1Blood 1600
                --relaxivity 0.0039
                --S0grad 15.0
                --hematocrit 0.4
                --aucTimeInterval 90
                --fTolerance 1e-4 
                --gTolerance 1e-4 
                --xTolerance 1e-5 
                --epsilon 1e-9 
                --maxIter 200 
                --outputKtrans ${TEMP}/${testName}-ktrans.nrrd 
                --outputVe ${TEMP}/${testName}-ve.nrrd
                --outputMaxSlope ${TEMP}/${testName}-maxslope.nrrd
                --outputAUC ${TEMP}/${testName}-auc.nrrd           
                --outputBAT ${TEMP}/${testName}-bat.nrrd           
                --fitted ${TEMP}/${testName}-fit.nrrd
                --concentrations ${TEMP}/${testName}-conc.nrrd
                --roiMask ${QINPROSTATE001}/Input/${testName}-phantom-ROI.nrrd
                --aifMask ${QINPROSTATE001}/Input/${testName}-phantom-AIF.nrrd
                ${QINPROSTATE001}/Input/${testName}-phantom.nrrd                   
                )
  set_property(TEST ${testName} PROPERTY LABELS ${CLP})
endif()


#-----------------------------------------------------------------------------
# Regression tests
#-----------------------------------------------------------------------------
include(RegressionTests.cmake)



#-----------------------------------------------------------------------------
#if(DUKEDATA_DIR_DIR)
#  set(testName ${CLP}TestDukeSyntheticData)
#  add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
#    ModuleEntryPoint        
#                --T1Tissue 1000
#                --T1Blood 1440
#                --relaxivity 0.0049 
#                --hematocrit 0.45
#                --aucTimeInterval 90
#                ${DUKEDATA_DIR}/SyntheticDuke.nrrd                   
#                ${DUKEDATA_DIR}/SyntheticDukeMask.nrrd
#                --outputKtrans ${TEMP}/ktransVolumeSyn.nrrd 
#                --outputVe ${TEMP}/veVolumeSyn.nrrd
#                --outputMaxSlope ${TEMP}/maxSlopeVolumeSyn.nrrd
#                --outputAUC ${TEMP}/AUCVolumeSyn.nrrd           
#                )
#  set_property(TEST ${testName} PROPERTY LABELS ${CLP})
#endif()


#-----------------------------------------------------------------------------
#if(DUKEDATA_DIR_DIR)
#  set(testName ${CLP}TestDukeSyntheticData)
#  add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
#    ModuleEntryPoint        
#                --T1Tissue 1000
#                --T1Blood 1440
#                --relaxivity 0.0049 
#                --hematocrit 0.45
#                --aucTimeInterval 90
#                ${DUKEDATA_DIR}/SyntheticDuke.nrrd                   
#                ${DUKEDATA_DIR}/SyntheticDukeMask.nrrd
#                --outputKtrans ${TEMP}/ktransVolumeSyn.nrrd 
#                --outputVe ${TEMP}/veVolumeSyn.nrrd
#                --outputMaxSlope ${TEMP}/maxSlopeVolumeSyn.nrrd
#                --outputAUC ${TEMP}/AUCVolumeSyn.nrrd           
#                )
#  set_property(TEST ${testName} PROPERTY LABELS ${CLP})
#endif()

#-----------------------------------------------------------------------------
# set(testName ${CLP}TestSmallVolumeOfDukeSyntheticData)
# add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
#                 ModuleEntryPoint
#                 --T1Tissue 1000
#                 --T1Blood 1440
#                 #--TR 5
#                 #--FA 25
#                 --relaxivity 0.0049 
#                 --hematocrit 0.45
#                 --aucTimeInterval 90
#                 ${DUKEDATA_DIR}/SyntheticDukeSmallVolume.nrrd   
#                 ${DUKEDATA_DIR}/SyntheticDukeSmallVolumeMask.nrrd
#                 --outputKtrans ${TEMP}/ktransVolumeSynSmallVolume.nrrd 
#                 --outputVe ${TEMP}/veVolumeSynSmallVolume.nrrd
#                 --outputMaxSlope ${TEMP}/maxSlopeVolumeSynSmallVolume.nrrd
#                 --outputAUC ${TEMP}/AUCVolumeSynSmallVolume.nrrd   
#   )
# set_property(TEST ${testName} PROPERTY LABELS ${CLP})
 
# set(testName ${CLP}TestDukeSyntheticData)
# add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
#   ModuleEntryPoint        
#                 --T1Tissue 1000
#                 --T1Blood 1440
#                 #--TR 5
#                 #--FA 25
#                 --relaxivity 0.0049 
#                 --hematocrit 0.45
#                 --aucTimeInterval 90
#                 ${DUKEDATA_DIR}/SyntheticDuke.nrrd                   
#                 ${DUKEDATA_DIR}/SyntheticDukeMask.nrrd
#                 --outputKtrans ${TEMP}/ktransVolumeSyn.nrrd 
#                 --outputVe ${TEMP}/veVolumeSyn.nrrd
#                 --outputMaxSlope ${TEMP}/maxSlopeVolumeSyn.nrrd
#                 --outputAUC ${TEMP}/AUCVolumeSyn.nrrd           
#   )
# set_property(TEST ${testName} PROPERTY LABELS ${CLP})
 
# set(testName ${CLP}TestDCEMRIData)
# add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
#                 ModuleEntryPoint
#                 --T1Tissue 1597
#                 --T1Blood 1600
#                 #--TR 3.984
#                 #--FA 15
#                 --relaxivity 0.0049 
#                 --hematocrit 0.4
#                 --aucTimeInterval 90
#                 ${DCEMRIDATA}/OriginalVolume.nrrd             
#                 ${DCEMRIDATA}/maskVolume.nrrd
#                 --outputKtrans ${TEMP}/ktransVolume.nrrd 
#                 --outputVe ${TEMP}/veVolume.nrrd
#                 --outputMaxSlope ${TEMP}/maxSlopeVolume.nrrd
#                 --outputAUC ${TEMP}/AUCVolume.nrrd             
#   )
# set_property(TEST ${testName} PROPERTY LABELS ${CLP})
