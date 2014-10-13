set(Launcher_Command ${Slicer_LAUNCH_COMMAND})

set(TEMP ${PkModeling_BINARY_DIR}/Testing/Temporary)

#set(DCEMRIDATA ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/DCEMRIData)
#set(DUKEDATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/DukeData)
#set(SYNDATA ${CMAKE_CURRENT_SOURCE_DIR}/../../Data/SyntheticData)
set(QINPROSTATE001 ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/SampledPhantoms/QINProstate001)

set(CLP ${MODULE_NAME})

#-----------------------------------------------------------------------------
add_executable(${CLP}Test ${CLP}Test.cxx)
target_link_libraries(${CLP}Test ${CLP}Lib)
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testname ${CLP}TestHelp)
add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ModuleEntryPoint --help
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
if(QINPROSTATE001)
  set(testname QINProstate001)
  add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
    --compare ${QINPROSTATE001}/Baseline/phantom-ktrans.nrrd
    ${TEMP}/${testname}-ktrans.nrrd
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
                --outputKtrans ${TEMP}/${testname}-ktrans.nrrd 
                --outputVe ${TEMP}/${testname}-ve.nrrd
                --outputMaxSlope ${TEMP}/${testname}-maxslope.nrrd
                --outputAUC ${TEMP}/${testname}-auc.nrrd           
                --outputBAT ${TEMP}/${testname}-bat.nrrd           
                --fitted ${TEMP}/${testname}-fit.nrrd
                --concentrations ${TEMP}/${testname}-conc.nrrd
                --roiMask ${QINPROSTATE001}/Input/${testname}-phantom-ROI.nrrd
                --aifMask ${QINPROSTATE001}/Input/${testname}-phantom-AIF.nrrd
                ${QINPROSTATE001}/Input/${testname}-phantom.nrrd                   
                )
  set_property(TEST ${testname} PROPERTY LABELS ${CLP})
endif()

#-----------------------------------------------------------------------------
#if(DUKEDATA_DIR_DIR)
#  set(testname ${CLP}TestDukeSyntheticData)
#  add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
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
#  set_property(TEST ${testname} PROPERTY LABELS ${CLP})
#endif()


#-----------------------------------------------------------------------------
#if(DUKEDATA_DIR_DIR)
#  set(testname ${CLP}TestDukeSyntheticData)
#  add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
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
#  set_property(TEST ${testname} PROPERTY LABELS ${CLP})
#endif()

#-----------------------------------------------------------------------------
# set(testname ${CLP}TestSmallVolumeOfDukeSyntheticData)
# add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
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
# set_property(TEST ${testname} PROPERTY LABELS ${CLP})
 
# set(testname ${CLP}TestDukeSyntheticData)
# add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
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
# set_property(TEST ${testname} PROPERTY LABELS ${CLP})
 
# set(testname ${CLP}TestDCEMRIData)
# add_test(NAME ${testname} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
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
# set_property(TEST ${testname} PROPERTY LABELS ${CLP})
