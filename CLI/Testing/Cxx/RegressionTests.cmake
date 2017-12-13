#-----------------------------------------------------------------------------
# Convenience functions to get command snippets needed for most tests
#-----------------------------------------------------------------------------
function(set_compareArgs computeFpv)
  set(compareArgs --compareIntensityTolerance 1e-4
                  --compare ${referenceDataBaseDir}-conc.nrrd
                  ${tempOutDataBaseName}-conc.nrrd
                  --compare ${referenceDataBaseName}-ktrans.nrrd
                  ${tempOutDataBaseName}-ktrans.nrrd
                  --compare ${referenceDataBaseName}-ve.nrrd
                  ${tempOutDataBaseName}-ve.nrrd
                  --compare ${referenceDataBaseName}-maxslope.nrrd
                  ${tempOutDataBaseName}-maxslope.nrrd
                  --compare ${referenceDataBaseName}-auc.nrrd
                  ${tempOutDataBaseName}-auc.nrrd
                  --compare ${referenceDataBaseName}-rsq.nrrd
                  ${tempOutDataBaseName}-rsq.nrrd
                  --compare ${referenceDataBaseName}-bat.nrrd
                  ${tempOutDataBaseName}-bat.nrrd
                  --compare ${referenceDataBaseName}-fit.nrrd
                  ${tempOutDataBaseName}-fit.nrrd
                  --compare ${referenceDataBaseName}-diag.nrrd
                  ${tempOutDataBaseName}-diag.nrrd)
  if(${computeFpv})
    set(compareArgs ${compareArgs}
                    --compare ${referenceDataBaseName}-fpv.nrrd
                    ${tempOutDataBaseName}-fpv.nrrd)
  endif()
  # make result available in parent scope
  set(compareArgs ${compareArgs} PARENT_SCOPE)
endfunction()


#-----------------------------------------------------------------------------
function(set_paramsArgs computeFpv)
  set(paramsArgs --T1Tissue 1597
                 --T1Blood 1600
                 --relaxivity 0.0039
                 --S0grad 15.0
                 --hematocrit 0.4
                 --aucTimeInterval 90
                 --fTolerance 1e-4 
                 --gTolerance 1e-4 
                 --xTolerance 1e-5 
                 --epsilon 1e-9 
                 --maxIter 200)
  if(${computeFpv})
    set(paramsArgs ${paramsArgs}
                   --computeFpv)
  endif()
  # make result available in parent scope
  set(paramsArgs ${paramsArgs} PARENT_SCOPE)
endfunction()


#-----------------------------------------------------------------------------
function(set_outputParamsArgs computeFpv)
  set(outputParamsArgs --concentrations ${tempOutDataBaseName}-conc.nrrd
                       --outputKtrans ${tempOutDataBaseName}-ktrans.nrrd 
                       --outputVe ${tempOutDataBaseName}-ve.nrrd
                       --outputMaxSlope ${tempOutDataBaseName}-maxslope.nrrd
                       --outputAUC ${tempOutDataBaseName}-auc.nrrd      
                       --outputRSquared ${tempOutDataBaseName}-rsq.nrrd
                       --outputBAT ${tempOutDataBaseName}-bat.nrrd           
                       --fitted ${tempOutDataBaseName}-fit.nrrd
                       --outputDiagnostics ${tempOutDataBaseName}-diag.nrrd)
  if(${computeFpv})
    set(outputParamsArgs ${outputParamsArgs}
                         --outputFpv ${tempOutDataBaseName}-fpv.nrrd)
  endif()
  # make result available in parent scope
  set(outputParamsArgs ${outputParamsArgs} PARENT_SCOPE)  
endfunction()



#-----------------------------------------------------------------------------
# Regression Tests QINProstate001
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_SOURCE_DIR}/Data/RegressionTests/QINProstate001/Input/QINProstate001-phantom)
set(referenceDataBaseDir ${CMAKE_SOURCE_DIR}/Data/RegressionTests/QINProstate001/Reference/)

#-----------------------------------------------------------------------------
set(testName QINProstate001_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set_paramsArgs(FALSE)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName QINProstate001_AllOutputsInclFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(TRUE)
set_paramsArgs(TRUE)
set_outputParamsArgs(TRUE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})


#-----------------------------------------------------------------------------
# Regression Tests QINBreast001
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/QINBreast001/Input/QINBreast001-phantom)
set(referenceDataBaseDir ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/QINBreast001/Reference/)

#-----------------------------------------------------------------------------
set(testName QINBreast001_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set_paramsArgs(FALSE)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    --usePopAif
    ${outputParamsArgs}
    ${inputDataBaseName}.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName QINBreast001_ConstantBat)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set_paramsArgs(FALSE)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    --usePopAif
    --BATCalculationMode UseConstantBAT
    --constantBAT 4
    ${outputParamsArgs}
    ${inputDataBaseName}.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})


#-----------------------------------------------------------------------------
# Regression Tests DROs
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/DROs/Input/DRO)
set(referenceDataBaseDir ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/DROs/Reference/)

#-----------------------------------------------------------------------------
set(testName DRO5min1secinf_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set(paramsArgs --T1Tissue 1434
               --T1Blood 1600
               --relaxivity 0.0037
               --S0grad 15.0
               --hematocrit 0.45
               --aucTimeInterval 90
               --fTolerance 1e-4 
               --gTolerance 1e-4 
               --xTolerance 1e-5 
               --epsilon 1e-9 
               --maxIter 200)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}5min1secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName DRO5min3secinf_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set(paramsArgs --T1Tissue 1434
               --T1Blood 1600
               --relaxivity 0.0037
               --S0grad 15.0
               --hematocrit 0.45
               --aucTimeInterval 90
               --fTolerance 1e-4 
               --gTolerance 1e-4 
               --xTolerance 1e-5 
               --epsilon 1e-9 
               --maxIter 200)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}5min3secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName DRO5min5secinf_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set(paramsArgs --T1Tissue 1434
               --T1Blood 1600
               --relaxivity 0.0037
               --S0grad 15.0
               --hematocrit 0.45
               --aucTimeInterval 90
               --fTolerance 1e-4 
               --gTolerance 1e-4 
               --xTolerance 1e-5 
               --epsilon 1e-9 
               --maxIter 200)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}5min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName DRO3min3secinf_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set(paramsArgs --T1Tissue 1434
               --T1Blood 1600
               --relaxivity 0.0037
               --S0grad 15.0
               --hematocrit 0.45
               --aucTimeInterval 90
               --fTolerance 1e-4 
               --gTolerance 1e-4 
               --xTolerance 1e-5 
               --epsilon 1e-9 
               --maxIter 200)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}3min3secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName DRO3min5secinf_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}${testName})
set_compareArgs(FALSE)
set(paramsArgs --T1Tissue 1434
               --T1Blood 1600
               --relaxivity 0.0037
               --S0grad 15.0
               --hematocrit 0.45
               --aucTimeInterval 90
               --fTolerance 1e-4 
               --gTolerance 1e-4 
               --xTolerance 1e-5 
               --epsilon 1e-9 
               --maxIter 200)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
# Regression Tests DROs with T1Map
#-----------------------------------------------------------------------------
set(testName DRO3min5secinf_AllOutputsExceptFpv_WithT1Map)
set(tempOutDataBaseName ${TEMP}/${testName})
set(referenceDataBaseName ${referenceDataBaseDir}DRO3min5secinf_AllOutputsExceptFpv)
set_compareArgs(FALSE)
set(paramsArgs --T1Tissue 1111
               --T1Blood 1600
               --T1Map ${inputDataBaseName}-T1Map.nrrd
               --relaxivity 0.0037
               --S0grad 15.0
               --hematocrit 0.45
               --aucTimeInterval 90
               --fTolerance 1e-4 
               --gTolerance 1e-4 
               --xTolerance 1e-5 
               --epsilon 1e-9 
               --maxIter 200)
set_outputParamsArgs(FALSE)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  ${compareArgs}
  ModuleEntryPoint
    ${paramsArgs}
    ${outputParamsArgs}
    --roiMask ${inputDataBaseName}-ROI.nrrd
    --aifMask ${inputDataBaseName}-AIF.nrrd
    ${inputDataBaseName}3min5secinf.nrrd                   
)
set_property(TEST ${testName} PROPERTY LABELS ${CLP})






