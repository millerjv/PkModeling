#-----------------------------------------------------------------------------
# Convenience functions to get command snippets needed for most tests
#-----------------------------------------------------------------------------
function(set_compareArgs computeFpv)
  set(compareArgs --compareIntensityTolerance 1e-5
                  --compare ${referenceDataBaseDir}${testName}-conc.nrrd
                  ${tempOutDataBaseName}-conc.nrrd
                  --compare ${referenceDataBaseDir}${testName}-ktrans.nrrd
                  ${tempOutDataBaseName}-ktrans.nrrd
                  --compare ${referenceDataBaseDir}${testName}-ve.nrrd
                  ${tempOutDataBaseName}-ve.nrrd
                  --compare ${referenceDataBaseDir}${testName}-maxslope.nrrd
                  ${tempOutDataBaseName}-maxslope.nrrd
                  --compare ${referenceDataBaseDir}${testName}-auc.nrrd
                  ${tempOutDataBaseName}-auc.nrrd
                  --compare ${referenceDataBaseDir}${testName}-rsq.nrrd
                  ${tempOutDataBaseName}-rsq.nrrd
                  --compare ${referenceDataBaseDir}${testName}-bat.nrrd
                  ${tempOutDataBaseName}-bat.nrrd
                  --compare ${referenceDataBaseDir}${testName}-fit.nrrd
                  ${tempOutDataBaseName}-fit.nrrd
                  --compare ${referenceDataBaseDir}${testName}-diag.nrrd
                  ${tempOutDataBaseName}-diag.nrrd)
  if(${computeFpv})
    set(compareArgs ${compareArgs}
                    --compare ${referenceDataBaseDir}${testName}-fpv.nrrd
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








