#-----------------------------------------------------------------------------
# Regression Tests QINProstate001
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_SOURCE_DIR}/Data/RegressionTests/QINProstate001/Input/QINProstate001)
set(referenceDataBaseDir ${CMAKE_SOURCE_DIR}/Data/RegressionTests/QINProstate001/Reference/)

#-----------------------------------------------------------------------------
set(testName QINProstate001_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-5
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
  ${tempOutDataBaseName}-diag.nrrd
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
              --concentrations ${tempOutDataBaseName}-conc.nrrd
              --outputKtrans ${tempOutDataBaseName}-ktrans.nrrd 
              --outputVe ${tempOutDataBaseName}-ve.nrrd
              --outputMaxSlope ${tempOutDataBaseName}-maxslope.nrrd
              --outputAUC ${tempOutDataBaseName}-auc.nrrd      
              --outputRSquared ${tempOutDataBaseName}-rsq.nrrd
              --outputBAT ${tempOutDataBaseName}-bat.nrrd           
              --fitted ${tempOutDataBaseName}-fit.nrrd
              --outputDiagnostics ${tempOutDataBaseName}-diag.nrrd
              --roiMask ${inputDataBaseName}-phantom-ROI.nrrd
              --aifMask ${inputDataBaseName}-phantom-AIF.nrrd
              ${inputDataBaseName}-phantom.nrrd                   
              )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName QINProstate001_AllOutputsInclFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-5
  --compare ${referenceDataBaseDir}${testName}-conc.nrrd
  ${tempOutDataBaseName}-conc.nrrd
  --compare ${referenceDataBaseDir}${testName}-ktrans.nrrd
  ${tempOutDataBaseName}-ktrans.nrrd
  --compare ${referenceDataBaseDir}${testName}-ve.nrrd
  ${tempOutDataBaseName}-ve.nrrd
  --compare ${referenceDataBaseDir}${testName}-fpv.nrrd
  ${tempOutDataBaseName}-fpv.nrrd
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
  ${tempOutDataBaseName}-diag.nrrd
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
              --computeFpv
              --concentrations ${tempOutDataBaseName}-conc.nrrd
              --outputKtrans ${tempOutDataBaseName}-ktrans.nrrd 
              --outputVe ${tempOutDataBaseName}-ve.nrrd
              --outputFpv ${tempOutDataBaseName}-fpv.nrrd
              --outputMaxSlope ${tempOutDataBaseName}-maxslope.nrrd
              --outputAUC ${tempOutDataBaseName}-auc.nrrd      
              --outputRSquared ${tempOutDataBaseName}-rsq.nrrd
              --outputBAT ${tempOutDataBaseName}-bat.nrrd           
              --fitted ${tempOutDataBaseName}-fit.nrrd
              --outputDiagnostics ${tempOutDataBaseName}-diag.nrrd
              --roiMask ${inputDataBaseName}-phantom-ROI.nrrd
              --aifMask ${inputDataBaseName}-phantom-AIF.nrrd
              ${inputDataBaseName}-phantom.nrrd                   
              )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
# Regression Tests QINBreast001
#-----------------------------------------------------------------------------
set(inputDataBaseName ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/QINBreast001/Input/QINBreast001)
set(referenceDataBaseDir ${CMAKE_CURRENT_SOURCE_DIR}/../../../Data/RegressionTests/QINBreast001/Reference/)

#-----------------------------------------------------------------------------
set(testName QINBreast001_AllOutputsExceptFpv)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-5
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
  ${tempOutDataBaseName}-diag.nrrd
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
              --usePopAif
              --concentrations ${tempOutDataBaseName}-conc.nrrd
              --outputKtrans ${tempOutDataBaseName}-ktrans.nrrd 
              --outputVe ${tempOutDataBaseName}-ve.nrrd
              --outputMaxSlope ${tempOutDataBaseName}-maxslope.nrrd
              --outputAUC ${tempOutDataBaseName}-auc.nrrd      
              --outputRSquared ${tempOutDataBaseName}-rsq.nrrd
              --outputBAT ${tempOutDataBaseName}-bat.nrrd           
              --fitted ${tempOutDataBaseName}-fit.nrrd
              --outputDiagnostics ${tempOutDataBaseName}-diag.nrrd
              ${inputDataBaseName}-phantom.nrrd                   
              )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName QINBreast001_ConstantBat)
set(tempOutDataBaseName ${TEMP}/${testName})
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-5
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
  ${tempOutDataBaseName}-diag.nrrd
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
              --usePopAif
              --BATCalculationMode UseConstantBAT
              --constantBAT 4
              --concentrations ${tempOutDataBaseName}-conc.nrrd
              --outputKtrans ${tempOutDataBaseName}-ktrans.nrrd 
              --outputVe ${tempOutDataBaseName}-ve.nrrd
              --outputMaxSlope ${tempOutDataBaseName}-maxslope.nrrd
              --outputAUC ${tempOutDataBaseName}-auc.nrrd      
              --outputRSquared ${tempOutDataBaseName}-rsq.nrrd
              --outputBAT ${tempOutDataBaseName}-bat.nrrd           
              --fitted ${tempOutDataBaseName}-fit.nrrd
              --outputDiagnostics ${tempOutDataBaseName}-diag.nrrd
              ${inputDataBaseName}-phantom.nrrd                   
              )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})













