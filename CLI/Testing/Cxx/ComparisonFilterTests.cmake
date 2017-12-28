#-----------------------------------------------------------------------------
# Unit/Integration Tests ComparisonFilter
#-----------------------------------------------------------------------------
set(referenceDataBaseDir ${CMAKE_SOURCE_DIR}/Data/TestData/MiniVolumes/)

#-----------------------------------------------------------------------------
set(testName ComparisonFilter_BasicEqual)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-5
  --compare ${referenceDataBaseDir}TestBasic01.nrrd
  ${referenceDataBaseDir}TestBasic01.nrrd
  DoNothingAndPass
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName ComparisonFilter_BasicDiff)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --expectFail
  --compareIntensityTolerance 1e-5
  --compare ${referenceDataBaseDir}TestBasic01.nrrd
  ${referenceDataBaseDir}TestBasic02.nrrd
  DoNothingAndPass
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName ComparisonFilter_VecImgEqual)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --compareIntensityTolerance 1e-5
  --compare ${referenceDataBaseDir}TestConc01.nrrd
  ${referenceDataBaseDir}TestConc01.nrrd
  DoNothingAndPass
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})

#-----------------------------------------------------------------------------
set(testName ComparisonFilter_VecImgDiff)
add_test(NAME ${testName} COMMAND ${Launcher_Command} $<TARGET_FILE:${CLP}Test>
  --expectFail
  --compareIntensityTolerance 1e-5
  --compare ${referenceDataBaseDir}TestConc01.nrrd
  ${referenceDataBaseDir}TestConc02.nrrd
  DoNothingAndPass
  )
set_property(TEST ${testName} PROPERTY LABELS ${CLP})







