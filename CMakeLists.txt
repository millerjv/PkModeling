#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.4)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
#-----------------------------------------------------------------------------
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#-----------------------------------------------------------------------------
# Module name
# - Usually identique to the directory name)
# - You shouldn't call 'project(...)'. See SlicerConfig.cmake for more details.
#-----------------------------------------------------------------------------
set(EXTENSION_NAME PkModeling)

#-----------------------------------------------------------------------------
# Extension metadata
#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "https://github.com/millerjv/PkModeling/wiki/PkModeling")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_ICONURL "")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_CONTRIBUTORS "Yingxuan Zhu (GE), Jim Miller (GE), Ming-ching Chang (GE), Mahnaz Maddah (SRI)")
set(EXTENSION_DESCRIPTION "PkModeling is a Slicer4 Extension that provides pharmacokinetic modeling for dynamic contrast enhanced MRI (DCE MRI).")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_ACKNOWLEDGEMENTS "This work is part of the National Alliance for Medical Image Computing (NAMIC), funded by the National Institutes of Health through the NIH Roadmap for Medical Research, Grant U54 EB005149. Information on the National Centers for Biomedical Computing can be obtained from http://nihroadmap.nih.gov/bioinformatics .")
set(EXTENSION_MAJOR_VERSION 0)
set(EXTENSION_MINOR_VERSION 1)
set(EXTENSION_PATCH_VERSION 0)

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)

  # Additional C/CXX flags - Should be defined before including Slicer_USE_FILE
  set(ADDITIONAL_C_FLAGS "" CACHE STRING "Additional ${EXTENSION_NAME} C Flags")
  set(ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional ${EXTENSION_NAME} CXX Flags")

  include(${Slicer_USE_FILE})

  set(EXTENSION_DEPENDS NA) # Specified as a space separated list or 'NA' if any
  set(EXTENSION_LICENSE_FILE ${Slicer_LICENSE_FILE})
  set(EXTENSION_README_FILE ${Slicer_README_FILE})
  set(EXTENSION_LICENSE_SHORT_DESCRIPTION ${Slicer_LICENSE_SHORT_DESCRIPTION})

  include(SlicerEnableExtensionTesting)
endif()

# Configure XML description
set(SignalIntensitiesToConcentrationValuesCLI_NAME SignalIntensitiesToConcentrationValues)
configure_file(
  ${SignalIntensitiesToConcentrationValuesCLI_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${SignalIntensitiesToConcentrationValuesCLI_NAME}.xml
  @ONLY
  )
  
# Add subdirectories
add_subdirectory(PkSolver)
add_subdirectory(CLI)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)  
endif()

#-----------------------------------------------------------------------------
# Generate extension description file '<EXTENSION_NAME>.s4ext'
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(SlicerFunctionGenerateExtensionDescription)
  slicerFunctionGenerateExtensionDescription(
    EXTENSION_NAME ${EXTENSION_NAME}
    EXTENSION_CATEGORY ${EXTENSION_CATEGORY}
    EXTENSION_ICONURL ${EXTENSION_ICONURL}
    EXTENSION_STATUS ${EXTENSION_STATUS}
    EXTENSION_HOMEPAGE ${EXTENSION_HOMEPAGE}
    EXTENSION_CONTRIBUTORS ${EXTENSION_CONTRIBUTORS}
    EXTENSION_DESCRIPTION ${EXTENSION_DESCRIPTION}
    EXTENSION_SCREENSHOTURLS ${EXTENSION_SCREENSHOTURLS}
    EXTENSION_DEPENDS ${EXTENSION_DEPENDS}
    DESTINATION_DIR ${CMAKE_CURRENT_BINARY_DIR}
    SLICER_WC_REVISION ${Slicer_WC_REVISION}
    SLICER_WC_ROOT ${Slicer_WC_ROOT}
    )
endif()

#-----------------------------------------------------------------------------
# Packaging
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

