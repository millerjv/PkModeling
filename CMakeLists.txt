cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
set(MODULE_NAME PkModeling)

#-----------------------------------------------------------------------------
set(MODULE_HOMEPAGE "https://github.com/millerjv/PkModeling/wiki/PkModeling")
set(MODULE_CATEGORY "Quantification")
set(MODULE_ICONURL "http://cloud.github.com/downloads/millerjv/PkModeling/PkModeling.png")
set(MODULE_STATUS "Beta")
set(MODULE_CONTRIBUTORS "Yingxuan Zhu (GE), Jim Miller (GE), Ming-ching Chang (GE), Mahnaz Maddah (SRI)")
set(MODULE_DESCRIPTION "PkModeling is a Slicer4 Extension that provides pharmacokinetic modeling for dynamic contrast enhanced MRI (DCE MRI).")
set(MODULE_SCREENSHOTURLS " ")
set(MODULE_ACKNOWLEDGEMENTS "This work is part of the National Alliance for Medical Image Computing (NAMIC), funded by the National Institutes of Health through the NIH Roadmap for Medical Research, Grant U54 EB005149. Information on the National Centers for Biomedical Computing can be obtained from http://nihroadmap.nih.gov/bioinformatics .")
set(MODULE_MAJOR_VERSION 0)
set(MODULE_MINOR_VERSION 1)
set(MODULE_PATCH_VERSION 0)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

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
# Packaging
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

