cmake_minimum_required(VERSION 2.8.7)

if (NOT Slicer_SOURCE_DIR)
  set(EXTENSION_NAME PkModeling)
  set(EXTENSION_HOMEPAGE "https://github.com/millerjv/PkModeling/wiki/PkModeling")
  set(EXTENSION_CATEGORY "Quantification")
  set(EXTENSION_CONTRIBUTORS "Yingxuan Zhu (GE), Jim Miller (GE), Ming-ching Chang (GE), Mahnaz Maddah (SRI)")
  set(EXTENSION_DESCRIPTION "PkModeling is a Slicer4 Extension that provides pharmacokinetic modeling for dynamic contrast enhanced MRI (DCE MRI).")
  set(EXTENSION_ICONURL "http://wiki.slicer.org/slicerWiki/images/3/34/PkModeling.png")
  set(EXTENSION_STATUS "Beta")
  set(EXTENSION_SCREENSHOTURLS " ")
  set(EXTENSION_ACKNOWLEDGEMENTS "This work is part of the National Alliance for Medical Image Computing (NAMIC), funded by the National Institutes of Health through the NIH Roadmap for Medical Research, Grant U54 EB005149. Information on the National Centers for Biomedical Computing can be obtained from http://nihroadmap.nih.gov/bioinformatics .")
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
add_subdirectory(PkSolver)
add_subdirectory(CLI)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

