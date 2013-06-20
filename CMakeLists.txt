cmake_minimum_required(VERSION 2.8.7)

if (NOT Slicer_SOURCE_DIR)
  set(EXTENSION_NAME PkModeling)
  set(EXTENSION_HOMEPAGE "https://github.com/millerjv/PkModeling/wiki/PkModeling")
  set(EXTENSION_CATEGORY "Quantification")
  set(EXTENSION_CONTRIBUTORS "Yingxuan Zhu (GE), Jim Miller (GE), Andriy Fedorov (BWH), Ming-ching Chang (GE), Mahnaz Maddah (SRI)")
  set(EXTENSION_DESCRIPTION "PkModeling is a Slicer4 Extension that provides pharmacokinetic modeling for dynamic contrast enhanced MRI (DCE MRI).")
  set(EXTENSION_ICONURL "http://wiki.slicer.org/slicerWiki/images/3/34/PkModeling.png")
  set(EXTENSION_STATUS "Beta")
  set(EXTENSION_SCREENSHOTURLS " ")
  set(EXTENSION_ACKNOWLEDGEMENTS "This work is part of the National Alliance for Medical Image Computing (NAMIC), funded by the National Institutes of Health through the NIH Roadmap for Medical Research, Grant U54 EB005149. Information on the National Centers for Biomedical Computing can be obtained from http://nihroadmap.nih.gov/bioinformatics . This work was also funded in part by the Quantitative Imaging Network (QIN) Grant U01CA151261.")
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
# PkModeling has a support library that needs to be installed as part
# of the packaging
#
# Disabling.  Changed the PkSolver library to be build statiically so it does not need to be installed
# option(PkModeling_BUILD_SLICER_EXTENSION "Build as a Slicer extension" OFF)
#
# set(PkModeling_INSTALL_SUPPORT_RUNTIME_DESTINATION bin)
# set(PkModeling_INSTALL_SUPPORT_LIBRARY_DESTINATION lib)
# set(PkModeling_INSTALL_SUPPORT_ARCHIVE_DESTINATION lib)
# if(PkModeling_BUILD_SLICER_EXTENSION)
#  set(PkModeling_INSTALL_SUPPORT_RUNTIME_DESTINATION ${Slicer_INSTALL_BIN_DIR})
#  set(PkModeling_INSTALL_SUPPORT_LIBRARY_DESTINATION ${Slicer_INSTALL_LIB_DIR})
#  set(PkModeling_INSTALL_SUPPORT_ARCHIVE_DESTINATION ${Slicer_INSTALL_LIB_DIR})
# endif()

#-----------------------------------------------------------------------------
add_subdirectory(PkSolver)
add_subdirectory(CLI)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

