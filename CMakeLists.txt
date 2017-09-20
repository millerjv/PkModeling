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


set(PkModeling_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

#-----------------------------------------------------------------------------
# Superbuild
#
set(CMAKE_MODULE_PATH "${PkModeling_CMAKE_DIR}" ${CMAKE_MODULE_PATH})
set(EXTERNAL_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CMakeExternals)
set(EXTERNAL_PROJECT_FILE_PREFIX "")
include(ExternalProject)
include(ExternalProjectDependency)

#-----------------------------------------------------------------------------
if(APPLE)
  # Note: By setting CMAKE_OSX_* variables before any enable_language() or project() calls,
  #       we ensure that the bitness will be properly detected.
  include(PkModelingBlockSetCMakeOSXVariables)
  mark_as_superbuild(
    VARS CMAKE_OSX_ARCHITECTURES:STRING CMAKE_OSX_SYSROOT:PATH CMAKE_OSX_DEPLOYMENT_TARGET:STRING
    ALL_PROJECTS
    )
endif()

project(PkModeling)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LC)


#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Configuring with build type '${CMAKE_BUILD_TYPE}'")
  mark_as_superbuild(VARS CMAKE_BUILD_TYPE ALL_PROJECTS)
endif()

mark_as_superbuild(
  VARS
    CMAKE_PREFIX_PATH:STRING
    CMAKE_DEBUG_POSTFIX:STRING
  ALL_PROJECTS
  )

#-----------------------------------------------------------------------------
# Option default values

set(superbuild_default ON)
if(Slicer_SOURCE_DIR)
  # Default values used when PkModeling is integrated into Slicer as Remote project.
  set(superbuild_default OFF)
endif()

#-----------------------------------------------------------------------------
# Options

include(CMakeDependentOption)

option(PkModeling_BUILD_APPS "Build ${PROJECT_NAME} applications." ON)
mark_as_superbuild(PkModeling_BUILD_APPS)

#option(DCMQI_BUILD_DOC "Build ${PROJECT_NAME} documentation." ${build_doc_default})
#mark_as_superbuild(DCMQI_BUILD_DOC)


#-----------------------------------------------------------------------------
# Standalone vs Slicer extension option
#

# This option should be named after the project name, it corresponds to the
# option set to ON when the project is build by the Slicer Extension build
# system.

set(_default OFF)
set(_reason "PkModeling_BUILD_SLICER_EXTENSION is ON")
if(NOT DEFINED PkModeling_BUILD_SLICER_EXTENSION AND DEFINED Slicer_DIR)
  set(_default ON)
  set(_reason "Slicer_DIR is SET")
endif()

option(PkModeling_BUILD_SLICER_EXTENSION "Build as a Slicer Extension" ${_default})

set(_msg "Checking if building as a Slicer extension")
message(STATUS ${_msg})
if(PkModeling_BUILD_SLICER_EXTENSION)
  message(STATUS "${_msg} - yes (${_reason})")
else()
  message(STATUS "${_msg} - no (PkModeling_BUILD_SLICER_EXTENSION is OFF)")
endif()

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#
option(PkModeling_SUPERBUILD "Build ${PROJECT_NAME} and the projects it depends on." ${superbuild_default})
mark_as_advanced(PkModeling_SUPERBUILD)

#-----------------------------------------------------------------------------
# Output directories.
#
foreach(type LIBRARY RUNTIME ARCHIVE)
  # Make sure the directory exists
  if(DEFINED PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY
     AND NOT EXISTS ${PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY})
    message(FATAL_ERROR "PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY is set to a non-existing directory [${PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY}]")
  endif()

  if(PkModeling_SUPERBUILD)
    set(output_dir ${PkModeling_BINARY_DIR}/bin)
    if(NOT DEFINED PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY)
      set(PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY ${PkModeling_BINARY_DIR}/${PROJECT_NAME_LC}-build/bin)
    endif()
    mark_as_superbuild(PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY:PATH)
  else()
    if(NOT DEFINED PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY)
      set(output_dir ${PkModeling_BINARY_DIR}/bin)
    else()
      set(output_dir ${PkModeling_CMAKE_${type}_OUTPUT_DIRECTORY})
    endif()
  endif()
  set(CMAKE_${type}_OUTPUT_DIRECTORY ${output_dir} CACHE INTERNAL "Single output directory for building all libraries.")
endforeach()

#-----------------------------------------------------------------------------
# PkModeling version number
#
set(PkModeling_VERSION_MAJOR 1)
set(PkModeling_VERSION_MINOR 0)
set(PkModeling_VERSION_PATCH 5)
set(PkModeling_VERSION
    "${PkModeling_VERSION_MAJOR}.${PkModeling_VERSION_MINOR}.${PkModeling_VERSION_PATCH}")

set(PkModeling_LATEST_TAG "latest")

# Append the library version information to the library target
# properties.  A parent project may set its own properties and/or may
# block this.
if(NOT PkModeling_NO_LIBRARY_VERSION)
  set(PkModeling_LIBRARY_PROPERTIES ${PkModeling_LIBRARY_PROPERTIES}
    VERSION "${PkModeling_VERSION}"
    SOVERSION "${PkModeling_VERSION_MAJOR}.${PkModeling_VERSION_MINOR}"
    )
endif()

#-----------------------------------------------------------------------------
# Install directories, used for install rules.
#
if(NOT PkModeling_INSTALL_BIN_DIR)
  set(PkModeling_INSTALL_BIN_DIR "bin")
endif()
if(NOT PkModeling_INSTALL_LIB_DIR)
  set(PkModeling_INSTALL_LIB_DIR "lib/PkModeling-${PkModeling_VERSION_MAJOR}.${PkModeling_VERSION_MINOR}")
endif()
if(NOT PkModeling_INSTALL_CMAKE_DIR)
  set(PkModeling_INSTALL_CMAKE_DIR "lib/PkModeling-${PkModeling_VERSION_MAJOR}.${PkModeling_VERSION_MINOR}/CMake")
endif()
if(NOT PkModeling_INSTALL_INCLUDE_DIR)
  set(PkModeling_INSTALL_INCLUDE_DIR "include/PkModeling-${PkModeling_VERSION_MAJOR}.${PkModeling_VERSION_MINOR}")
endif()

mark_as_superbuild(
  PkModeling_INSTALL_BIN_DIR:STRING
  PkModeling_INSTALL_LIB_DIR:STRING
  PkModeling_INSTALL_CMAKE_DIR:STRING
  PkModeling_INSTALL_INCLUDE_DIR:STRING
  )

#-----------------------------------------------------------------------------
# Other options

# Git protocol option
option(PkModeling_USE_GIT_PROTOCOL "If behind a firewall turn this OFF to use http instead." ON)
mark_as_advanced(PkModeling_USE_GIT_PROTOCOL)
set(git_protocol "git")
if(NOT PkModeling_USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

#-----------------------------------------------------------------------------
# Dependencies
#

# Slicer
if(PkModeling_BUILD_SLICER_EXTENSION)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  mark_as_superbuild(Slicer_DIR)
endif()


#-----------------------------------------------------------------------------
# Superbuild script
#

if(PkModeling_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Dependencies
#

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(ITK NO_MODULE REQUIRED)
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
# Testing
#
include(CTest)
set(TEMP_DIR ${CMAKE_BINARY_DIR}/Testing/Temporary)

#-----------------------------------------------------------------------------
# PkModeling src code
#
add_subdirectory(PkSolver)
add_subdirectory(CLI)

#-----------------------------------------------------------------------------
# PkModelingTargets
#

if(NOT Slicer_DIR)

  get_property(PkModeling_TARGETS GLOBAL PROPERTY PkModeling_TARGETS)
  export(TARGETS ${PkModeling_TARGETS} FILE ${CMAKE_BINARY_DIR}/PkModelingTargets.cmake)

endif()

#-----------------------------------------------------------------------------
# Version information
#
# Add target extracting version information. The target configures
# 'pkmodelingCPackOptions.cmake' and 'pkmodelingVersionConfigure.h'
#
include(pkmodelingConfigureVersionTarget)
add_dependencies(PkModeling pkmodelingConfigureVersion)


#-----------------------------------------------------------------------------
# PkModelingConfig
#

if(NOT Slicer_DIR)

  # Configure 'pkmodelingConfig.cmake'
  include(CMakePackageConfigHelpers)

  # For a build tree
  set(CONFIG_DIR_CONFIG ${CMAKE_BINARY_DIR})
  set(DCMTK_DIR_CONFIG ${DCMTK_DIR})
  set(SlicerExecutionModel_DIR_CONFIG ${SlicerExecutionModel_DIR})
  set(${PROJECT_NAME}_CONFIG_CODE "####### Expanded from \@DCMQI_CONFIG_CODE\@ #######\n")
  set(${PROJECT_NAME}_CONFIG_CODE "find_package(DCMTK NO_MODULE REQUIRED)\n")
  set(${PROJECT_NAME}_CONFIG_CODE "${${PROJECT_NAME}_CONFIG_CODE}##################################################")
  set(_config ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
  configure_package_config_file(
      CMake/${PROJECT_NAME}Config.cmake.in
      ${_config}
      INSTALL_DESTINATION ${CMAKE_BINARY_DIR}
      PATH_VARS
        CONFIG_DIR_CONFIG
        DCMTK_DIR_CONFIG
        SlicerExecutionModel_DIR_CONFIG
      NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )

endif()

#-----------------------------------------------------------------------------
# Packaging
#

if(DEFINED Slicer_DIR)
  include(${Slicer_EXTENSION_CPACK})
else()
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${PROJECT_NAME};RuntimeLibraries;/")
  if(EXISTS "${PkModeling_DIR}/cmake_install.cmake")
    set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${DCMTK_DIR};DCMTK;shlib;/")
  endif()
  if(EXISTS "${ITK_DIR}/cmake_install.cmake")
    set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${ITK_DIR};ITK;RuntimeLibraries;/")
  endif()
  set(CPACK_PACKAGE_NAME "PkModeling")
  set(CPACK_PACKAGE_VENDOR "Quantitative Image Informatics for Cancer Research (QIICR)")
  set(CPACK_PACKAGE_VERSION_MAJOR ${PkModeling_VERSION_MAJOR})
  set(CPACK_PACKAGE_VERSION_MINOR ${PkModeling_VERSION_MINOR})
  set(CPACK_PACKAGE_VERSION_PATCH ${PkModeling_VERSION_PATCH})
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(CPACK_SYSTEM_NAME "mac")
  elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CPACK_SYSTEM_NAME "linux")
  elseif(CMAKE_SYSTEM_NAME MATCHES "Windows")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(CPACK_SYSTEM_NAME win64)
    else()
      set(CPACK_SYSTEM_NAME win32)
    endif()
  endif()
  set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${PkModeling_VERSION}-${CPACK_SYSTEM_NAME}")
  set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${EXTENSION_DESCRIPTION})
  if(UNIX)
    set(CPACK_GENERATOR "TGZ")
  elseif(WIN32)
    set(CPACK_GENERATOR "ZIP")
  endif()
  if(NOT DEFINED PkModeling_CPACK_PROJECT_CONFIG_FILE)
    message(FATAL_ERROR "PkModeling_CPACK_PROJECT_CONFIG_FILE is not defined")
  endif()
  set(CPACK_PROJECT_CONFIG_FILE "${PkModeling_CPACK_PROJECT_CONFIG_FILE}")
  include(CPack)
endif()
