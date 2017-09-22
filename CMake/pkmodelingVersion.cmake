
#.rst:
# PkModelingVersion
# ------------
#
# This module will set these variables in the including scope:
#
# ::
#
#   PkModeling_VERSION_IS_RELEASE
#   PkModeling_VERSION_QUALIFIER
#
# It will also set all variables describing the SCM associated
# with PkModeling_SOURCE_DIR.
#
# It has been designed to be included in the build system of PkModeling.
#
# The following variables are expected to be defined in the including
# scope:
#
# ::
#
#   GIT_EXECUTABLE
#   PkModeling_CMAKE_DIR
#   PkModeling_LATEST_TAG
#   PkModeling_SOURCE_DIR
#
#
# Details description:
#
# .. variable:: PkModeling_VERSION_IS_RELEASE
#
#   Indicate if the current build is release. A particular build is a
#   release if the HEAD of the associated Git checkout corresponds to
#   a tag (either lightweight or annotated).
#
# .. variable:: PkModeling_CMAKE_DIR
#
#   Directory containing all CMake modules included in this module.
#
# .. variable:: PkModeling_LATEST_TAG
#
#   Name of the latest tag. It is used to reference the "bleeding edge"
#   version on GitHub. This lightweight tag is automatically updated
#   by a script running on CI services like Appveyor, CircleCI or TravisCI.
#
# .. variable:: PkModeling_SOURCE_DIR
#
#   PkModeling Git checkout
#
#
# Dependent CMake modules:
#
# ::
#
#   PkModelingMacroExtractRepositoryInfo
#   FindGit.cmake
#

# --------------------------------------------------------------------------
# Sanity checks
# --------------------------------------------------------------------------
set(expected_defined_vars
  GIT_EXECUTABLE
  PkModeling_CMAKE_DIR
  PkModeling_LATEST_TAG
  PkModeling_SOURCE_DIR
  PkModeling_VERSION
  )
foreach(var ${expected_defined_vars})
  if(NOT DEFINED ${var})
    message(FATAL_ERROR "${var} is mandatory")
  endif()
endforeach()

#-----------------------------------------------------------------------------
# Update CMake module path
#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${PkModeling_CMAKE_DIR}
  ${CMAKE_MODULE_PATH}
  )

find_package(Git REQUIRED)

#-----------------------------------------------------------------------------
# PkModeling version number
#-----------------------------------------------------------------------------

set(msg "Checking if building a release")
message(STATUS "${msg}")

include(pkmodelingFunctionExtractGitTagNames)
PkModelingFunctionExtractGitTagNames(
  REF "HEAD"
  SOURCE_DIR ${PkModeling_SOURCE_DIR}
  OUTPUT_VARIABLE TAG_NAMES
  )

# If at least one tag (different from "${PkModeling_LATEST_TAG}") is
# associated with HEAD, we are building a release
list(LENGTH TAG_NAMES tag_count)
list(FIND TAG_NAMES "${PkModeling_LATEST_TAG}" latest_index)
if(
    (tag_count GREATER 0 AND latest_index EQUAL -1) OR
    (tag_count GREATER 1 AND NOT latest_index EQUAL -1)
  )
  set(is_release 1)
  set(is_release_answer "yes (found tags ${TAG_NAMES})")
else()
  set(is_release 0)
  if(NOT latest_index EQUAL -1)
    set(is_release_answer "no (found only '${PkModeling_LATEST_TAG}' tag)")
  else()
    set(is_release_answer "no (found 0 tags)")
  endif()
endif()

message(STATUS "${msg} - ${is_release_answer}")

GIT_WC_INFO(${PkModeling_SOURCE_DIR} PkModeling)

string(REGEX MATCH ".*([0-9][0-9][0-9][0-9])\\-([0-9][0-9])\\-([0-9][0-9]).*" _out
  "${PkModeling_WC_LAST_CHANGED_DATE}")
set(PkModeling_VERSION_DATE "${CMAKE_MATCH_1}${CMAKE_MATCH_2}${CMAKE_MATCH_3}") # YYYYMMDD

if(NOT is_release)
  set(_version_qualifier "-${PkModeling_VERSION_DATE}-${PkModeling_WC_REVISION}")
endif()

set(PkModeling_VERSION_QUALIFIER "${_version_qualifier}")
set(PkModeling_VERSION_IS_RELEASE ${is_release})
