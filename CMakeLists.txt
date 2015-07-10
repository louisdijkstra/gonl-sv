# Copyright (C) 2015 Louis Dijkstra
# 
# This file is part of gonl-sv
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


cmake_minimum_required (VERSION 2.6)
project (gonl_sv)

# The version number: 
set (gonl_tag_deletions_VERSION_MAJOR 1)
set (gonl_tag_deletions_VERSION_MINOR 0)

set (gonl_create_pairs_VERSION_MAJOR 1)
set (gonl_create_pairs_VERSION_MINOR 0)

set (gonl_imputation_VERSION_MAJOR 1)
set (gonl_imputation_VERSION_MINOR 0)

# Get the source code from src
file(GLOB SOURCES src/*_main.c)

# Set the output directory for the build executables
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/src")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}) 

FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

FIND_PACKAGE(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})
set(LIBS ${LIBS} ${GMP_LIBRARIES})

file(GLOB LIB_TAG_SOURCES src/haplotypes.c src/tagsnps_input.c)
file(GLOB LIB_TAG_HEADERS src/haplotypes.h src/tagsnps_input.h)
add_library( tag_library ${LIB_TAG_SOURCES} ${LIB_TAG_HEADERS})

file(GLOB LIB_PAIRS_SOURCES src/haplotypes.c src/snp_deletion_pair_input.c)
file(GLOB LIB_PAIRS_HEADERS src/haplotypes.h src/snp_deletion_pair_input.h)
add_library( pairs_library ${LIB_PAIRS_SOURCES} ${LIB_PAIRS_HEADERS})

file(GLOB LIB_IMPUTATION_SOURCES src/haplotypes.c src/imputation_input.c)
file(GLOB LIB_IMPUTATION_HEADERS src/haplotypes.h src/imputation_input.h)
add_library( imputation_library ${LIB_IMPUTATION_SOURCES} ${LIB_IMPUTATION_HEADERS})

add_executable(gonl_tag_deletions src/tagsnps_main.c)
target_link_libraries(gonl_tag_deletions tag_library ${LIBS})

add_executable(gonl_create_pairs src/snp_deletion_pair_main.c)
target_link_libraries(gonl_create_pairs pairs_library ${LIBS})

add_executable(gonl_imputation src/imputation_main.c)
target_link_libraries(gonl_imputation imputation_library ${LIBS})

#install (TARGETS sm_caller DESTINATION "${PROJECT_SOURCE_DIR}/bin")
install (FILES "${PROJECT_BINARY_DIR}/config.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
install (FILES "src/snp_deletion_pair_input.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
install (FILES "src/tagsnps_input.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
install (FILES "src/imputation_input.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
install (FILES "src/haplotypes.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
