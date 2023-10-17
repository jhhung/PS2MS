##########################################################################

# Mass Spec Prediction and Identification of Metabolites

#

# CMakeLists.txt

#

# Copyright (c) 2013, Felicity Allen

# All rights reserved.

#

# This file is part of the cfm-id project.

# The contents are covered by the terms of the GNU Lesser General Public

# License, which is included in the file license.txt, found at the root

# of the cfm source tree.

#########################################################################*/

cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall --std=c++17 -Ofast" )


project (cfm)



set (CFM_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE STRING "Where do you want the output binaries to install to?" )



#It is recommended, but not strictly necessary, to build and run the tests

set( INCLUDE_TESTS OFF CACHE BOOL "Do you want to build the tests?" )



#Building the training code requires extra dependencies, e.g. MPI

set( INCLUDE_TRAIN OFF CACHE BOOL "Do you want to build the training code?" )



# BOOST
SET (BOOST_ROOT "/mnt/ec/ness/brian/spectrum/boost_1_72_0")
SET (BOOST_INCLUDEDIR "/mnt/ec/ness/brian/spectrum/boost_1_72_0/include")
SET (BOOST_LIBRARYDIR "/mnt/ec/ness/brian/spectrum/boost_1_72_0/lib")
set (Boost_NO_BOOST_CMAKE ON)

set (Boost_ADDITIONAL_VERSIONS "1.72")

set (Boost_USE_STATIC_LIBS ON)

include (FindBoost)

find_package ( Boost 1.72 COMPONENTS filesystem system serialization iostreams REQUIRED )

include_directories( ${Boost_INCLUDE_DIRS} )

link_directories ( ${Boost_LIBRARY_DIR} )

message(STATUS " Includes - ${Boost_INCLUDE_DIRS}")
message(STATUS " Library  - ${Boost_LIBRARY_DIRS}")

# RDKit

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} )

find_package( MyRDKit COMPONENTS ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs MolDraw2D GraphMol DistGeometry DistGeomHelpers
  MolAlign Optimizer ForceField ForceFieldHelpers Alignment
  ForceField MolTransforms EigenSolvers Descriptors FragCatalog PartialCharges RDInchiLib ChemTransforms )

if (RDKIT_FOUND)

    include_directories( ${RDKIT_INCLUDE_DIR} ${RDKIT_INCLUDE_EXT_DIR} )

    link_directories ( ${RDKIT_LIBRARY_DIR} )

else ()

    message (FATAL_ERROR "Could not find RDKit" )

endif ()



# LPSolve

set ( LPSOLVE_INCLUDE_DIR "" CACHE STRING "Where are the LPSolve headers?" )

set ( LPSOLVE_LIBRARY_DIR "" CACHE STRING "Where are the LPSolve libraries?" )

include_directories ( ${LPSOLVE_INCLUDE_DIR} )

link_directories ( ${LPSOLVE_LIBRARY_DIR} )

find_library(LPSOLVE_LIBRARY NAMES lpsolve55 PATHS ${LPSOLVE_LIBRARY_DIR})

if (NOT LPSOLVE_LIBRARY)

  MESSAGE(FATAL_ERROR "Cannot find LPSOLVE libraries!")

endif (NOT LPSOLVE_LIBRARY)

mark_as_advanced( LPSOLVE_LIBRARY )



find_package( Threads )



#  MPI

if ( INCLUDE_TRAIN OR INCLUDE_TESTS )

    include (FindMPI)

    if ( MPI_FOUND )

        include_directories ( ${MPI_INCLUDE_PATH} )

    else ()

        message (FATAL_ERROR "Could not find MPI" )

    endif ()

endif ()



# LBFGS

if ( INCLUDE_TRAIN OR INCLUDE_TESTS )

	set ( LBFGS_INCLUDE_DIR "" CACHE STRING "Where are the lbfgs headers?" )

	set ( LBFGS_LIBRARY_DIR "" CACHE STRING "Where are the lbfgs libraries?" )

	include_directories ( ${LBFGS_INCLUDE_DIR} )

	link_directories ( ${LBFGS_LIBRARY_DIR} )

	find_library(LBFGS_LIBRARY NAMES lbfgs PATHS ${LBFGS_LIBRARY_DIR})

	if (NOT LBFGS_LIBRARY)

	  MESSAGE(FATAL_ERROR "Cannot find LBFGS libraries!")

	endif (NOT LBFGS_LIBRARY)

	mark_as_advanced( LBFGS_LIBRARY )

endif ()



#Common Library

add_subdirectory( cfm-code )



#Executable for generating fragmentation graphs

include_directories( ${CMAKE_MODULE_PATH}/cfm-code )

# add_subdirectory( fraggraph-gen )



#Executable for running tests

# if (INCLUDE_TESTS)

#     add_subdirectory( cfm-test )

# endif ()



# if (INCLUDE_TRAIN)

#     add_subdirectory( cfm-train )

# endif ()      



#Spectrum Predicter

# add_subdirectory( cfm-predict )



#Identifier

# add_subdirectory( cfm-id )



#Annotater

# add_subdirectory( cfm-annotate )



#Compute Statistics (for computed spectra)

# add_subdirectory( compute-stats )



# add_subdirectory( cfm-id-precomputed )

#drug-detection

add_subdirectory( drug-comparation )
