# FindRDKit.cmake
# Placed in the public domain by NextMove Software in 2013
# Try to find RDKit headers and libraries
# Defines:
#
#  RDKIT_FOUND - system has RDKit
#  RDKIT_INCLUDE_DIR - the RDKit include directory
#  RDKIT_INCLUDE_EXT_DIR - the RDKit external directory when including Inchi support
#  RDKIT_LIBRARIES - Link these to use RDKit
#
# References:
#  
#  http://nextmovesoftware.com/blog/2013/02/04/looking-for-a-c-cheminformatics-toolkit/
#  https://github.com/timvdm/MolDB/blob/master/cmake/modules/FindRDKit.cmake

# include(FindPackageHandleStandardArgs)

# if(RDKIT_INCLUDE_DIR AND RDKIT_LIBRARIES)
#   # in cache already or user-specified
#   find_package_handle_standard_args(RDKit  DEFAULT_MSG
#                                   RDKIT_INCLUDE_DIR RDKIT_INCLUDE_EXT_DIR RDKIT_LIBRARIES)
# else()

#   if(NOT RDKIT_INCLUDE_DIR)
#     if(WIN32)
#       find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
#         PATHS
#         ${RDBASE}\\Code
#         $ENV{RDKIT_INCLUDE_DIR}
#         $ENV{RDKIT_INCLUDE_PATH}
#         $ENV{RDKIT_BASE}\\Code
#         $ENV{RDBASE}\\Code
#         C:\\RDKit\\include
#         C:\\RDKit\\Code
#       )
#       find_path(RDKIT_INCLUDE_EXT_DIR INCHI-API/inchi.h
#         PATHS
#         ${RDBASE}\\External
#         $ENV{RDKIT_INCLUDE_EXT_DIR}
#         $ENV{RDKIT_INCLUDE_EXT_PATH}
#         $ENV{RDKIT_BASE}\\External
#         $ENV{RDBASE}\\External
#         C:\\RDKit\\include
#         C:\\RDKit\\External
#       )
#     else()
#       find_path(RDKIT_INCLUDE_DIR GraphMol/RDKitBase.h
#         PATHS
#           ${RDBASE}/Code
#           $ENV{RDKIT_INCLUDE_DIR}
#           $ENV{RDKIT_INCLUDE_PATH}
#           $ENV{RDKIT_BASE}/Code
#           $ENV{RDBASE}/Code
# 	  /usr/include/rdkit
#           /usr/local/include/rdkit
#           /usr/local/rdkit/include/Code
#           /usr/local/rdkit/include
#           /usr/local/rdkit/Code
#           ~/rdkit/Code
#       )
#       find_path(RDKIT_INCLUDE_EXT_DIR INCHI-API/inchi.h
#         PATHS
#           ${RDBASE}/External
#           $ENV{RDKIT_INCLUDE_EXT_DIR}
#           $ENV{RDKIT_INCLUDE_EXT_PATH}
#           $ENV{RDKIT_BASE}/External
#           $ENV{RDBASE}/External
#           /usr/local/rdkit/include/External
#           /usr/local/rdkit/include
#           /usr/local/rdkit/External
#           ~/rdkit/External
#       )      
#     endif()
#     if(RDKIT_INCLUDE_DIR)
#        message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_DIR}")
#     endif()
#     if(RDKIT_INCLUDE_EXT_DIR)
#        message(STATUS "Found RDKit include files at ${RDKIT_INCLUDE_EXT_DIR}")
#     endif()    
#   endif()

#   if(NOT RDKIT_LIBRARIES)
#       find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers
#       PATHS
#         ${RDBASE}/lib
#         $ENV{RDKIT_LIB_DIR}
#         $ENV{RDKIT_LIB_PATH}
#         $ENV{RDKIT_LIBRARIES}
#         $ENV{RDKIT_BASE}/lib
#         $ENV{RDBASE}/lib
#         /usr/local/rdkit/lib
#         ~/rdkit/lib
#         ${RDKIT_LIBRARY_DIR}
#         $ENV{LD_LIBRARY_PATH}

#        #ignore default path, so search starts with above paths
#        NO_DEFAULT_PATH
#     )

#     #run with default paths this time
#     find_library(FILEPARSERS_LIB NAMES FileParsers RDKitFileParsers)

#     if(FILEPARSERS_LIB)
#        GET_FILENAME_COMPONENT(RDKIT_LIBRARY_DIR ${FILEPARSERS_LIB} PATH)
#        message(STATUS "Found RDKit libraries at ${RDKIT_LIBRARY_DIR}")

#       # Note that the order of the following libraries is significant!!
#       find_library(SMILESPARSE_LIB NAMES SmilesParse RDKitSmilesParse
#                                    HINTS ${RDKIT_LIBRARY_DIR})
#       find_library(DEPICTOR_LIB NAMES Depictor RDKitDepictor
#                                 HINTS ${RDKIT_LIBRARY_DIR})
#       find_library(CHEMTRANS_LIB NAMES ChemTransforms RDKitChemTransforms
#                                 HINTS ${RDKIT_LIBRARY_DIR})								
#       find_library(GRAPHMOL_LIB NAMES GraphMol RDKitGraphMol
#                                 HINTS ${RDKIT_LIBRARY_DIR})
#       find_library(RDGEOMETRYLIB_LIB NAMES RDGeometryLib RDKitRDGeometryLib
#                                 HINTS ${RDKIT_LIBRARY_DIR})
#       find_library(RDGENERAL_LIB NAMES RDGeneral RDKitRDGeneral
#                                  HINTS ${RDKIT_LIBRARY_DIR})
#       find_library(GASTEIGER_LIB NAMES PartialCharges RDKitPartialCharges
#                                  HINTS ${RDKIT_LIBRARY_DIR})  
#       find_library(DATASTRUCT_LIB NAMES DataStructs RDKitDataStructs
#                                  HINTS ${RDKIT_LIBRARY_DIR}) 
#       find_library(SUBGRAPH_LIB NAMES Subgraphs RDKitSubgraphs                              
#                                  HINTS ${RDKIT_LIBRARY_DIR})
#       find_library(FINGERPRINT_LIB NAMES Fingerprints RDKitFingerprints
#                                  HINTS ${RDKIT_LIBRARY_DIR})  
#       find_library(INCHI_LIB NAMES Inchi RDKitInchi
#                                  HINTS ${RDKIT_LIBRARY_DIR})      
#       find_library(RDINCHI_LIB NAMES RDInchiLib RDKitRDInchiLib
#                                  HINTS ${RDKIT_LIBRARY_DIR}) 
#       find_library(OPT NAMES Optimizer RDKitOptimizer
# 	                 HINTS ${RDKIT_LIBRARY_DIR}) 								 
#       find_library(FF NAMES ForceField RDKitForceField
# 	                 HINTS ${RDKIT_LIBRARY_DIR})  
#       find_library(FFHELP NAMES ForceFieldHelpers RDKitForceFieldHelpers
# 	                 HINTS ${RDKIT_LIBRARY_DIR}) 
#       find_library(CATALOG NAMES Catalogs RDKitCatalogs RDKit
# 	                 HINTS ${RDKIT_LIBRARY_DIR})  
#       find_library(FRAGCAT NAMES FragCatalog RDKitFragCatalog
# 	                 HINTS ${RDKIT_LIBRARY_DIR})  
#       find_library(SUBSTRUCTMATCH_LIB NAMES SubstructMatch RDKitSubstructMatch
#                                  HINTS ${RDKIT_LIBRARY_DIR})

#       find_library(DATASTRUCTS_LIB NAMES DataStructs RDKitDataStructs
#                                  HINTS ${RDKIT_LIBRARY_DIR})

#       set (RDKIT_LIBRARIES ${FILEPARSERS_LIB} ${SMILESPARSE_LIB} ${SUBSTRUCTMATCH_LIB}
#               ${DEPICTOR_LIB} ${CHEMTRANS_LIB} ${GRAPHMOL_LIB} ${RDGEOMETRYLIB_LIB}
#               ${RDGENERAL_LIB} ${SUBSTRUCT_LIB} ${GASTEIGER_LIB} 
#               ${DATASTRUCT_LIB} ${SUBGRAPH_LIB} ${FINGERPRINT_LIB} 
#               ${INCHI_LIB} ${RDINCHI_LIB} ${OPT} ${FF} ${FFHELP} ${CATALOG} ${FRAGCAT})
#     endif()
#     if(RDKIT_LIBRARIES)
#             message(STATUS "Found RDKit library files at ${RDKIT_LIBRARIES}")
#     endif()
#   endif()

#   find_package_handle_standard_args(RDKit  DEFAULT_MSG
#                                   RDKIT_INCLUDE_DIR RDKIT_INCLUDE_EXT_DIR RDKIT_LIBRARIES)

# mark_as_advanced(RDINCHI_LIB INCHI_LIB GASTEIGER_LIB SUBSTRUCT_LIB RDGENERAL_LIB RDGEOMETRYLIB_LIB GRAPHMOL_LIB DEPICTOR_LIB SMILESPARSE_LIB FILEPARSERS_LIB)
#   mark_as_advanced(RDKIT_INCLUDE_DIR RDKIT_INCLUDE_EXT_DIR RDKIT_LIBRARIES)
# endif()

############

# set(RDKIT_DIR $ENV{RDBASE})
set(RDKIT_DIR ${RDBASE})
if(NOT RDKIT_DIR)
  message (FATAL_ERROR "RDBASE didn't specified" )
endif(NOT RDKIT_DIR)

set(RDKIT_INCLUDE_DIR ${RDKIT_DIR}/Code)
set(RDKIT_INCLUDE_EXT_DIR ${RDKIT_DIR}/External)
if (EXISTS ${RDKIT_DIR}/build/lib)
  set(RDKIT_LIBRARY_DIR ${RDKIT_DIR}/build/lib)
else()
  set(RDKIT_LIBRARY_DIR ${RDKIT_DIR}/lib)
endif()

set(RDKIT_FOUND "MyRDKit_FOUND")
# libraries, as specified in the COMPONENTS
foreach(component ${MyRDKit_FIND_COMPONENTS})
  message( "Looking for RDKit component ${component}" )
  find_file( MyRDKit_LIBRARY_${component}
    libRDKit${component}.so
    PATH ${RDKIT_LIBRARY_DIR} NO_DEFAULT_PATH)
  message("MyRDKit_LIBRARY_${component} : ${MyRDKit_LIBRARY_${component}}")
  if(${MyRDKit_LIBRARY_${component}} MATCHES "-NOTFOUND$")
    message(FATAL_ERROR "Didn't find RDKit ${component} library.")
  endif(${MyRDKit_LIBRARY_${component}} MATCHES "-NOTFOUND$")
  set(RDKIT_LIBRARIES ${RDKIT_LIBRARIES} ${MyRDKit_LIBRARY_${component}})
endforeach(component)

message("RDKIT_FOUND : ${RDKIT_FOUND}")
message("RDKIT_INCLUDE_DIR : ${RDKIT_INCLUDE_DIR}")
message("RDKIT_INCLUDE_EXT_DIR : ${RDKIT_INCLUDE_EXT_DIR}")
message("RDKIT_LIBRARY_DIR : ${RDKIT_LIBRARY_DIR}")
message("RDKIT_LIBRARIES : ${RDKIT_LIBRARIES}")