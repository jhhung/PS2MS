##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-id/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( SRC_FILES  main.cpp )
                                           
add_executable ( cfm-id-precomputed ${SRC_FILES} )
target_link_libraries ( cfm-id-precomputed cfm-code ${REQUIRED_LIBS} )

install ( TARGETS cfm-id-precomputed 
          DESTINATION ${CFM_OUTPUT_DIR} )
          
if (WIN32)
    install ( FILES ${LPSOLVE_INCLUDE_DIR}/lpsolve55.dll
              DESTINATION ${CFM_OUTPUT_DIR} )
endif (WIN32)              
          