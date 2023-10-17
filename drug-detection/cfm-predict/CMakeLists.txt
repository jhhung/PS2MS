##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-predict/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( SRC_FILES  test.cpp )
                                           
add_executable ( compare_test ${SRC_FILES} )
target_link_libraries ( compare_test cfm-code ${REQUIRED_LIBS} )

install ( TARGETS compare_test
          DESTINATION ${CFM_OUTPUT_DIR} )

set ( SRC_FILES  main.cpp )
                                           
add_executable ( cfm-predict ${SRC_FILES} )
target_link_libraries ( cfm-predict cfm-code ${REQUIRED_LIBS} )

install ( TARGETS cfm-predict
          DESTINATION ${CFM_OUTPUT_DIR} )
          
if (WIN32)
    install ( FILES ${LPSOLVE_INCLUDE_DIR}/lpsolve55.dll
              DESTINATION ${CFM_OUTPUT_DIR} )
endif (WIN32)              
          