##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# compute-stats/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( SRC_FILES  main.cpp )
                                           
add_executable ( compute-stats ${SRC_FILES} )
target_link_libraries ( compute-stats cfm-code ${REQUIRED_LIBS} )

install ( TARGETS compute-stats 
          DESTINATION ${CFM_OUTPUT_DIR} )
          
          
          