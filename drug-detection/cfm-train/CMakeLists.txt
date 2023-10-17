##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( SRC_FILES  main.cpp )
                                           
add_executable ( cfm-train ${SRC_FILES} )
target_link_libraries ( cfm-train cfm-code ${REQUIRED_LIBS} )

install ( TARGETS cfm-train 
          DESTINATION ${CFM_OUTPUT_DIR} )
          
          
          