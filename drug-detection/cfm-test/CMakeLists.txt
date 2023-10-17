##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-test/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( SRC_FILES  test.cpp 
                 tests/comms_test.cpp tests/comms_test.h 
                 tests/em_test.cpp tests/em_test.h 
                 tests/features_test.cpp tests/features_test.h
				 tests/fraggen_test.cpp tests/fraggen_test.h
				 tests/inference_test.cpp tests/inference_test.h   
                 tests/ipfp_test.cpp tests/ipfp_test.h 
                 tests/message_test.cpp tests/message_test.h 
				 tests/msp_reader_test.cpp tests/msp_reader_test.h 
                 tests/param_test.cpp tests/param_test.h
				 tests/nn_param_test.cpp tests/nn_param_test.h
				 tests/spectrum_test.cpp tests/spectrum_test.h
				 tests/comparator_test.cpp tests/comparator_test.h )
                                           
add_executable ( cfm-test ${SRC_FILES} )
target_link_libraries ( cfm-test cfm-code ${REQUIRED_LIBS} )

install ( TARGETS cfm-test 
          DESTINATION ${CFM_OUTPUT_DIR} )
          
install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_data
          DESTINATION ${CFM_OUTPUT_DIR}/tests )          
          
          
          