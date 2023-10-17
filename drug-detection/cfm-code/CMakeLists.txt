##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-code/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
set ( BASE_HEADERS FragmentGraph.h Features.h Config.h Spectrum.h
FragmentGraphGenerator.h FragmentTreeNode.h Inference.h IPFP.h MILP.h Isotope.h
MolData.h MspReader.h Param.h NNParam.h Util.h Message.h Identifier.h Comparators.h
FunctionalGroups.h )

set ( TRAIN_HEADERS Comms.h EM.h EM_NN.h )
 
set ( BASE_SRC_FILES  FragmentGraph.cpp 
                      Features.cpp
                      Config.cpp
                      FragmentGraphGenerator.cpp
                      FragmentTreeNode.cpp
                      FunctionalGroups.cpp
		      Identifier.cpp
                      Inference.cpp
                      IPFP.cpp
		      Isotope.cpp
		      Message.cpp 
                      MILP.cpp
                      MolData.cpp
		      MspReader.cpp
                      Param.cpp
					  NNParam.cpp
		      Spectrum.cpp
                      Util.cpp
		      Comparators.cpp )
                     
set ( TRAIN_SRC_FILES Comms.cpp 
                      EM.cpp
                      EM_NN.cpp )
                      
if (INCLUDE_TRAIN)
    set ( SRC_FILES ${BASE_SRC_FILES} ${TRAIN_SRC_FILES} )
    set ( HEADER_FILES ${BASE_HEADERS} ${TRAIN_HEADERS} )
    set ( REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES}
${LPSOLVE_LIBRARY} ${LBFGS_LIBRARY} ${MPI_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})    
else ()
    set ( SRC_FILES ${BASE_SRC_FILES} )
    set ( HEADER_FILES ${BASE_HEADERS} )
    set ( REQUIRED_LIBS ${RDKIT_LIBRARIES} ${Boost_LIBRARIES}
${LPSOLVE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})       
endif () 

add_library ( cfm-code STATIC ${SRC_FILES} ${HEADER_FILES} )
target_link_libraries ( cfm-code ${REQUIRED_LIBS} )

install ( FILES ${CMAKE_CURRENT_SOURCE_DIR}/ISOTOPE.DAT
          DESTINATION ${CFM_OUTPUT_DIR} )    