##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# cfm-predict/CMakeLists.txt
#
# Author: Felicity Allen
# Created: August 2013
##########################################################################
 
#add_library(planar SHARED graphs/graphIsolator.c graphs/graphStructure.c graphs/graphEmbed.c graphs/graphNonplanar.c graphs/graphTests.c graphs/stack.c graphs/graphIO.c graphs/graphPreprocess.c graphs/listcoll.c)

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} ${RDKit_THREAD_LIBS} z /mnt/ec/ness/brian/spectrum/filter/cpp/libplanar.so )

include_directories(${RDKIT_INCLUDE_DIR} graphs)
LINK_directories("/usr/local/lib")

# set(SRC_FILES build_database fetch_database cfm_validate indicate_modification filt_spectrum build_from_smiles predict_smiles_peak)
set(SRC_FILES cfm_comparation cfm_comparation_without_answer)
foreach(SRC ${SRC_FILES})
    add_executable ( ${SRC} ${SRC}.cpp )
    target_link_libraries ( ${SRC} cfm-code ${REQUIRED_LIBS} ${LIBS} profiler ) #tcmalloc

    install ( TARGETS ${SRC}
            DESTINATION ${CFM_OUTPUT_DIR} )
endforeach(SRC)
          
if (WIN32)
    install ( FILES ${LPSOLVE_INCLUDE_DIR}/lpsolve55.dll
              DESTINATION ${CFM_OUTPUT_DIR} )
endif (WIN32)              
          
