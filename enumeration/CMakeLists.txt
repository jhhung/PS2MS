cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall --std=c++17 -Ofast" )

# Boost
find_package ( Boost 1.72 REQUIRED COMPONENTS filesystem system serialization iostreams  )

# RDKIT
find_package( RDKit REQUIRED COMPONENTS
    ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
    RDGeneral SubstructMatch Subgraphs MolDraw2D GraphMol DistGeometry DistGeomHelpers
    MolAlign Optimizer ForceField ForceFieldHelpers Alignment
    ForceField MolTransforms EigenSolvers Descriptors FragCatalog PartialCharges ChemTransforms)
set(RDKit_LIBS
    RDKitChemReactions RDKitFileParsers RDKitSmilesParse RDKitDepictor
    RDKitRDGeometryLib RDKitRDGeneral RDKitSubstructMatch RDKitSubgraphs
    RDKitMolDraw2D RDKitGraphMol RDKitDistGeometry RDKitDistGeomHelpers
    RDKitMolAlign RDKitOptimizer RDKitForceField RDKitForceFieldHelpers
    RDKitAlignment RDKitForceField  RDKitMolTransforms RDKitEigenSolvers RDKitDescriptors)

message(STATUS " Include  - ${RDKit_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}")
message(STATUS " link  - ${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS} ${RDKit_INCLUDE_DIRS} ${CMAKE_CURRENT_LIST_DIR}/graphs)
link_directories(${Boost_LIBRARY_DIRS})

# drug-detection
add_library(planar SHARED
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphIsolator.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphStructure.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphEmbed.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphNonplanar.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphTests.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/stack.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphIO.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/graphPreprocess.c
    ${CMAKE_CURRENT_LIST_DIR}/graphs/listcoll.c)

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)
set(LIBS ${RDKit_THREAD_LIBS} ${Boost_LIBRARIES} z planar ${RDKit_LIBS})


set(SRC_FILES build_database)
foreach(SRC ${SRC_FILES})
    add_executable ( ${SRC} ${SRC}.cpp )
    target_link_libraries ( ${SRC} ${REQUIRED_LIBS} ${LIBS} )
    install ( TARGETS ${SRC}
            DESTINATION ${CMAKE_CURRENT_LIST_DIR}/build )
endforeach(SRC)
          
if (WIN32)
    install ( FILES ${LPSOLVE_INCLUDE_DIR}/lpsolve55.dll
              DESTINATION ${CFM_OUTPUT_DIR} )
endif (WIN32)              
          
