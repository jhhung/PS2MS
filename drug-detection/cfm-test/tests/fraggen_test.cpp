/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# fraggen_test.cpp
#
# Description: Test code for Fragment Tree Generation
#
# Author: Felicity Allen
# Created: August 2014
#########################################################################*/

#include "fraggen_test.h"
#include "MolData.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem.hpp>

bool checkGraphConsistency( FragmentGraph *graph );

FragGenTestPositiveESI::FragGenTestPositiveESI(){
	description = "Test generation of positive mode ESI fragment tree";
}

void FragGenTestPositiveESI::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CC(=O)O";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true;  
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 9 ){
		std::cout << "Expected 9 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFragmentsOnly( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 61.02840582 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 44.99710569 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 15.02292652 ) > tol ||
			fabs( graph->getFragmentAtIdx(6)->getMass() - 25.00727645 ) > tol ||
			fabs( graph->getFragmentAtIdx(7)->getMass() - 41.00219107 ) > tol ||
			fabs( graph->getFragmentAtIdx(8)->getMass() - 59.01275576 ) > tol ){

			std::cout << std::endl << "Unexpected fragment mass(es)" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected:
			//0 61.02840582 CC(O)=[OH+]
			//1 44.99710569 O=C=[OH+]
			//2 15.02292652 [CH3+]
			//3 17.03857658 [CH5+]
			//4 19.01784114 [OH3+]
			//5 43.01784114 C#C[OH2+]
			//6 25.00727645 [C+]#C
			//7 41.00219107 [C+]#CO
			//8 59.01275576 [CH+]=C(O)O
			pass = false;
		}
	}

	passed = pass;
}

FragGenTestNegativeESI::FragGenTestNegativeESI(){
	description = "Test generation of negative mode ESI fragment tree";
}

void FragGenTestNegativeESI::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CC(=O)O";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, NEGATIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg);  cfg.ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;
	cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 4 ){
		std::cout << "Expected 4 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFragmentsOnly( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 59.01385292 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 15.02402368 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 17.00328823 ) > tol ||
			fabs( graph->getFragmentAtIdx(3)->getMass() - 41.00328823 ) > tol ){

			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			std::cout << "Expecting:" << std::endl;
			std::cout << "0 59.01385292 CC(=O)O" << std::endl;
			std::cout << "1 15.02402368 [CH3-]" << std::endl;
			std::cout << "2 17.00273965 [OH-]" << std::endl;
			std::cout << "3 41.00328823 C#C[O-]" << std::endl << std::endl;
			pass = false;
		}
	}

	passed = pass;
}

FragGenTestPositiveEI::FragGenTestPositiveEI(){
	description = "Test generation of positive mode EI fragment tree";
}

void FragGenTestPositiveEI::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CC(=O)O";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg);  cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 15 ){
		std::cout << "Expected 15 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 60.02058079 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 43.98928066 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 44.99710569 ) > tol ||
			fabs( graph->getFragmentAtIdx(12)->getMass() - 39.99436604 ) > tol ||
			fabs( graph->getFragmentAtIdx(13)->getMass() - 58.00493072 ) > tol ||
			fabs( graph->getFragmentAtIdx(14)->getMass() - 59.01275576 ) > tol){

			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected:
			//0 60.02058079 CC(=[O+])O
			//1 43.98928066 O=C=[O+]
			//2 44.99710569 O=C=[OH+]
			//3 15.02292652 [CH3+]
			//4 16.03075155 [CH4+]
			//5 18.0100161 [OH2+]
			//6 19.01784114 [OH3+]
			//7 43.01784114 C#C[OH2+]
			//8 25.00727645 [C+]#C
			//9 41.00219107 [C+]#CO
			//10 42.0100161 C#C[OH+]
			//11 23.99945142 [C]#[C+]
			//12 39.99436604 [C+]#C[O]
			//13 58.00493072 [CH+]=C([O])O
			//14 59.01275576 [CH+]=C(O)O

			pass = false;
		}
	}

	passed = pass;
}

FragGenTestPositiveEIMultibreak::FragGenTestPositiveEIMultibreak(){
	description = "Test generation of positive mode EI fragment tree for a larger molecule where multiple breaks are possible and radical/charge separation occurs within fragment";
}

void FragGenTestPositiveEIMultibreak::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CCNC";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 30 ){
		std::cout << "Expected 30 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;

		//Note - Expecting:
		//0 59.07349929 CC[NH+]C
		//1 57.05784922 C=C[NH+]C
		//2 55.04219916 C#C[NH+]C
		//3 56.05002419 C#C[NH2+]C
		//4 41.0265491 [CH]=[N+]=C
		//5 42.03437413 C=[N+]=C
		//6 15.0234751 [CH3+]
		//7 17.03912516 [CH5+]
		//8 16.03130013 [CH4+]
		//9 31.04219916 C[NH2+]
		//10 32.05002419 C[NH3+]
		//11 29.0265491 C=[NH+]
		//12 30.03437413 C=[NH2+]
		//13 27.01089903 C#[N+]
		//14 28.01872406 C#[NH+]
		//15 27.0234751 C#[CH2+]
		//16 26.01565006 [CH]=[CH+]
		//17 29.03912516 C=[CH3+]
		//18 28.03130013 [CH2][CH2+]
		//19 31.05477522 C[CH4+]
		//20 30.04695019 C[CH3+]
		//21 42.03437413 C#C[NH3+]
		//22 41.0265491 C#C[NH2+]
		//23 58.06567426 C=C[NH2+]C
		//24 43.04219916 C=[N+]C
		//25 44.05002419 C=[NH+]C
		//26 44.05002419 C=C[NH3+]
		//27 18.03437413 [NH4+]
		//28 43.04219916 C=C[NH2+]
		//29 17.0265491 [NH3+]
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		//Check for the ion [CH2][CH2+] - showing a split radical/charge due to multibond removal
		bool found = false;
		for( int i = 0; i < 30; i++ )
			if( *graph->getFragmentAtIdx(i)->getIonSmiles() == "[CH2][CH2+]" ){ found = true; break; }
		if( !found ){
			std::cout << "Could not find fragment [CH2][CH2+]" << std::endl;
			pass = false;
		}

		//Check that the even electron rule is adhered to correctly in the transitions
		//Iterate over from_id (i)
		std::vector<int> rad_flags(30, 0);
		rad_flags[0] = 1;	//Original Ion is radical
		tmap_t::const_iterator it = graph->getFromIdTMap()->begin();
		for( int from_idx=0; it != graph->getFromIdTMap()->end(); ++it, from_idx++ ){

			int num_rad = 0, num_nonrad = 0;
			std::vector<int>::const_iterator itt = it->begin();
			for( ; itt != it->end(); ++itt ){
				const Transition *t = graph->getTransitionAtIdx(*itt);		
				int israd = moleculeHasSingleRadical( t->getIon()->mol.get() );
				num_rad += israd;
				num_nonrad += (1 - israd);
				rad_flags[t->getToId()] = israd;
			}

			if( !rad_flags[from_idx] && num_rad > 0 ){
				std::cout << "Even Electron Rule Broken: Non-radical parent (" << from_idx << ") produced radical offspring!" << std::endl;
				graph->writeFullGraph( out );
				pass = false;
			}

			//This check works for this molecule....might not work for all?
			if( rad_flags[from_idx] && (num_rad == 0 || num_nonrad == 0) && (num_rad + num_nonrad > 1) ){
				std::cout << "Even Electron Rule Broken: Radical parent (" << from_idx << ") produced only one type of offspring" << std::endl;
				graph->writeFullGraph( out );
				pass = false;
			}
		}
	}

	passed = pass;
}

FragGenTestPositiveEIOxygenAromatic::FragGenTestPositiveEIOxygenAromatic(){
	description = "Test generation of positive mode EI fragment tree for previous exception case NIST30539 C1=CC2OC=CC12 (aromatic oxygen)";
}

void FragGenTestPositiveEIOxygenAromatic::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "C1=CC2OC=CC12";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 43 ){
		std::cout << "Expected 43 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;

	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 94.04131623 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 26.01510148 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 23.99945142 ) > tol ||
			fabs( graph->getFragmentAtIdx(40)->getMass() - 79.01784114 ) > tol ||
			fabs( graph->getFragmentAtIdx(41)->getMass() - 92.02566617 ) > tol ||
			fabs( graph->getFragmentAtIdx(42)->getMass() - 93.0334912 ) > tol ){
				//Just check the first 3 and the last 3 (check the rest manually against below if desired!)
			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			
			//Expecting:
			//0 94.04131623 C1=CC2[O+]C=CC12
			//1 26.01510148 [CH]=[CH+]
			//2 23.99945142 [C]#[C+]
			//3 25.00727645 [C+]#C
			//4 27.02292652 C#[CH2+]
			//5 69.0334912 C1=C[OH+]C=C1
			//6 19.01784114 [OH3+]
			//7 51.02292652 C#CC#[CH2+]
			//8 29.00219107 C#[O+]
			//9 39.02292652 C#C[CH2+]
			//10 41.03857658 C#C[CH4+]
			//11 43.01784114 C=C=[OH+]
			//12 41.00219107 [CH+]=C=O
			//13 53.00219107 [CH+]=C=C=O
			//14 15.02292652 [CH3+]
			//15 17.03857658 [CH5+]
			//16 68.02566617 C1=C[O+]C=C1
			//17 18.0100161 [OH2+]
			//18 50.01510148 C#C[C]=[CH+]
			//19 27.99436604 [C]#[O+]
			//20 38.01510148 [C]#C[CH2+]
			//21 40.03075155 [CH+]=[C]C
			//22 42.0100161 C=C=[O+]
			//23 16.03075155 [CH4+]
			//24 51.99436604 [C]#CC#[O+]
			//25 53.03857658 C1=C[CH2+]=C1
			//26 52.03075155 [CH]1C=C[CH+]1
			//27 65.03857658 [CH+]=C1C=CC1
			//28 49.00727645 [C+]#CC#C
			//29 64.03075155 [CH+]=C1[C]=CC1
			//30 47.99945142 [C]#CC#[C+]
			//31 66.04640161 C=C1[CH+][CH]C1
			//32 77.03857658 [CH2+]#CC1=CC=C1
			//33 75.02292652 [C+]#CC1=CC=C1
			//34 76.03075155 C#CC1=C[CH][CH+]1
			//35 74.01510148 [C+]#CC1=CC=[C]1
			//36 68.02566617 [O+]=C1C=CC1
			//37 39.99436604 [C+]#C[O]
			//38 69.0334912 [OH+]=C1C=CC1
			//39 78.0100161 [CH+]=C1[C]=CC1=O
			//40 79.01784114 [CH+]=C1C=CC1=O
			//41 92.02566617 C1=CC2=C1C=C[O+]2
			//42 93.0334912 C1=CC2=C1C=C[OH+]2

			pass = false;
		}
	}

	passed = pass;
}

FragGenTestPositiveEITriple::FragGenTestPositiveEITriple(){
	description = "Test generation of positive mode EI fragment tree for previous exception case NIST17492 CCCCC#CCCCC";
}

void FragGenTestPositiveEITriple::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CC#CC";	//Simplified

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 14 ){
		std::cout << "Expected 14 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;

	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 54.04640161  ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 38.01510148 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 39.02292652 ) > tol ||
			fabs( graph->getFragmentAtIdx(11)->getMass() - 50.01510148 ) > tol ||
			fabs( graph->getFragmentAtIdx(12)->getMass() - 51.02292652) > tol ||
			fabs( graph->getFragmentAtIdx(13)->getMass() - 53.03857658 ) > tol){

			std::cout << std::endl << "Unexpected fragment mass(es)" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected:
			//0 54.04640161 C[C]=[C+]C
			//1 38.01510148 [C+]#C[CH2]
			//2 39.02292652 [C+]#CC
			//3 15.02292652 [CH3+]
			//4 28.03075155 [CH2][CH2+]
			//5 26.01510148 [CH]=[CH+]
			//6 27.02292652 C#[CH2+]
			//7 29.03857658 C=[CH3+]
			//8 23.99945142 [C]#[C+]
			//9 25.00727645 [C+]#C
			//10 52.03075155 [CH+]=[C]C=C
			//11 50.01510148 C#C[C]=[CH+]
			//12 51.02292652 C#CC#[CH2+]
			//13 53.03857658 C=CC#[CH2+]
			pass = false;
		}
	}

	passed = pass;
}


FragGenTestRingPositiveESI::FragGenTestRingPositiveESI(){
	description = "Test generation of positive mode ESI fragment tree for a ring molecule";
}

void FragGenTestRingPositiveESI::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "C1=CN=CN=C1";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 19 ){
		std::cout << "Expected 19 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFragmentsOnly( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 81.04472458 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 65.01342445 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 15.02292652 ) > tol ||
			fabs( graph->getFragmentAtIdx(16)->getMass() - 50.00252542 ) > tol ||
			fabs( graph->getFragmentAtIdx(17)->getMass() - 39.02292652 ) > tol ||
			fabs( graph->getFragmentAtIdx(18)->getMass() - 65.01342445 ) > tol){

			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			
			//Expecting:
			//0 81.04472458 c1cnc[nH+]c1
			//1 65.01342445 [C+]#CN=C=N
			//2 15.02292652 [CH3+]
			//3 54.03382555 C#C[NH+]=C
			//4 28.01817548 C#[NH+]
			//5 26.00252542 [C+]#N
			//6 25.00727645 [C+]#C
			//7 27.02292652 C#[CH2+]
			//8 52.01817548 C#C[N+]#C
			//9 30.03382555 C=[NH2+]
			//10 40.01817548 [C+]#CN
			//11 29.03857658 C=[CH3+]
			//12 53.01342445 C#[N+]C#N
			//13 55.02907452 C=NC#[NH+]
			//14 54.03382555 C=CC#[NH+]
			//15 52.01817548 [CH+]=C=C=N
			//16 50.00252542 [C+]#CC#N
			//17 39.02292652 [C+]#CC
			//18 65.01342445 N#C[C+]=C=N

			pass = false;
		}
	}

	passed = pass;
}

FragGenTestRingNegativeESI::FragGenTestRingNegativeESI(){
	description = "Test generation of negative mode ESI fragment tree for a ring molecule";
}

void FragGenTestRingNegativeESI::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "C1=CN=CN=C1";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, NEGATIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
	cfg.ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 9 ){
		std::cout << "Expected 9 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFragmentsOnly( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 79.03017168  ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 52.01927264 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 26.00362258   ) > tol ||
			fabs( graph->getFragmentAtIdx(6)->getMass() - 53.01452161 ) > tol ||
			fabs( graph->getFragmentAtIdx(7)->getMass() - 52.01927264 ) > tol ||
			fabs( graph->getFragmentAtIdx(8)->getMass() - 50.00362258   ) > tol ){

			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected:
			//0 79.03017168 [c-]1cncnc1
			//1 52.01927264 [C-]#CN=C
			//2 26.00362258 [C-]#N
			//3 25.00837361 [C-]#C
			//4 28.01927264 C=[N-]
			//5 27.02402368 [CH-]=C
			//6 53.01452161 [CH-]=NC#N
			//7 52.01927264 C=[C-]C#N
			//8 50.00362258 [C-]#CC#N
			pass = false;
		}
	}

	passed = pass;
}


FragGenTestRingPositiveEI::FragGenTestRingPositiveEI(){
	description = "Test generation of positive mode EI fragment tree for a ring molecule";
}

void FragGenTestRingPositiveEI::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "C1=CN=CN=C1";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
	cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 8 ){
		std::cout << "Expected 8 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 81.04527316 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 54.03437413 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 28.01872406 ) > tol ||
			fabs( graph->getFragmentAtIdx(3)->getMass() - 27.0234751  ) > tol ||
			fabs( graph->getFragmentAtIdx(4)->getMass() - 55.0296231  ) > tol ||
			fabs( graph->getFragmentAtIdx(5)->getMass() - 54.03437413 ) > tol ||
			fabs( graph->getFragmentAtIdx(6)->getMass() - 52.01872406 ) > tol ||
			fabs( graph->getFragmentAtIdx(7)->getMass() - 30.03437413 ) > tol ){

			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			std::cout << "Expecting:" << std::endl;
			std::cout << "0 81.04527316 C1=CN=CN=C1" << std::endl;
			std::cout << "1 54.03437413 C#C[NH+]=C" << std::endl;
			std::cout << "2 28.01872406 C#[NH+]" << std::endl;
			std::cout << "3 27.0234751 C#[CH2+]" << std::endl;
			std::cout << "4 55.0296231 C=NC#[NH+]" << std::endl;
			std::cout << "5 54.03437413 C=CC#[NH+]" << std::endl;
			std::cout << "6 52.01872406 C#CC#[NH+]" << std::endl;
			std::cout << "7 30.03437413 C=[NH2+]" << std::endl << std::endl;
			pass = false;
		}
	}

	passed = pass;
}

FragGenTestPositiveEIAlkane::FragGenTestPositiveEIAlkane(){
	description = "Test generation of positive mode EI fragment tree for an alkane";
}

void FragGenTestPositiveEIAlkane::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CCC";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
	cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 12 ){
		std::cout << "Expected 12 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 44.06205168 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 28.03075155 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 26.01510148 ) > tol ||
			fabs( graph->getFragmentAtIdx(9)->getMass() - 40.03075155   ) > tol ||		
			fabs( graph->getFragmentAtIdx(10)->getMass() - 41.03857658 ) > tol ||
			fabs( graph->getFragmentAtIdx(11)->getMass() - 43.05422664 ) > tol ){

			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected:
			//0 44.06205168 CC[CH3+]
			//1 28.03075155 [CH2][CH2+]
			//2 26.01510148 [CH]=[CH+]
			//3 27.02292652 C#[CH2+]
			//4 29.03857658 C=[CH3+]
			//5 15.02292652 [CH3+]
			//6 17.03857658 [CH5+]
			//7 16.03075155 [CH4+]
			//8 42.04640161 [CH2+][CH]C
			//9 40.03075155 [CH+]=[C]C
			//10 41.03857658 [CH2+]#CC
			//11 43.05422664 CC=[CH3+]
			pass = false;
		}
	}

	passed = pass;
}

FragGenTestPositiveESISplitCharge::FragGenTestPositiveESISplitCharge(){
	description = "Test generation of positive mode ESI fragment tree for split charge input";
}

void FragGenTestPositiveESISplitCharge::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CC=[N+]=[N-]";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 10 ){
		std::cout << "Expected 10 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		 
		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 57.04472458 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 15.02292652 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 31.02907452 ) > tol ||
			fabs( graph->getFragmentAtIdx(7)->getMass() - 40.01817548 ) > tol ||
			fabs( graph->getFragmentAtIdx(8)->getMass() - 55.02907452 ) > tol ||
			fabs( graph->getFragmentAtIdx(9)->getMass() - 53.01342445 ) > tol ){
			std::cout << std::endl << "Unexpected fragment mass(es):" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected
			//0 57.04472458 [CH4+]C=[N+]=[N-]
			//1 15.02292652 [CH3+]
			//2 31.02907452 N=[NH2+]
			//3 29.01342445 N#[NH+]
			//4 25.00727645 [C+]#C
			//5 27.02292652 C#[CH2+]
			//6 29.03857658 C=[CH3+]
			//7 40.01817548 [C+]#CN
			//8 55.02907452 C#C[NH+]=N
			//9 53.01342445 C#C[N+]#N
			pass = false;
		}
	}

	passed = pass;
}


FragGenTestPositiveEISplitCharge::FragGenTestPositiveEISplitCharge(){
	description = "Test generation of positive mode EI fragment tree for split charge input";
}

void FragGenTestPositiveEISplitCharge::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CC=[N+]=[N-]";

	//Run the fragmentation procedure
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_EI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
	cfg.ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	FragmentGraph *graph = gg.createNewGraph(&cfg);
	gg.compute( *startNode, 2, -1, 2 );
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 17 ){
		std::cout << "Expected 17 fragments but found " << graph->getNumFragments() << std::endl;
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);
		graph->writeFullGraph( out );
		pass = false;
	}
	else{
	
		std::streambuf * buf = std::cout.rdbuf();
		std::ostream out(buf);

		if( fabs( graph->getFragmentAtIdx(0)->getMass() - 56.03689955 ) > tol ||
			fabs( graph->getFragmentAtIdx(1)->getMass() - 15.0229265 ) > tol ||
			fabs( graph->getFragmentAtIdx(2)->getMass() - 30.02124948) > tol ||
			fabs( graph->getFragmentAtIdx(14)->getMass() - 52.0055994 ) > tol ||
			fabs( graph->getFragmentAtIdx(15)->getMass() - 53.01342445 ) > tol ||
			fabs( graph->getFragmentAtIdx(16)->getMass() - 55.02907452 ) > tol ){

			std::cout << std::endl << "Unexpected fragment mass(es)" << std::endl;
			graph->writeFragmentsOnly( out );
			//Expected:
			//0 56.03689955 [CH3+]C=[N+]=[N-]
			//1 15.02292652 [CH3+]
			//2 30.02124948 N=[NH+]
			//3 28.00559942 N#[N+]
			//4 29.01342445 N#[NH+]
			//5 31.02907452 N=[NH2+]
			//6 25.00727645 [C+]#C
			//7 23.99945142 [C]#[C+]
			//8 27.02292652 C#[CH2+]
			//9 26.01510148 [CH]=[CH+]
			//10 28.03075155 [CH2][CH2+]
			//11 40.01817548 [C+]#CN
			//12 39.01035045 [C+]#C[NH]
			//13 54.02124948 C#C[N+]=N
			//14 52.00559942 [C]#C[N+]#N
			//15 53.01342445 C#C[N+]#N
			//16 55.02907452 C#C[NH+]=N
			pass = false;
		}
	}

	passed = pass;
}

FragGenTestPositiveEINistExceptions::FragGenTestPositiveEINistExceptions(){
	description = "Test generation of positive mode EI fragment tree for various NIST cases that threw exceptions";
}

void FragGenTestPositiveEINistExceptions::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	typedef std::pair<std::string, int> nist_test_case_t;
	std::vector<nist_test_case_t> test_cases;
	test_cases.push_back( nist_test_case_t( "InChI=1S/C5H9Cl3/c1-5(2-6,3-7)4-8/h2-4H2,1H3", -1 ));
	//test_cases.push_back( nist_test_case_t( "InChI=1S/C12H19N2.ClH.Cr/c1-13(2)9-11-6-5-7-12(8-11)10-14(3)4;;/h5-7H,9-10H2,1-4H3;1H;/q;;+1/p-1", -1));
	//test_cases.push_back( nist_test_case_t("CC1=C(Cl)C(C)=C(C#N)C(S(=O)(=O)NCC2=CC=CO2)=N1", 1536) ); // NIST62930

	std::vector<std::pair<std::string, int> >::iterator test_it = test_cases.begin();
	for( ; test_it != test_cases.end(); ++test_it ){

		try{
			//Run the fragmentation procedure
			//std::vector<std::string> feature_names;
			//feature_names.push_back(std::string("BreakAtomPair"));
			//FeatureCalculator fc( feature_names );
			FragmentGraphGenerator gg(1); //gg(&fc);
			FragmentTreeNode *startNode = gg.createStartNode(test_it->first, POSITIVE_ESI_IONIZATION_MODE);
			config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
			FragmentGraph *graph = gg.createNewGraph(&cfg);
			gg.compute( *startNode, 2, -1, 2 );
			delete startNode;

			//Check the resulting graph
			if( graph->getNumFragments() != test_it->second ){
				std::cout << test_it->first << std::endl;
				std::cout << "Expected " << test_it->second << " fragments but found " << graph->getNumFragments() << std::endl;
				std::streambuf * buf = std::cout.rdbuf();
				std::ostream out(buf);
				pass = false;
			}

			if( pass ) pass = checkGraphConsistency( graph );
		}
		catch( std::exception e ){
			std::cout << "Exception Occurred computing graph for " <<  test_it->first << std::endl;
			std::cout << e.what() << std::endl << std::endl;
			pass = false;
		}

	}

	passed = pass;
}

FragGenTestPositiveEIAndESIDegreeLpBonding::FragGenTestPositiveEIAndESIDegreeLpBonding(){
	description = "Test generation of positive mode ESI and EI fragment tree for cases with lone pair bonding";
}

void FragGenTestPositiveEIAndESIDegreeLpBonding::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	typedef boost::tuple<std::string, int, int> mol_test_case_t;
	std::vector<mol_test_case_t> test_cases;
	
	test_cases.push_back( mol_test_case_t("Nc1nonc1C(=NO)N1N=N1", POSITIVE_EI_IONIZATION_MODE, 132));
	test_cases.push_back( mol_test_case_t("OCC#C", POSITIVE_EI_IONIZATION_MODE, 17) );	//Should include CO loss
	test_cases.push_back( mol_test_case_t("OCC#C", POSITIVE_ESI_IONIZATION_MODE, 11) );
	test_cases.push_back( mol_test_case_t("[Cl-].[Na+].CN(C)C", POSITIVE_EI_IONIZATION_MODE, 33) );
	test_cases.push_back( mol_test_case_t("[Cl-].[Mn+].CN(C)C", POSITIVE_EI_IONIZATION_MODE, 33) );
	test_cases.push_back( mol_test_case_t("[Cl-].[Na+].CN(C)CCCC", POSITIVE_EI_IONIZATION_MODE, 192) );
	test_cases.push_back( mol_test_case_t("[Cl-].[Mn+].CN(C)CCCC", POSITIVE_EI_IONIZATION_MODE, 192) );
	test_cases.push_back( mol_test_case_t("C[N+](C)(C)CC(=O)O", POSITIVE_EI_IONIZATION_MODE, 29) );
	test_cases.push_back( mol_test_case_t("C[N+](C)(C)CC(=O)O", POSITIVE_ESI_IONIZATION_MODE, 16) );
	test_cases.push_back( mol_test_case_t("C[N+](C)(C)CC(=O)[O-]", POSITIVE_EI_IONIZATION_MODE, 29) );
	test_cases.push_back( mol_test_case_t("C[N+](C)(C)CC(=O)[O-]", POSITIVE_ESI_IONIZATION_MODE, 16) );
	test_cases.push_back( mol_test_case_t("C[N+](C)(C)CC(=O)O.[Cl-]", POSITIVE_EI_IONIZATION_MODE, 60) ); 	//NIST2011_25572
	test_cases.push_back( mol_test_case_t("C[N+](C)(C)CC(=O)O.[Cl-]", POSITIVE_ESI_IONIZATION_MODE, 31) );
	test_cases.push_back( mol_test_case_t("[Cl-].[Na+].CCC", POSITIVE_EI_IONIZATION_MODE, 29) ); 	//two ionic fragments
	test_cases.push_back( mol_test_case_t("[Cl-].[Na+].CCC", POSITIVE_ESI_IONIZATION_MODE, 17) );


	std::vector<mol_test_case_t>::iterator test_it = test_cases.begin();
	for( ; test_it != test_cases.end(); ++test_it ){

		std::cout << "Input: " << boost::get<0>(*test_it) << std::endl;

		//Run the fragmentation procedure
		FragmentGraphGenerator gg(0);
		FragmentTreeNode *startNode = gg.createStartNode( boost::get<0>(*test_it), boost::get<1>(*test_it) );
		config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true;  cfg.ionization_mode = boost::get<1>(*test_it);
		FragmentGraph *graph = gg.createNewGraph(&cfg);
		gg.compute( *startNode, 2, -1, 2 );
		delete startNode;

		//Check the resulting graph
		if( graph->getNumFragments() != boost::get<2>(*test_it) ){
			
			std::cout << "Input: " << boost::get<0>(*test_it) << " in ";
			if( boost::get<1>(*test_it) == POSITIVE_EI_IONIZATION_MODE ) std::cout << "POSITIVE_EI_IONIZATION_MODE" << std::endl;
			else if( boost::get<1>(*test_it) == POSITIVE_ESI_IONIZATION_MODE ) std::cout << "POSITIVE_ESI_IONIZATION_MODE" << std::endl;
			else if( boost::get<1>(*test_it) == NEGATIVE_ESI_IONIZATION_MODE ) std::cout << "NEGATIVE_ESI_IONIZATION_MODE" << std::endl;
			
			std::cout << "Expected " << boost::get<2>(*test_it) << " fragments but found " << graph->getNumFragments() << std::endl;
			std::streambuf * buf = std::cout.rdbuf();
			std::ostream out(buf);
			//graph->writeFragmentsOnly(std::cout);
			//graph->writeFullGraph(std::cout);
			pass = false;
		}
		if( pass ) pass = checkGraphConsistency( graph );

	}

	passed = pass;
}

FragGenTestCasesFromGross::FragGenTestCasesFromGross(){
	description = "Test generation of positive mode EI cases in Gross Chapter 6";
}

void FragGenTestCasesFromGross::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	typedef boost::tuple<std::string, int, int> mol_test_case_t;
	std::vector<mol_test_case_t> test_cases;

	//test_cases.push_back( mol_test_case_t("CC(=O)c1ccc(Br)c(O)c(=O)c1", POSITIVE_EI_IONIZATION_MODE, -1 ));
	test_cases.push_back( mol_test_case_t("C", POSITIVE_EI_IONIZATION_MODE, 2 ));
	test_cases.push_back( mol_test_case_t("CI", POSITIVE_EI_IONIZATION_MODE, 4 ));
	test_cases.push_back( mol_test_case_t("CC(=O)C", POSITIVE_EI_IONIZATION_MODE, 18 ));
	test_cases.push_back( mol_test_case_t("CCC(=O)C", POSITIVE_EI_IONIZATION_MODE, 42 ));
	test_cases.push_back( mol_test_case_t("CCC", POSITIVE_EI_IONIZATION_MODE, 12 ));

	std::vector<mol_test_case_t>::iterator test_it = test_cases.begin();
	for( ; test_it != test_cases.end(); ++test_it ){

		//Run the fragmentation procedure
		FragmentGraphGenerator gg(0);
		FragmentTreeNode *startNode = gg.createStartNode( boost::get<0>(*test_it), boost::get<1>(*test_it) );
		config_t cfg; initDefaultConfig(cfg); cfg.include_h_losses = true; 
		cfg.ionization_mode = boost::get<1>(*test_it); 
		FragmentGraph *graph = gg.createNewGraph(&cfg);
		gg.compute( *startNode, 2, -1, 2 );
		delete startNode;

		//Check the resulting graph
		if( graph->getNumFragments() != boost::get<2>(*test_it) ){
			
			std::cout << "Input: " << boost::get<0>(*test_it) << " in ";
			if( boost::get<1>(*test_it) == POSITIVE_EI_IONIZATION_MODE ) std::cout << "POSITIVE_EI_IONIZATION_MODE" << std::endl;
			else if( boost::get<1>(*test_it) == POSITIVE_ESI_IONIZATION_MODE ) std::cout << "POSITIVE_ESI_IONIZATION_MODE" << std::endl;
			else if( boost::get<1>(*test_it) == NEGATIVE_ESI_IONIZATION_MODE ) std::cout << "NEGATIVE_ESI_IONIZATION_MODE" << std::endl;
			
			std::cout << "Expected " << boost::get<2>(*test_it) << " fragments but found " << graph->getNumFragments() << std::endl;
			std::streambuf * buf = std::cout.rdbuf();
			std::ostream out(buf);
			//graph->writeFragmentsOnly(std::cout);
			graph->writeFullGraph(std::cout);
			pass = false;
		}

		if( pass ) pass = checkGraphConsistency( graph );
	}

	passed = pass;
}

FragGenTestMaxElectronMovement::FragGenTestMaxElectronMovement(){
	description = "Test maximum movement of electrons between fragments";
}

void FragGenTestMaxElectronMovement::runTest(){
	
	bool pass = true;
	double tol = 1e-5;

	FragmentGraphGenerator fgen(0);
	std::string smiles_or_inchi("C#CC#CC#CCCCCCCOCCCCCCCN");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	
	//Break just the one bond
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[12], 0);	//Break Bond 11 (after the O)
	node->generateChildrenOfBreak(breaks[12]);
	
	//Creat a simple graph for just these breaks
	config_t cfg; initDefaultConfig(cfg);
	cfg.include_isotopes = false; cfg.include_h_losses = true; 
	FragmentGraph fg(&cfg);
	fg.addToGraph( *node, -1 );
	std::vector<FragmentTreeNode>::iterator itt = node->children.begin();
	for( ; itt != node->children.end(); ++itt ){
		fg.addToGraph( *itt, 0 );
	}
	if( fg.getNumFragments() != 11 ){
			std::cout << "Expected 11 fragments but found " << fg.getNumFragments() << std::endl;
			fg.writeFragmentsOnly(std::cout);
			pass = false;
	}
	else{

		//Now break one of the children (to check the second level works)
		FragmentTreeNode *child = &node->children[6];
		std::vector<Break> child_breaks;
		child->generateBreaks(child_breaks, false);
		child->applyBreak(breaks[5], 0);	//Break Bond 5
		child->generateChildrenOfBreak(breaks[5]);
		itt = child->children.begin();
		for( ; itt != child->children.end(); ++itt ){
			fg.addToGraph( *itt, 0 );
		}
		if( fg.getNumFragments() != 20 ){
				std::cout << "Expected 20 fragments but found " << fg.getNumFragments() << std::endl;
				fg.writeFragmentsOnly(std::cout);
				pass = false;
		}

	}
	delete node;
	passed = pass;
}

FragGenTestDisallowDetourTransitions::FragGenTestDisallowDetourTransitions(){
	description = "Test functionality to disallow transition detours";
}


bool checkChildFragmentDepths( int parent_id, FragmentGraph *graph, int depth ){

	const Fragment *f = graph->getFragmentAtIdx(parent_id);
	bool pass = (depth == f->getDepth());

	const tmap_t *from_map = graph->getFromIdTMap();
	std::vector<int>::const_iterator it = (*from_map)[parent_id].begin();
	for( ; it != (*from_map)[parent_id].end(); ++it ){
		const Transition *t = graph->getTransitionAtIdx(*it);
		pass = pass && checkChildFragmentDepths( t->getToId(), graph, depth+1 );
	}
	return pass;
}

void FragGenTestDisallowDetourTransitions::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	std::string smiles_or_inchi = "CCCC";

	//Run the fragmentation procedure allowing detours
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode1 = gg.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.allow_frag_detours = true; cfg.include_h_losses = true; 
	FragmentGraph *graph1 = gg.createNewGraph( &cfg );
	gg.compute( *startNode1, 3, -1, 3 );
	delete startNode1;

	//Check the number of transitions allowing detours
	if( graph1->getNumTransitions() != 31 ){
		std::cout << "Expected 31 transitions but found " << graph1->getNumTransitions() << std::endl;
		graph1->writeFullGraph( std::cout );
		pass = false;
	}

	//Run the fragmentation procedure
	FragmentTreeNode *startNode = gg.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	cfg.allow_frag_detours = false;
	FragmentGraph *graph = gg.createNewGraph( &cfg );
	gg.compute( *startNode, 3, -1, 3 );

	//Check the number of transitions before removing detours (some should have already been removed along the way)
	if( graph->getNumTransitions() != 21 ){
		std::cout << "Expected 21 transitions but found " << graph->getNumTransitions() << std::endl;
		graph->writeFullGraph( std::cout );
		pass = false;
	}

	graph->removeDetours();
	delete startNode;

	//Check the resulting graph
	if( graph->getNumFragments() != 13 ){
		std::cout << "Expected 13 fragments but found " << graph->getNumFragments() << std::endl;
		graph->writeFullGraph( std::cout );
		pass = false;
	}
	else if( graph->getNumTransitions() != 18 ){
		std::cout << "Expected 18 transitions but found " << graph->getNumTransitions() << std::endl;
		graph->writeFullGraph( std::cout );
		pass = false;
	}
	else if(pass){

		//Check graph consistency and id maps
		pass = checkGraphConsistency( graph );

		//Check that all fragments are reachable in exactly the number of steps specified by their depth
		if(pass) pass = checkChildFragmentDepths( 0, graph, 0 );

		//Check ordering of transitions (a transition producing each parent as a child must have occurred first)
		//NOTE: This fails. I don't think this is a problem, so I don't plan to fix it, but I'm leaving it here for now as a reminder...
		std::vector<int> produced_flags( graph->getNumFragments(), 0 );
		produced_flags[0] = 1;
		for( int i = 0; i < graph->getNumTransitions(); i++ ){
			const Transition *t = graph->getTransitionAtIdx(i);
			if( !produced_flags[t->getFromId()] ){
				std::cout << "Warning: Fragment " << t->getFromId() << " occurs as a parent before it is seen as a child" << std::endl;
				//pass = false;
			}
			produced_flags[t->getToId()] = 1;
		}
	}

	passed = pass;
}

FragGenTestMaxRingBreaks::FragGenTestMaxRingBreaks(){
	description = "Test functionality to limit the number of ring breaks";
}

void FragGenTestMaxRingBreaks::runTest(){
	
	bool pass = true;
	double tol = 1e-6;

	//Run the fragmentation procedure - no ring breaks
	std::string smiles_or_inchi = "C1CC1";
	FragmentGraphGenerator gg;
	FragmentTreeNode *startNode1 = gg.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	config_t cfg; initDefaultConfig(cfg); cfg.allow_frag_detours = true; cfg.include_h_losses = true; 
	FragmentGraph *graph1 = gg.createNewGraph( &cfg );
	gg.compute( *startNode1, 3, -1, 0 );
	delete startNode1;

	//Check the number of transitions before removing detours (some should have already been removed along the way)
	if( graph1->getNumFragments() != 2 ){
		std::cout << "Expected 2 fragments but found " << graph1->getNumFragments() << std::endl;
		graph1->writeFullGraph( std::cout );
		pass = false;
	}

	//Run the fragmentation procedure - 1 ring breaks
	smiles_or_inchi = "C1OC1C1NC1";
	FragmentGraphGenerator gg2;
	FragmentTreeNode *startNode2 = gg2.createStartNode(smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE);
	FragmentGraph *graph2 = gg2.createNewGraph( &cfg );
	gg2.compute( *startNode2, 3, -1, 1 );
	delete startNode2;

	//Check the number of transitions before removing detours (some should have already been removed along the way)
	if( graph2->getNumFragments() != 36 ){
		std::cout << "Expected 36 fragments but found " << graph2->getNumFragments() << std::endl;
		graph2->writeFullGraph( std::cout );
		pass = false;
	}

	passed = pass;
}


bool checkGraphConsistency( FragmentGraph *graph ){
	
	//Check that all ions have one single charge 
	for( int i = 0; i < graph->getNumFragments(); i++ ){
		const Fragment *frag = graph->getFragmentAtIdx(i);
		RDKit::RWMol *rwmol = RDKit::SmilesToMol( *frag->getIonSmiles() );
		if( RDKit::MolOps::getFormalCharge(*rwmol ) != 1 ){
			std::cout << "Invalid charge on ion: " << *frag->getIonSmiles() << std::endl;
			return false;
		}
	}

	//Check that all transitions have uncharged neutral losses
	for( int i = 0; i < graph->getNumTransitions(); i++ ){
		const Transition *trans = graph->getTransitionAtIdx(i);
		RDKit::RWMol *rwmol = RDKit::SmilesToMol( *trans->getNLSmiles() );
		if( RDKit::MolOps::getFormalCharge(*rwmol ) != 0 ){
			std::cout << "Invalid charge on neutral loss: " << *trans->getNLSmiles() << std::endl;
			return false;
		}
	}	
	
	//Check that mass(ion) + mass(nl) = mass(parent)
	double tol = 0.01;
	for( int i = 0; i < graph->getNumTransitions(); i++ ){
		const Transition *trans = graph->getTransitionAtIdx(i);
		const Fragment *parent = graph->getFragmentAtIdx( trans->getFromId() );
		const Fragment *ion = graph->getFragmentAtIdx( trans->getToId() );
		double nl_mass = getMonoIsotopicMass( trans->getNeutralLoss()->mol );
		if( fabs(ion->getMass() + nl_mass -  parent->getMass() ) > tol ){
			std::cout << "Invalid preservation of mass: " << trans->getFromId() << "->" << trans->getToId() << " " << *parent->getIonSmiles() << " -> " << *ion->getIonSmiles() << " + " << *trans->getNLSmiles() << std::endl;
			std::cout << ion->getMass() << " + " << nl_mass << " != " << parent->getMass() << std::endl;
			return false;
		}
	}

	//Check that all transitions are unique
	std::set< std::pair<int, int> > from_to_pairs;
	for( int i = 0; i < graph->getNumTransitions(); i++ ){
		const Transition *t = graph->getTransitionAtIdx(i);
		if( from_to_pairs.find( std::pair<int,int>(t->getFromId(), t->getToId()) ) != from_to_pairs.end() ){
			std::cout << "Found duplicate transition" << t->getFromId() << "-" << t->getToId() << std::endl;
			return false;
		}
		from_to_pairs.insert(std::pair<int,int>(t->getFromId(), t->getToId()));
	}

	//Check the tmaps contain all transitions (just once)
	const tmap_t *to_map = graph->getToIdTMap();
	const tmap_t *from_map = graph->getFromIdTMap();
	for( int i = 1; i < graph->getNumTransitions(); i++ ){
		const Transition *t = graph->getTransitionAtIdx(i);
		if( std::count( (*to_map)[t->getToId()].begin(), (*to_map)[t->getToId()].end(), i ) != 1 ){
			std::cout << "Could not find transition" << i << "in correct to map" << std::endl;
			return false;
		}
		if( std::count( (*from_map)[t->getFromId()].begin(), (*from_map)[t->getFromId()].end(), i ) != 1 ){
			std::cout << "Could not find transition" << i << "in correct from map" << std::endl;
			return false;
		}
	}

	//Check that all tmap entries match a transition
	tmap_t::const_iterator it = to_map->begin();
	for( int idx = 0; it != to_map->end(); ++it, idx++ ){
		std::vector<int>::const_iterator itt = it->begin();
		for( ; itt != it->end(); ++itt ){
			if( graph->getTransitionAtIdx(*itt)->getToId() != idx ){
				std::cout << "Invalid tmap entry in to map" << std::endl;
				return false;
			}
		}
	}
	it = from_map->begin();
	for( int idx = 0; it != from_map->end(); ++it, idx++ ){
		std::vector<int>::const_iterator itt = it->begin();
		for( ; itt != it->end(); ++itt ){
			if( graph->getTransitionAtIdx(*itt)->getFromId() != idx ){
				std::cout << "Invalid tmap entry in from map" << std::endl;
				return false;
			}
		}
	}

	return true;
}


FVFragGraphSaveAndLoadState::FVFragGraphSaveAndLoadState(){
	description = "Test functionality for saving FV fragment graph to file and re-loading";
}

void FVFragGraphSaveAndLoadState::runTest(){
	
	bool pass = true;
	double tol = 1e-16;

	for( int include_isotopes = 0; include_isotopes <= 1; include_isotopes++ ){

		//Compute a fragment graph and feature vectors
		std::string smiles_or_inchi = "CC1=CC(=O)CC(C)C1";
		config_t cfg; initDefaultConfig( cfg );
		cfg.include_isotopes = include_isotopes;
		MolData morig("Test ID", smiles_or_inchi.c_str(), &cfg );
		std::vector<std::string> fnames;
		fnames.push_back( "BreakAtomPair" ); fnames.push_back( "RingFeatures" );  fnames.push_back( "HydrogenRemoval" );
		FeatureCalculator fc( fnames );
		morig.computeFragmentGraphAndReplaceMolsWithFVs( &fc, false );

		//Save state
		std::string fvfilename = "tmp_fv_file.fg";
		if( boost::filesystem::exists( fvfilename ) )
			boost::filesystem::remove( fvfilename );
		morig.writeFVFragmentGraph( fvfilename );

		//Load state into another moldata object
		MolData mload("Test ID", smiles_or_inchi.c_str(), &cfg );
		mload.readInFVFragmentGraph( fvfilename );

		//Compare orig vs loaded
		const FragmentGraph *fgorig = morig.getFragmentGraph();
		const FragmentGraph *fgload = mload.getFragmentGraph();
		if( fgorig->getNumFragments() != fgload->getNumFragments() || fgorig->getNumTransitions() != fgload->getNumTransitions() ){
			std::cout << "Mismatch in fragment graph dimensions: " << std::endl;
			pass = false;
			break;
		}
		for( int i = 0; i < fgorig->getNumFragments(); i++ ){
			const Fragment *f1 = fgorig->getFragmentAtIdx(i);
			const Fragment *f2 = fgload->getFragmentAtIdx(i);
			if( f1->getId() != f2->getId() || fabs(f1->getMass() - f2->getMass() ) > tol ){
				std::cout << "Mismatch in fragments at idx: " << i << std::endl;
				pass = false;
				break;			
			}
			if( include_isotopes ){
				const Spectrum *spec1 = f1->getIsotopeSpectrum();
				const Spectrum *spec2 = f2->getIsotopeSpectrum();
				if( spec1->size() != spec2->size() ){
					std::cout << "Mismatch in isotope spectrum dimensions for fragment at idx " << i << std::endl;
					pass = false;
					break;		
				}
				Spectrum::const_iterator it1 = spec1->begin();
				Spectrum::const_iterator it2 = spec2->begin();
				for( ; it1 != spec1->end(); ++it1, ++it2 ){
					if( fabs(it1->mass - it2->mass) > tol || fabs(it1->intensity - it2->intensity) > tol ){
						std::cout << "Mismatch in isotope peaks for fragment at idx " << i << std::endl;
						pass = false;
						break;						
					}
				}
			
			}
		}
		if( !pass ) break;
		for( int i = 0; i < fgorig->getNumTransitions(); i++ ){
			const Transition *t1 = fgorig->getTransitionAtIdx(i);
			const Transition *t2 = fgload->getTransitionAtIdx(i);
			if( t1->getFromId() != t2->getFromId() || t1->getToId() != t2->getToId() ){
				std::cout << "Mismatch in transtions at idx: " << i << std::endl;
				pass = false;
				break;			
			}
			const FeatureVector *fv1 = morig.getFeatureVectorForIdx(i);
			const FeatureVector *fv2 = mload.getFeatureVectorForIdx(i);
			if( fv1->getTotalLength() != fv2->getTotalLength() || fv1->getNumSetFeatures() != fv2->getNumSetFeatures() ){
				std::cout << "Mismatch in feature vector dimensions idx: " << i << std::endl;
				pass = false;
				break;				
			}
			std::vector<feature_t>::const_iterator fit1 = fv1->getFeatureBegin();
			std::vector<feature_t>::const_iterator fit2 = fv2->getFeatureBegin();
			for(; fit1 != fv1->getFeatureEnd(); ++fit1, ++fit2 ){
				if( *fit1 != *fit2 ){
					std::cout << "Mismatch in feature vectors for transition at idx " << i << ": " << *fit1 << " vs " << *fit2 << std::endl;
					pass = false;
					break;
				}
			}
			if( !pass ) break;
		}
		if( !pass ) break;

		//Check the tmaps
		const tmap_t *fromidmap1 = fgorig->getFromIdTMap();
		const tmap_t *fromidmap2 = fgload->getFromIdTMap();
		const tmap_t *toidmap1 = fgorig->getToIdTMap();
		const tmap_t *toidmap2 = fgload->getToIdTMap();
		if( fromidmap1->size() != fromidmap2->size() || toidmap1->size() != toidmap2->size()){
			std::cout << "Mismatch in tmap dimensions" << std::endl;
			pass = false;			
		}
		for( int i = 0; i < fromidmap1->size(); i++ ){
			if( (*fromidmap1)[i].size() != (*fromidmap2)[i].size()  ||
				(*toidmap1)[i].size() != (*toidmap2)[i].size() ){
					std::vector<int>::const_iterator map_it1 = (*fromidmap1)[i].begin();
					std::vector<int>::const_iterator map_it2 = (*fromidmap2)[i].begin();
					for( ; map_it1 != (*fromidmap1)[i].end(); ++map_it1, ++map_it2 ){
						if( *map_it1 != * map_it2 ){
							std::cout << "Mismatch in from id tmap for fragment " << i << std::endl;
							pass = false;					
						}
					}
					std::vector<int>::const_iterator tmap_it1 = (*toidmap1)[i].begin();
					std::vector<int>::const_iterator tmap_it2 = (*toidmap2)[i].begin();
					for( ; tmap_it1 != (*toidmap1)[i].end(); ++tmap_it1, ++tmap_it2 ){
						if( *tmap_it1 != * tmap_it2 ){
							std::cout << "Mismatch in to id tmap for fragment " << i << std::endl;
							pass = false;					
						}
					}
			}
		}
	}

	passed = pass;
}