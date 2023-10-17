/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.cpp
#
# Description: Test code for features.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "features_test.h"
#include "Util.h"
#include "MolData.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <boost/filesystem.hpp>

void initMolProps( romol_ptr_t &mol ){
	RDKit::ROMol::AtomIterator ai; 
	for( ai = mol.get()->beginAtoms(); ai != mol.get()->endAtoms(); ++ai ){	
		(*ai)->setProp("Root",0);
		(*ai)->setProp("OtherRoot",0);
	}
	mol.get()->setProp("IsRingBreak", 0);
}

FeaturesTestInit::FeaturesTestInit(){
	description = "Test initialisation of feature class";
}

void FeaturesTestInit::runTest(){
	
	bool pass = true;
    
	//Valid Config
	std::string config_filename = "tests/test_data/valid_feature_config.txt";
	FeatureCalculator *fc = new FeatureCalculator( config_filename );
	std::vector<std::string> fnames = fc->getFeatureNames();
	if( fnames.size() != 2 ){
		std::cout << "Unexpected number of feature names" << std::endl;
		pass = false;
	}
	else{
		if( fnames[0] != "BreakAtomPair" ||
			fnames[1] != "RingFeatures" ){
			std::cout << "Unexpected feature names" << std::endl;
			pass = false;		
		}
		if( fc->getNumFeatures() != 85 ){
			std::cout << "Unexpected feature count";
			pass = false;			
		}
	}
	delete fc;

	//Invalid Config
	bool except_thrown = false;
	try{
		std::string config_filename = "tests/test_data/invalid_feature_config.txt";
		FeatureCalculator *fc = new FeatureCalculator( config_filename );
	}
	catch( InvalidConfigException e ){
		except_thrown = true;
	}
	if( !except_thrown ){
		std::cout << "Allowed invalid feature configuration" << std::endl;
		pass = false;		
	}
	passed = pass;

}

FeaturesTestBreakAtomPair::FeaturesTestBreakAtomPair(){
	description = "Test BreakAtomPair feature";
}

void FeaturesTestBreakAtomPair::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("BreakAtomPair");
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Non-Ring Break C-N
	RDKit::Atom *null_atom = NULL;
	romol_ptr_t ion = createMolPtr("C");
	initMolProps(ion);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom);
	romol_ptr_t nl =  createMolPtr("N");
	initMolProps(nl);
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom);

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 3 ){
			std::cout << "Unexpected value for non-ring C-N root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Non-Ring Break X-C
	ion = createMolPtr("B");
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), null_atom );
	nl = createMolPtr("C");
	initMolProps(nl);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), null_atom );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 61 ){
			std::cout << "Unexpected value for non-ring X-C root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring Break double C-N
	ion = createMolPtr("CC");
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );
	nl = romol_ptr_t(RDKit::SmilesToMol("NN"));
	initMolProps(nl);
	nl.get()->setProp("IsRingBreak", 1);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 4 ){
			std::cout << "Unexpected value for ring double C-N root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring Break C-N, X-X
	ion = createMolPtr("CB");
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );
	nl = createMolPtr("NB");
	initMolProps(nl);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );
	nl.get()->setProp("IsRingBreak", 1);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 4 ){
			std::cout << "Unexpected value for ring C-N, X-X root pair" << std::endl;
			pass = false;			
		}
		if( fv->getFeature(2) != 72 ){
			std::cout << "Unexpected value for ring C-N, X-X root pair" << std::endl;
			pass = false;			
		}
	}
	delete fv;
	delete fc;
	passed = pass;
}

FeaturesTestRootPairs::FeaturesTestRootPairs(){
	description = "Test IonRootPairs and NLRootPairs feature";
}

void FeaturesTestRootPairs::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonRootPairs"); //NLRootPairs should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Non-ring "C-N,C-N,C-X"
	RDKit::Atom *null_atom = NULL;
	romol_ptr_t ion = createMolPtr("C(N)(B)N");
	romol_ptr_t nl = createMolPtr("C");
	initMolProps(ion);
	initMolProps(nl);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 6 ) ok = false;
		if( fv->getFeature(2) != 7 ) ok = false;
		if( fv->getFeature(3) != 22 ) ok = false;
		if(!ok){
			std::cout << "Unexpected value for non-ring C-N,C-N,C-X" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Non-ring "X-X,X-N"
	ion = createMolPtr("B(B)N");
	initMolProps(ion);
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), null_atom );
	
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 126 ) ok = false;	
		if( fv->getFeature(2) != 142 ) ok = false;	
		if(!ok){
			std::cout << "Unexpected value for non-ring X-X,X-N" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring "C-N,C-N,X-X,C-X,X-N"
	ion = createMolPtr("C(N)(B)NNBB");
	initMolProps(ion);	
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(5) );
	nl.get()->setProp("IsRingBreak", 1);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 6 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 8 ) ok = false;
		if( fv->getFeature(2) != 9 ) ok = false;
		if( fv->getFeature(3) != 24 ) ok = false;
		if( fv->getFeature(4) != 128 ) ok = false;	
		if( fv->getFeature(5) != 144 ) ok = false;	
		if(!ok){
			std::cout << "Unexpected value for ring C-N,C-N,X-X,C-X,X-N" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Empty pairs (set feature indicating no pairs)
	ion = createMolPtr("C");
	initMolProps(ion);	
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), null_atom );
	nl.get()->setProp("IsRingBreak", 0);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 ){ 
			std::cout << "Missing feature flagging no pairs" << std::endl;
			pass = false;
		}
	}
	delete fv;

	delete fc;
	passed = pass;

}

FeaturesTestRootTriples::FeaturesTestRootTriples(){
	description = "Test IonRootTriples and NLRootTriples feature";
}

void FeaturesTestRootTriples::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("NLRootTriples"); //IonRootTriples should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Non-ring "C-C-N,C-C-N,C-X-C"
	romol_ptr_t ion = createMolPtr("C");
	romol_ptr_t nl = createMolPtr("C(BC)C(N)N");
	initMolProps(ion);
	initMolProps(nl);
	RDKit::Atom *null_atom = NULL;
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom ); 

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 6 ) ok = false;
		if( fv->getFeature(2) != 7 ) ok = false;
		if( fv->getFeature(3) != 122 ) ok = false;
		if(!ok){
			std::cout << "Unexpected value for non-ring C-C-N,C-C-N,C-X-C" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Non-ring "N-X-X"
	nl = createMolPtr("NBB");
	initMolProps(nl);	
	rtd_nl = RootedROMolPtr(nl, nl.get()->getAtomWithIdx(0), null_atom );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 286 ){
			std::cout << "Unexpected value for non-ring N-X-X" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Ring "C-C-N,C-C-N,C-X-C,N-X-X"
	ion = createMolPtr("C");
	nl = createMolPtr("C(BC)C(N)NBBN");
	initMolProps(ion);
	initMolProps(nl);	
	rtd_nl = RootedROMolPtr(nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(8) );
	nl.get()->setProp("IsRingBreak", 1);

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 5 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		bool ok = true;
		if( fv->getFeature(1) != 8 ) ok = false;
		if( fv->getFeature(2) != 9 ) ok = false;
		if( fv->getFeature(3) != 124 ) ok = false;
		if( fv->getFeature(4) != 288 ) ok = false;	
		if(!ok){
			std::cout << "Unexpected value for non-ring C-C-N,C-C-N,C-X-C,N-X-X" << std::endl;
			pass = false;			
		}
	}
	delete fv;

	//Empty triples (set feature indicating no pairs)
	nl = createMolPtr("C");
	initMolProps(nl);
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), null_atom );
	
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 ){ 
			std::cout << "Missing feature flagging no triples" << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;
}

FeaturesTestGasteigerCharges::FeaturesTestGasteigerCharges(){
	description = "Test Gasteiger Charges feature";
}

void FeaturesTestGasteigerCharges::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("GasteigerCharges"); //IonRootTriples should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	RDKit::Atom *null_atom = NULL;

	//Non-ring (-0.33042488,-0.00652530)
	romol_ptr_t ion = createMolPtr("C");
	initMolProps(ion);
	ion.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.33042488 );
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	romol_ptr_t nl = createMolPtr("C");
	initMolProps(nl);	
	nl.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.00652530 );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );
	nl.get()->setProp("IsRingBreak",0);

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 9 ){		
			std::cout << "Unexpected value for ion gasteiger charge (non-ring) " << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Ring Retain Order (-0.01,0.9), (0.05,-0.4)
	ion = createMolPtr("CC");
	initMolProps(ion);	
	ion.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", -0.01 );
	ion.get()->getAtomWithIdx(1)->setProp<double>("OrigGasteigerCharge", 0.05 );
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );	
	nl = createMolPtr("CC");
	initMolProps(nl);
	nl.get()->getAtomWithIdx(0)->setProp<double>("OrigGasteigerCharge", 0.9 );
	nl.get()->getAtomWithIdx(1)->setProp<double>("OrigGasteigerCharge", -0.4 );
	nl.get()->setProp("IsRingBreak",1);	
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );
	

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 54 ){		
			std::cout << "Unexpected value for ion gasteiger charge ring " << fv->getFeature(1) << std::endl;
			pass = false;
		}
		if( fv->getFeature(2) != 56 ){		
			std::cout << "Unexpected value for nl gasteiger charge ring " << fv->getFeature(2) << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}

FeaturesTestHydrogenMovement::FeaturesTestHydrogenMovement(){
	description = "Test Hydrogen Movement feature";
}

void FeaturesTestHydrogenMovement::runTest(){
	
	bool pass = true;
	RDKit::Atom *null_atom = NULL;
	std::vector<std::string> fnames;
	fnames.push_back("HydrogenMovement"); //IonRootTriples should work the same
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Positive within range
	romol_ptr_t ion = createMolPtr("C");
	initMolProps(ion);
	double h_movement = 3.00452;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	romol_ptr_t nl = createMolPtr("C");
	initMolProps(nl);
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 8 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Negative within range
	h_movement = -1.115;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 4 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Zero
	h_movement = 0.0;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 5 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Out of range positive
	h_movement = 5.0;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 10 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Out of range negative
	h_movement = -5.0;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 10 ){		
			std::cout << "Unexpected idx for hydrogen movement" << fv->getFeature(1) << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}

FeaturesTestFunctionalGroups::FeaturesTestFunctionalGroups(){
	description = "Test Functional Groups features";
}

void FeaturesTestFunctionalGroups::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonFunctionalGroupFeatures");
	fnames.push_back("NLFunctionalGroupFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("CCCCC(O)=O");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[2], 0);	//Break Bond 2
	node->generateChildrenOfBreak(breaks[2]);

	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 9 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(1) != 8 ||		//Ion: 7,9 and 86 are at the root atom
			fv->getFeature(2) != 10 ||
			fv->getFeature(3) != 87 ||
			fv->getFeature(4) != 163 ||		//Ion: 0,10,87 and 114 are one away
			fv->getFeature(5) != 173 ||
			fv->getFeature(6) != 246 ||
			fv->getFeature(7) != 270 ||
			fv->getFeature(8) != 486 ){		//No functional group on NL side
			std::cout << "Unexpected Functional Group features" << std::endl;
			pass = false;
		}
	}

	passed = pass;
}

FeaturesTestExtraFunctionalGroups::FeaturesTestExtraFunctionalGroups(){
	description = "Test Extra Functional Groups features";
}

void FeaturesTestExtraFunctionalGroups::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonExtraFunctionalGroupFeatures");
	fnames.push_back("NLExtraFunctionalGroupFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("C1CC1CC(C)(C)C");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[0], 0);	//Break Bond 2
	node->generateChildrenOfBreak(breaks[0]);

	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 6 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(1) != 3 ||
			fv->getFeature(2) != 15 ||
			fv->getFeature(3) != 30 ||
			fv->getFeature(4) != 35 ||		
			fv->getFeature(5) != 45  ){		
			std::cout << "Unexpected Extra Functional Group features" << std::endl;
			pass = false;
		}
	}

	passed = pass;
}


FeaturesTestFunctionalGroupsRootOnly::FeaturesTestFunctionalGroupsRootOnly(){
	description = "Test Functional Groups root only features";
}

void FeaturesTestFunctionalGroupsRootOnly::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonFunctionalGroupRootOnlyFeatures");
	fnames.push_back("NLFunctionalGroupRootOnlyFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("CCCCC(O)=O");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[2], 0);	//Break Bond 2
	node->generateChildrenOfBreak(breaks[2]);

	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 5 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(1) != 8 ||		//Ion: 7,9 and 86 are at the root atom
			fv->getFeature(2) != 10 ||
			fv->getFeature(3) != 87 ||
			fv->getFeature(4) != 324 ){		//No functional group on NL side
			std::cout << "Unexpected Functional Group features" << std::endl;
			pass = false;
		}
	}

	passed = pass;
}


FeaturesTestRadicalFeatures::FeaturesTestRadicalFeatures(){
	description = "Test Radical features";
}

void FeaturesTestRadicalFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("RadicalFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//Ion Radical, NL Non-radical (1)
	romol_ptr_t ion = createMolPtr("CC[CH3+]");
	romol_ptr_t nl = createMolPtr("CC(=O)O");
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(2) );
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(2) );

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 ){		
			std::cout << "Unexpected features for ion radical specification " << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Ion Non-Radical, NL Radical (2)
	ion = createMolPtr("CCC[CH4+]");
	nl = createMolPtr("[CH]=C");
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 2 ){		
			std::cout << "Unexpected features for nl radical specification " << std::endl;
			pass = false;
		}
	}
	delete fv;

	//Neither radical (3)
	ion = createMolPtr("CCC[CH4+]");
	nl = createMolPtr("CCC");
	rtd_ion = RootedROMolPtr( ion, ion.get()->getAtomWithIdx(0), ion.get()->getAtomWithIdx(1) );
	rtd_nl = RootedROMolPtr( nl, nl.get()->getAtomWithIdx(0), nl.get()->getAtomWithIdx(1) );

	fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 2 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 3 ){		
			std::cout << "Unexpected features for non-radical specification " << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}

FeaturesTestRingFeatures::FeaturesTestRingFeatures(){
	description = "Test Ring features";
}

void FeaturesTestRingFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("RingFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	//1,1,3,6
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("c1ccc2ccc(CCC)cc2c1");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[4], 0);	//Break Ring with two bonds at distance 3
	node->generateChildrenOfBreak(breaks[4]);

	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 5 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 2 || fv->getFeature(2) != 3 || fv->getFeature(3) != 6 || fv->getFeature(4) != 11 ){		
			std::cout << "Unexpected features for ring break 1,1,3,6 " << std::endl;
			pass = false;
		}
	}
	delete fv;

	//0,0,3,7
	std::string smiles_or_inchi2("C1CCCCCC1");
	FragmentTreeNode *node2 = fgen.createStartNode( smiles_or_inchi2, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks2;
	node2->generateBreaks(breaks2, false);
	node2->applyBreak(breaks2[3], 0);	//Break Ring with two bonds at distance 3
	node2->generateChildrenOfBreak(breaks2[3]);

	child = &(node2->children[0]);
	Transition tmp_t2( -1, -1, child->nl, child->ion ); 
	fv = fc->computeFV( tmp_t2.getIon(), tmp_t2.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 || fv->getFeature(2) != 6 || fv->getFeature(3) != 12 ){		
			std::cout << "Unexpected features for ring break 0,0,3,7 " << std::endl;
			pass = false;
		}
	}
	delete fv;
	delete fc;
	passed = pass;

}

FeaturesTestExtraRingFeatures::FeaturesTestExtraRingFeatures(){
	description = "Test Extra Ring features";
}

void FeaturesTestExtraRingFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("ExtraRingFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );

	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("c1ccccc1-C2CCCCC2-CCC");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);

	//Break single bond between two rings
	node->applyBreak(breaks[0], 0);
	node->generateChildrenOfBreak(breaks[0]);
	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 || fv->getFeature(2) != 2 || fv->getFeature(3) != 3 ){		
			std::cout << "Unexpected extra ring features " << std::endl;
			pass = false;
		}
	}
	delete fv;
	node->undoBreak(breaks[0], 0);
	node->children = std::vector<FragmentTreeNode>();

	//Break single bond between ring and non-ring
	node->applyBreak(breaks[1], 0);	
	node->generateChildrenOfBreak(breaks[1]);
	child = &(node->children[0]);
	Transition tmp_t2( -1, -1, child->nl, child->ion ); 
	fv = fc->computeFV( tmp_t2.getIon(), tmp_t2.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	else{
		if( fv->getFeature(1) != 1 || fv->getFeature(2) != 3 ){		
			std::cout << "Unexpected extra ring features " << std::endl;
			pass = false;
		}
	}
	delete fv;
	node->undoBreak(breaks[1], 0);
	node->children = std::vector<FragmentTreeNode>();

	//Break Ring (no features set)
	node->applyBreak(breaks[5], 0);	
	node->generateChildrenOfBreak(breaks[5]);
	child = &(node->children[0]);
	Transition tmp_t3( -1, -1, child->nl, child->ion ); 
	fv = fc->computeFV( tmp_t3.getIon(), tmp_t3.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 1 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;		
	} 
	delete fv;

	passed = pass;
}


FeaturesTestRootMMFFAtomType::FeaturesTestRootMMFFAtomType(){
	description = "Test root MMFF Atom Type features";
}

void FeaturesTestRootMMFFAtomType::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonRootMMFFAtomType");
	fnames.push_back("NLRootMMFFAtomType");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("CCCCC(O)=O");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[4], 0);	//Break Bond 4
	node->generateChildrenOfBreak(breaks[4]);

	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 3 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(1) != 6 || fv->getFeature(2) != 103 ){		
			std::cout << "Unexpected MMFF atom types: expecting 6,103" << std::endl;
			pass = false;
		}
	}

	passed = pass;
}

FeaturesTestNeighbourMMFFAtomType::FeaturesTestNeighbourMMFFAtomType(){
	description = "Test neighbour MMFF Atom Type features";
}

void FeaturesTestNeighbourMMFFAtomType::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonNeighbourMMFFAtomType");
	fnames.push_back("NLNeighbourMMFFAtomType");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("CCCCC(O)=O");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[4], 0);	//Break Bond 4
	node->generateChildrenOfBreak(breaks[4]);

	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 
	FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(1) != 101 || fv->getFeature(2) != 102 || fv->getFeature(3) != 108  ){		
			std::cout << "Unexpected MMFF atom types: expecting 101, 102, 108" << std::endl;
			pass = false;
		}
	}

	passed = pass;
}

FeaturesTestBrokenOrigBondType::FeaturesTestBrokenOrigBondType(){
	description = "Test orig broken bond type features";
}

void FeaturesTestBrokenOrigBondType::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("BrokenOrigBondType");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("[Na+].C#CCCC(=O)CC=CC=CC=CCCc1ccccc1");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, true);
	
	typedef std::pair<int, int> bond_test_spec_t;	//Break_idx, Expected Bond Type
	std::vector<bond_test_spec_t> test_cases;
	test_cases.push_back( bond_test_spec_t(0, 6) );	//IONIC
	test_cases.push_back( bond_test_spec_t(31, 7) );	//H ONLY
	test_cases.push_back( bond_test_spec_t(1, 3) );	//TRIPLE
	test_cases.push_back( bond_test_spec_t(2, 1) );	//SINGLE
	test_cases.push_back( bond_test_spec_t(5, 2) );	//DOUBLE
	test_cases.push_back( bond_test_spec_t(8, 5) );	//CONJUGATED
	test_cases.push_back( bond_test_spec_t(19, 4) );//AROMATIC
	
	std::vector<bond_test_spec_t>::iterator testit = test_cases.begin();
	for( ; testit != test_cases.end(); ++testit ){
	
		Break *brk = &breaks[testit->first];
		node->applyBreak(*brk, 0);	//Break specified bond
		node->generateChildrenOfBreak(*brk);

		FragmentTreeNode *child = &(node->children[0]);
		Transition tmp_t( -1, -1, child->nl, child->ion ); 
		FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
		if( fv->getNumSetFeatures() != 2 ){
			std::cout << "Unexpected number of non-zero features" << std::endl;
			pass = false;	
		}
		else{
			if( fv->getFeature(1) != testit->second ){		
				std::cout << "Unexpected broken bond type: expecting " << testit->second << " but found " << fv->getFeature(1) << std::endl;
				pass = false;
			}
		}
		node->children = std::vector<FragmentTreeNode>();
	}

	passed = pass;
}

FeaturesTestNeighbourOrigBondType::FeaturesTestNeighbourOrigBondType(){
	description = "Test neighbour orig bond type features";
}

void FeaturesTestNeighbourOrigBondType::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("NeighbourOrigBondTypes");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("CC(c1ccccc1)=CC(C)C#C");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, true);
	
	typedef std::pair<int, std::vector<int> > bond_test_spec_t;	//Break_idx, Expected FV
	std::vector<bond_test_spec_t> test_cases;
	int exp_fv1[] = { 6, 7 }; //Conjugated on ion side, no bonds on nl side
	test_cases.push_back( bond_test_spec_t(0, std::vector<int>( exp_fv1, exp_fv1+2 )) );	
	int exp_fv2[] = { 5, 8, 12 };	//Aromatic on ion side, single and conjugated on nl side
	test_cases.push_back( bond_test_spec_t(1, std::vector<int>( exp_fv2, exp_fv2+3 )) );	
	int exp_fv3[] = { 4, 8 };	//Triple on ion side, single on nl side
	test_cases.push_back( bond_test_spec_t(5, std::vector<int>( exp_fv3, exp_fv3+2 )) );	

	std::vector<bond_test_spec_t>::iterator testit = test_cases.begin();
	for( ; testit != test_cases.end(); ++testit ){
	
		Break *brk = &breaks[testit->first];
		node->applyBreak(*brk, 0);	//Break specified bond
		node->generateChildrenOfBreak(*brk);

		FragmentTreeNode *child = &(node->children[0]);	//Take first child of break
		Transition tmp_t( -1, -1, child->nl, child->ion ); 
		FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
		if( fv->getNumSetFeatures() != testit->second.size()+1 ){
			std::cout << "Unexpected number of non-zero features:" << fv->getNumSetFeatures() << std::endl;
			pass = false;	
		}
		else{
			std::vector<int>::iterator it = testit->second.begin();
			for( int idx = 1; it != testit->second.end(); ++it, idx++ ){
				if( fv->getFeature(idx) != *it ){		
					std::cout << "Unexpected feature in neighbour bond type: expecting " << *it << " but found " << fv->getFeature(idx) << std::endl;
					pass = false;
				}
			}
		}
		node->children = std::vector<FragmentTreeNode>();
		node->undoBreak(*brk,0);
	}

	passed = pass;
}

FeaturesTestRootAtom::FeaturesTestRootAtom(){
	description = "Test root atom features";
}

void FeaturesTestRootAtom::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonRootAtom");
	fnames.push_back("NLRootAtom");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("[Br]C([Cl])(F)N(I)OPS[Se][Si]([Na])CCCCCc1ccccc1");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, true);
	
	typedef std::pair<int, int> atom_pair_t;
	typedef std::pair<int, atom_pair_t> atom_test_spec_t;	//Break_idx, Expected Atom feature idxs
	std::vector<atom_test_spec_t> test_cases;
	test_cases.push_back( atom_test_spec_t(breaks.size()-1, atom_pair_t(1,11)) );	//C-H
	test_cases.push_back( atom_test_spec_t(0, atom_pair_t(1,0)) );		//C-Br
	test_cases.push_back( atom_test_spec_t(1, atom_pair_t(2,1)) );		//C-Cl
	test_cases.push_back( atom_test_spec_t(2, atom_pair_t(3,1)) );		//C-F
	test_cases.push_back( atom_test_spec_t(3, atom_pair_t(5,1)) );		//C-N
	test_cases.push_back( atom_test_spec_t(4, atom_pair_t(4,5)) );		//N-I
	test_cases.push_back( atom_test_spec_t(5, atom_pair_t(6,5)) );		//N-O
	test_cases.push_back( atom_test_spec_t(6, atom_pair_t(7,6)) );	//O-P
	test_cases.push_back( atom_test_spec_t(7, atom_pair_t(8,7)) );	//P-S
	test_cases.push_back( atom_test_spec_t(8, atom_pair_t(9,8)) );	//S-Se
	test_cases.push_back( atom_test_spec_t(9, atom_pair_t(10,9)) );	//Se-Si
	test_cases.push_back( atom_test_spec_t(10, atom_pair_t(12,10)) );	//Si-Na
	test_cases.push_back( atom_test_spec_t(19, atom_pair_t(1,1)) );	//Ring c-c c-c
	
	std::vector<atom_test_spec_t>::iterator testit = test_cases.begin();
	for( ; testit != test_cases.end(); ++testit ){
	
		Break *brk = &breaks[testit->first];
		node->applyBreak(*brk, 0);	//Break specified bond
		node->generateChildrenOfBreak(*brk);

		//if(brk->getBondIdx() >= 0 ){
		//	std::cout << brk->getBondIdx() << ": " << node->ion.get()->getBondWithIdx(brk->getBondIdx())->getBeginAtom()->getSymbol();
		//	std::cout << "-" << node->ion.get()->getBondWithIdx(brk->getBondIdx())->getEndAtom()->getSymbol() << std::endl;
		//}

		FragmentTreeNode *child = &(node->children[0]);
		Transition tmp_t( -1, -1, child->nl, child->ion ); 
		FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
		//std::cout << *(tmp_t.getNLSmiles()) << std::endl;
		if( fv->getNumSetFeatures() != 3 ){
			std::cout << "Unexpected number of non-zero features" << std::endl;
			pass = false;	
		}
		else{
			if( fv->getFeature(1) != testit->second.first + 1  ){				//Ion
				std::cout << "Unexpected Ion Root Atom Feature: expecting " << testit->second.first + 1 << " but found " << fv->getFeature(1) << std::endl;
				pass = false;
			}
			if( fv->getFeature(2) != (testit->second.second + 13 + 1)){	//NL
				std::cout << "Unexpected NL Root Atom Feature: expecting " << (testit->second.second + 13 + 1)<< " but found " << fv->getFeature(2) << std::endl;
				pass = false;			
			}
		}
		node->children = std::vector<FragmentTreeNode>();
		node->undoBreak(*brk, 0);
	}

	passed = pass;
}

FeaturesTestIonicFeatures::FeaturesTestIonicFeatures(){
	description = "Test ionic fragment features";
}

class IonicFeatureTestCase {
public:
	IonicFeatureTestCase(int a_break_idx, int a_ionic_idx, int a_child_idx, bool a_rebreak_child, std::vector<int> &a_expected_output) : 
	  break_idx( a_break_idx), ionic_idx( a_ionic_idx ), child_idx( a_child_idx ), expected_output( a_expected_output ), rebreak_child(a_rebreak_child) {};
	int break_idx;
	int ionic_idx;
	int child_idx;
	bool rebreak_child;
	std::vector<int> expected_output;
};

void FeaturesTestIonicFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");
	FeatureCalculator *fc = new FeatureCalculator( fnames );
	
	FragmentGraphGenerator fgen(fc);
	std::string smiles_or_inchi("[Na+].[Cl-].CC=CC");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, true);
	
	std::vector<IonicFeatureTestCase> test_cases;

	//Tests case with Na+ on ion side, Cl- on NL side
	std::vector<int> exp_result_t1; 
	exp_result_t1.push_back(2); exp_result_t1.push_back(3);
	test_cases.push_back( IonicFeatureTestCase(0,0,0,false,exp_result_t1) );

	//Tests case with Na+ and Cl- on NL side
	std::vector<int> exp_result_t2; 
	exp_result_t2.push_back(1); exp_result_t2.push_back(3);
	test_cases.push_back( IonicFeatureTestCase(0,1,0,false,exp_result_t2) );

	//Tests case with Na+ and Cl- on Ion side
	std::vector<int> exp_result_t3; 
	exp_result_t3.push_back(2); exp_result_t3.push_back(4);
	test_cases.push_back( IonicFeatureTestCase(2,0,0,false,exp_result_t3) );

	//Tests case with no ionic fragments (need to re-break a child from above)
	std::vector<int> exp_result_t4; 
	exp_result_t4.push_back(5);
	test_cases.push_back( IonicFeatureTestCase(0,1,0,true,exp_result_t4) );
	
	std::vector<IonicFeatureTestCase>::iterator testit = test_cases.begin();
	for( ; testit != test_cases.end(); ++testit ){
	
		Break *brk = &breaks[ testit->break_idx];

		node->applyBreak(*brk, testit->ionic_idx );	//Break specified bond
		node->generateChildrenOfBreak(*brk);
		FragmentTreeNode *child = &(node->children[testit->child_idx]);

		if( testit->rebreak_child ){
			std::vector<Break> child_breaks;
			child->generateBreaks(child_breaks, true);
			child->applyBreak(child_breaks[0], 0 );
			child->generateChildrenOfBreak(child_breaks[0]);
			child = &(child->children[0]);
		}

		Transition tmp_t( -1, -1, child->nl, child->ion ); 
		FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
		
		if( fv->getNumSetFeatures() != testit->expected_output.size()+1 ){
			std::cout << "Unexpected number of non-zero features. Expected " << testit->expected_output.size()+1 << " but found " <<  fv->getNumSetFeatures() << std::endl;
			pass = false;	
		}
		else{
			std::vector<int>::iterator it = testit->expected_output.begin();
			for( int idx = 1; it != testit->expected_output.end(); ++it, idx++ ){
				if( fv->getFeature(idx) != *it ){	
					std::cout << "Unexpected Ionic Feature at idx " << idx << ": expecting " << *it << " but found " << fv->getFeature(idx) << std::endl;
					pass = false;
				}
			}
		}

		node->children = std::vector<FragmentTreeNode>();
		node->undoBreak(*brk, testit->ionic_idx );
	}

	passed = pass;
}

FeaturesTestQuadraticFeatures::FeaturesTestQuadraticFeatures(){
	description = "Test quadratic features";
}

void FeaturesTestQuadraticFeatures::runTest(){
	
	bool pass = true;
	std::vector<std::string> fnames;
	fnames.push_back("BreakAtomPair"); 
	fnames.push_back("HydrogenMovement"); 
	fnames.push_back("QuadraticFeatures"); 
	FeatureCalculator *fc = new FeatureCalculator( fnames );	

	//Simple initial vector with 3 bits set (indexes: 0,3,80 )
	romol_ptr_t ion= createMolPtr("C");
	initMolProps(ion);	
	RDKit::Atom *null_atom = NULL;
	RootedROMolPtr rtd_ion( ion, ion.get()->getAtomWithIdx(0), null_atom );
	double h_movement = 3.00452;
	ion.get()->getAtomWithIdx(0)->setProp<double>("OriginalMass", 16.0 - h_movement);
	romol_ptr_t nl = createMolPtr("N");
	initMolProps(nl);	
	RootedROMolPtr rtd_nl( nl, nl.get()->getAtomWithIdx(0), null_atom );
	nl.get()->setProp("IsRingBreak",0);

	FeatureVector *fv = fc->computeFV(&rtd_ion, &rtd_nl);
	if( fv->getNumSetFeatures() != 4 ){
		std::cout << "Unexpected number of non-zero features" << std::endl;
		pass = false;	
	}
	else{
		if( fv->getFeature(0) != 0 || fv->getFeature(1) != 3 || fv->getFeature(2) != 80 ){
			std::cout << "Unexpected singular features" << std::endl;
			pass = false;				
		}
		if( fv->getFeature(3) != 3166 ){
			std::cout << "Unexpected quadratic feature:" << fv->getFeature(3) << std::endl;
			pass = false;				
		}
	}
	int total = fv->getTotalLength();
	if( total != 3404 ){
		std::cout << "Unexpected total feature count: " << total << std::endl;
		pass = false;		
	}
	int numfeatures = fc->getNumFeatures();
	if( numfeatures != 3404 ){
		std::cout << "Unexpected total feature calculation: " << numfeatures << std::endl;
		pass = false;		
	}
	delete fv;
	delete fc;
	passed = pass;
}


FeaturesTestLength::FeaturesTestLength(){
	description = "Test length of features";
}

void FeaturesTestLength::runTest(){

	bool pass = true;
	
	//Create the feature calculator
	std::vector<std::string> names = FeatureCalculator::getValidFeatureNames();

	//Create some aribitrary input data
	FeatureCalculator full_fc( names );
	FragmentGraphGenerator fgen(&full_fc);
	std::string smiles_or_inchi("CCCCC(O)=O");
	FragmentTreeNode *node = fgen.createStartNode( smiles_or_inchi, POSITIVE_ESI_IONIZATION_MODE );
	std::vector<Break> breaks;
	node->generateBreaks(breaks, false);
	node->applyBreak(breaks[2], 0);	//Break Bond 2
	node->generateChildrenOfBreak(breaks[2]);
	FragmentTreeNode *child = &(node->children[0]);
	Transition tmp_t( -1, -1, child->nl, child->ion ); 

	//Check all feature lengths
	std::vector<std::string>::const_iterator it = names.begin();
	for( ; it != names.end(); ++it ){
		
		std::vector<std::string> feature_list;
		feature_list.push_back( *it );
		FeatureCalculator fc( feature_list );

		//Compute the feature vector
		FeatureVector *fv = fc.computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );

		//Check the length
		if( fv->getTotalLength() != fc.getNumFeatures() ){ 
			std::cout << "Feature length incorrect for " << *it;
			std::cout << ": expecting " << fc.getNumFeatures();
			std::cout << " but found " << fv->getTotalLength() << std::endl; 
			pass = false;
		}
		delete fv;
	}
	passed = pass;
}

FeaturesTestMetlinExample::FeaturesTestMetlinExample(){
	description = "Test Metlin Example Features";
}

void FeaturesTestMetlinExample::runTest(){

	bool pass = true;
	std::string smiles_Metlin_21361 = "O=C(NC(CCC(N)=O)C(O)=O)C(C)NC(=O)C(N)CC(C)C";

	std::string config_filename = "tests/test_data/example_feature_config.txt";
	FeatureCalculator fc( config_filename );

	//Ingegration Test - compute the fragment tree and transitions
	std::string id = "Metlin_21361";
	config_t cfg; initDefaultConfig(cfg);
	cfg.fg_depth = 1; cfg.include_h_losses = true; 
	MolData mol( id, smiles_Metlin_21361, 0, &cfg );

	mol.computeFragmentGraph(&fc);
	mol.computeFeatureVectors(&fc);

	//Transition 0: Check that the features are as expected
	const FeatureVector *fv = mol.getFeatureVectorForIdx(0);
	if( fv->getNumSetFeatures() != 11 ){
		std::cout << "Unexpected number of non-zero features " << fv->getNumSetFeatures();
		pass = false;
	}
	else{
		if( fv->getFeature(0) != 0 ){ 
			std::cout << "Bias error: " << fv->getFeature(0) <<std::endl;
			pass = false;
		}
		if( fv->getFeature(1) != 5 ){ 
			std::cout << "Incorrect BreakAtomPair: " << fv->getFeature(1) << std::endl;
			pass = false;
		}
		if( fv->getFeature(2) != 98 ){
			std::cout << "Incorrect Gasteiger: " << fv->getFeature(2) << std::endl;
			pass = false;
		}
		if( fv->getFeature(3) != 148 ){
			std::cout << "Incorrect HMovement: " << fv->getFeature(3)<< std::endl;
			pass = false;
		}
		if( fv->getFeature(4) != 156){
			std::cout << "Incorrect IonRootPair C-C: " << fv->getFeature(4) << std::endl;
			pass = false;
		}
		if( fv->getFeature(5) != 160 ){
			std::cout << "Incorrect IonRootPair C-N: " << fv->getFeature(5) << std::endl;
			pass = false;
		}
		if( fv->getFeature(6) != 301 ){
			std::cout << "Incorrect IonRootTriple C-C-C: " << fv->getFeature(6) << std::endl;
			pass = false;
		}
		if( fv->getFeature(7) != 305 ){
			std::cout << "Incorrect IonRootTriple C-C-N: " << fv->getFeature(7) << std::endl;
			pass = false;
		}
		if( fv->getFeature(8) != 325 ){
			std::cout << "Incorrect IonRootTriple C-N-C: " << fv->getFeature(8)<< std::endl;
			pass = false;
		}
		if( fv->getFeature(9) != 1165 ){
			std::cout << "Incorrect NLRootPair (None): " << fv->getFeature(9) << std::endl;
			pass = false;
		}
		if( fv->getFeature(10) != 1310 ){
			std::cout << "Incorrect NLRootTriple (None): " << fv->getFeature(10) << std::endl;
			pass = false;
		}
	}

	passed = pass;
}
