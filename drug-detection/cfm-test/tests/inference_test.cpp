/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# ipfp_test.cpp
#
# Description: Test code for IPFP.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "inference_test.h"
#include "Config.h"
#include "FragmentTreeNode.h"
#include "IPFP.h"

#include <boost/filesystem.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>

std::vector<int> null_eloc;

class InferTestMol : public MolData {
public:
	InferTestMol(config_t *cfg) : MolData("Infer Test Mol", "", cfg){
	
		//Create a molecule based on what was previously in test_bn_transition_ipfp.txt
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		fg->addToGraph( FragmentTreeNode( createMolPtr("O=C(O)C[NH2+]C(=O)C(NC(=O)CN)Cc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc), -1 ); //id = 0
		fg->addToGraph( FragmentTreeNode( createMolPtr("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc), 0 ); // id = 1, 0 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("N=CC(=O)[NH2+]CCc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc), 0 ); //id = 2, 0 -> 2
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH2+]=CCc1ccccc1"), basic_nl, -1, -1, &fh, null_eloc), 2 ); // id = 3, 2 -> 3

		//Set thetas/transition probs to match matlab reference2
		thetas.resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			thetas[energy].resize( fg->getNumTransitions() );
			thetas[energy][0] = 1.386294361119891;	//0->1
			thetas[energy][1] = 2.456735772821304;	//0->2
			thetas[energy][2] = 0.0;	//2->3
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(120.082039, 35.914500 ) );
		spectra[0].push_back( Peak(177.103613, 50.000000 ) );
		spectra[0].push_back( Peak(280.133689, 50.146650) );
		spectra[1].push_back( Peak(120.083820, 100.00000 ) );
		spectra[1].push_back( Peak(177.106284, 33.000500) );
		spectra[2].push_back( Peak(120.081802, 100.00000) );
		for( int i = 0; i <=2; i++ ) 
			spectra[i].normalizeAndSort();

	}
};

class InferComplexTestMol : public MolData {
public:
	InferComplexTestMol(config_t *cfg) : MolData("NIST2011_1201", "CC(C)C(N=C(O)CN)C(=O)O", cfg){
		
		//Compulte a fragmentation graph
		std::vector<std::string> fnames;
		fnames.push_back("BreakAtomPair");
		FeatureCalculator fc(fnames);
		computeFragmentGraphAndReplaceMolsWithFVs(&fc);

		//Randomly set the theta values
		thetas.resize(1);
		thetas[0].resize( fg->getNumTransitions() );
		for( int i = 0; i < fg->getNumTransitions(); i++ )
			thetas[0][i] = (double(std::rand())/double(RAND_MAX) -  0.5)*3;
		computeTransitionProbabilities();

		//Load a very simple spectrum
		spectra.resize(1);
		spectra[0].push_back( Peak(114.0, 100.0 ) );
		spectra[0].push_back( Peak(174.0, 100.0 ) );
		spectra[0].normalizeAndSort();

	}
};

class SpectrumTestMol : public MolData {
public:
	SpectrumTestMol(config_t *cfg) : MolData("Spectrum Test Mol", "", cfg){
	
		//Create a molecule based on what was previously in test_bn_transition_ipfp.txt
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		//Root node with Mass1 (0), two fragments with Mass2 (1,2), three fragments with Mass3 (3,4,5)
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCCCCCO"), basic_nl, -1, -1, &fh, null_eloc), -1 ); 
		fg->addToGraph( FragmentTreeNode( createMolPtr("CC(C)C"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCCC"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCOC"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CC(C)O"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCCO"), basic_nl, -1, -1, &fh, null_eloc), 0 );

		//Set the spectra
		spectra.resize(1);
		spectra[0].push_back( Peak(102.10446506800000, 10.0 ) );
		spectra[0].push_back( Peak(58.078250319999995, 30.0 ) );
		spectra[0].push_back( Peak(60.057514875999999, 60.0 ) );
		spectra[0].normalizeAndSort();

	}
};

class SpectrumTestMolNoisePeak : public MolData {
public:
	SpectrumTestMolNoisePeak(config_t *cfg) : MolData("Spectrum Test Mol with Noise Peak", "", cfg){
	
		//Create a molecule  as above, but add a peak quite far from any fragment mass (noise peak)
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		//Root node with Mass1 (0), two fragments with Mass2 (1,2), three fragments with Mass3 (3,4,5)
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCCCCCO"), basic_nl, -1, -1, &fh, null_eloc), -1 ); 
		fg->addToGraph( FragmentTreeNode( createMolPtr("CC(C)C"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCCC"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCOC"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CC(C)O"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCCO"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("CCOO"), basic_nl, -1, -1, &fh, null_eloc), 0 );

		//Set the spectra
		spectra.resize(1);
		spectra[0].push_back( Peak(102.10446506800000, 10.0 ) );
		spectra[0].push_back( Peak(58.078250319999995, 30.0 ) );
		spectra[0].push_back( Peak(60.057514875999999, 60.0 ) );
		spectra[0].push_back( Peak(70.0, 50.0 ) );
		spectra[0].normalizeAndSort();

	}
};

class IsotopeSpectrumTestMolNoisePeak : public MolData {
public:
	IsotopeSpectrumTestMolNoisePeak(config_t *cfg) : MolData("Isotope Spectrum Test Mol with Noise Peak", "", cfg){
	
		//Create a molecule as for isotope test above, but add a peak quite far from any fragment mass (noise peak)
		fg = new FragmentGraph(cfg);
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		//Two molecules, (almost) same mass but different isotope patterns.
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH2-]C(C)NC1=CN=CC(=C1C#N)Cl"), basic_nl, -1, -1, &fh, null_eloc), -1 ); 
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH-](C1=NNN=N1)NC2=NC(=O)C(=O)N2"), basic_nl, -1, -1, &fh, null_eloc), 0 );
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH-](C1=NNN=N1)NC"), basic_nl, -1, -1, &fh, null_eloc), 0 );

		//Set the spectra
		spectra.resize(1);
		spectra[0].push_back( Peak(194.0490426799, 100.0000000000 ) );
		spectra[0].push_back( Peak(195.0518167914, 11.3124853142 ) );
		spectra[0].push_back( Peak(196.0462415099, 32.9809130321 ) );
		spectra[0].push_back( Peak(197.0489073693, 3.6831305376) );
		spectra[0].push_back( Peak(50.0, 100.0) );
		spectra[0].normalizeAndSort();

	}
};

class IsotopeSpectrumTestMol : public MolData {
public:
	IsotopeSpectrumTestMol(config_t *cfg) : MolData("Isotope Spectrum Test Mol", "", cfg){
	
		//Create a molecule based on what was previously in test_bn_transition_ipfp.txt
		fg = new FragmentGraph(cfg);
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		//Two molecules, (almost) same mass but different isotope patterns.
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH2-]C(C)NC1=CN=CC(=C1C#N)Cl"), basic_nl, -1, -1, &fh, null_eloc), -1 ); 
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH-](C1=NNN=N1)NC2=NC(=O)C(=O)N2"), basic_nl, -1, -1, &fh, null_eloc), 0 );

		//Set the spectra
		spectra.resize(1);
		spectra[0].push_back( Peak(194.0490426799, 100.0000000000 ) );
		spectra[0].push_back( Peak(195.0518167914, 11.3124853142 ) );
		spectra[0].push_back( Peak(196.0462415099, 32.9809130321 ) );
		spectra[0].push_back( Peak(197.0489073693, 3.6831305376) );
		spectra[0].normalizeAndSort();

	}
};


bool runInferenceTestCompareVsIPFP( MolData &moldata, config_t &cfg ){

	bool pass = true;
	double tol = 1e-8;

	//Run the inference code
	Inference infer( &moldata, &cfg ); 
    beliefs_t infer_beliefs;
	infer.calculateBeliefs(infer_beliefs);
    
	//Run code for IPFP_OSC
	cfg.ipfp_algorithm = IPFP_WITH_MOD_ALGORITHM;
	IPFP ipfp( &moldata, &cfg ); 
    beliefs_t *ipfp_beliefs = ipfp.calculateBeliefs();

	//Compare the two
	const FragmentGraph *frg = moldata.getFragmentGraph();

	//Compare the results
	for( unsigned int i = 0; i < frg->getNumTransitions(); i++ ){

		const Transition *t = frg->getTransitionAtIdx(i);

		//F1 Marginal Probabilities
		if( t->getFromId() == 0 ){	//main ion is always id = 0
			double infer_belief = exp( infer_beliefs.tn[i][0] );
			double ipfp_belief = exp( ipfp_beliefs->tn[i][0] );
			if( fabs( infer_belief - ipfp_belief ) > tol ){
				std::cout << "Mismatch for F1 marginal in transition " << i << " ";
				std::cout << infer_belief << " vs " << ipfp_belief << std::endl;
				pass = false;
			} 
		}

		//F1,F2 and F2,F3 Marginal Joint Probabilities
		for( unsigned int j = 2; j <= cfg.model_depth; j++ ){
			double infer_belief = exp( infer_beliefs.tn[i][j-1] );
			double ipfp_belief = exp( ipfp_beliefs->tn[i][j-1] );
			if( fabs( infer_belief - ipfp_belief ) > tol ){
				std::cout << "Mismatch for " << j << " joint marginal in transition " << i << " ";
				std::cout << infer_belief << " vs " << ipfp_belief << std::endl;
				pass = false;
			} 
		}

	}

	//Persistence
	for( unsigned int i = 0; i < frg->getNumFragments(); i++ ){
		
		//F1 Marginal Probabilities
		if( i == 0 ){	//main ion is always id = 0
			double infer_belief = exp( infer_beliefs.ps[i][0] );
			double ipfp_belief =  exp( ipfp_beliefs->ps[i][0] );
			if( fabs( infer_belief - ipfp_belief ) > tol ){
				std::cout << "Mismatch for F1 persistence marginal in transition " << i << " ";
				std::cout << infer_belief << " vs " << ipfp_belief << std::endl;
				pass = false;
			} 
		}

		//F1,F2 and F2,F3 Marginal Joint Probabilities
		for( unsigned int j = 2; j <= cfg.model_depth; j++ ){
			double infer_belief = exp( infer_beliefs.ps[i][j-1] );
			double ipfp_belief = exp( ipfp_beliefs->ps[i][j-1] );
			if( fabs( infer_belief - ipfp_belief ) > tol ){
				std::cout << "Mismatch for " << j << " persistence joint marginal in transition " << i << " ";
				std::cout << infer_belief << " vs " << ipfp_belief << std::endl;
				pass = false;
			} 
		}
	}

	return pass;

}


InferenceTestCompareVsIPFPSingleEnergyCase::InferenceTestCompareVsIPFPSingleEnergyCase(){
	description = "Test computing of inference beliefs to compare against IPFP for single energy case";
}

void InferenceTestCompareVsIPFPSingleEnergyCase::runTest(){

	bool pass = true;

	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,2,2};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	double tmp_array2[3] = {1.0,1.0,1.0};
	cfg.spectrum_weights.assign(tmp_array2, tmp_array2+3);
	initDerivedConfig(cfg);

	config_t se_cfg;
	initSingleEnergyConfig(se_cfg, cfg, 1);

	//Load some molecule data
	InferTestMol moldata(&se_cfg);
		
	//Run the test
	pass = runInferenceTestCompareVsIPFP( moldata, se_cfg );

	passed = pass;
}

InferenceTestCompareVsIPFPSharedMassCase::InferenceTestCompareVsIPFPSharedMassCase(){
	description = "Test computing of inference beliefs to compare against IPFP for a case with shared masses";
}

void InferenceTestCompareVsIPFPSharedMassCase::runTest(){

	bool pass = true;

	//Create the config
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 2;
	cfg.spectrum_depths.push_back(2);
	cfg.spectrum_weights.push_back(1);
	cfg.abs_mass_tol = 0.5;
	initDerivedConfig(cfg);

	//Load a test molecule (NIST2011_1201)
	InferComplexTestMol moldata(&cfg);
	
	//Run the test
	pass = runInferenceTestCompareVsIPFP( moldata, cfg );

	passed = pass;
}


InferenceTestSpectrumMessage::InferenceTestSpectrumMessage(){
	description = "Test computing of spectrum message (no isotopes)";
}

void InferenceTestSpectrumMessage::runTest(){

	bool pass = true;
	double tol = 0.00001;

	config_t cfg;
	initDefaultConfig(cfg);

	SpectrumTestMol moldata(&cfg);

	Inference infer( &moldata, &cfg);
	Message msg, down_msg;
	down_msg.reset(6); for(int i=0; i<6; i++ ) down_msg.addToIdx(i, 0.0);
	infer.createSpectrumMessage( msg, 0, down_msg );

	if( fabs( std::exp(msg.getIdx(0)) - 0.1 ) > tol ||
		fabs( std::exp(msg.getIdx(1)) - 0.15 ) > tol ||
		fabs( std::exp(msg.getIdx(2)) - 0.15 ) > tol ||
		fabs( std::exp(msg.getIdx(3)) - 0.2 ) > tol ||
		fabs( std::exp(msg.getIdx(4)) - 0.2 ) > tol ||
		fabs( std::exp(msg.getIdx(5)) - 0.2 ) > tol
	){
		std::cout << "Unexpected value for spectrum message" << std::endl;
			for( int i = 0; i < 6; i++ ) std::cout << std::exp(msg.getIdx(i)) << " ";
			pass = false;
	}
	passed = pass;
}

InferenceTestSpectrumMessageNoisePeak::InferenceTestSpectrumMessageNoisePeak(){
	description = "Test computing of spectrum message with the presence of a noise peak";
}

void InferenceTestSpectrumMessageNoisePeak::runTest(){

	bool pass = true;
	double tol = 0.00001;

	config_t cfg;
	initDefaultConfig(cfg);

	SpectrumTestMolNoisePeak moldata(&cfg);

	Inference infer( &moldata, &cfg);
	Message msg, down_msg;
	down_msg.reset(7); for(int i=0; i<7; i++ ) down_msg.addToIdx(i, 0.0);
	infer.createSpectrumMessage( msg, 0, down_msg );

	if( fabs( std::exp(msg.getIdx(0)) - 0.1 ) > tol ||
		fabs( std::exp(msg.getIdx(1)) - 0.15 ) > tol ||
		fabs( std::exp(msg.getIdx(2)) - 0.15 ) > tol ||
		fabs( std::exp(msg.getIdx(3)) - 0.2 ) > tol ||
		fabs( std::exp(msg.getIdx(4)) - 0.2 ) > tol ||
		fabs( std::exp(msg.getIdx(5)) - 0.2 ) > tol ||
		fabs( std::exp(msg.getIdx(6)) - 0.0 ) > tol
	){
		std::cout << "Unexpected value for spectrum message" << std::endl;
			for( int i = 0; i < 7; i++ ) std::cout << std::exp(msg.getIdx(i)) << " ";
			pass = false;
	}
	passed = pass;
}


InferenceTestSpectrumMessageWithIsotopes::InferenceTestSpectrumMessageWithIsotopes(){
	description = "Test computing of spectrum message with isotopes.";
}

void InferenceTestSpectrumMessageWithIsotopes::runTest(){

	bool pass = true;
	double tol = 0.00001;

	config_t cfg;
	initDefaultConfig(cfg);
	cfg.abs_mass_tol = 0.1;
	cfg.ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;
	cfg.include_isotopes = true; cfg.isotope_thresh = 0.01;

	IsotopeSpectrumTestMol moldata(&cfg);

	Inference infer( &moldata, &cfg);
	Message msg, down_msg;
	down_msg.reset(2); for(int i=0; i<2; i++ ) down_msg.addToIdx(i, 0.0);
	infer.createSpectrumMessage( msg, 0, down_msg );

	if( fabs( std::exp(msg.getIdx(0)) - 0.567333 ) > tol ||
		fabs( std::exp(msg.getIdx(1)) - 0.432667 ) > tol
		){
			std::cout << "Unexpected value for isotope spectrum message" << std::endl;
				for( int i = 0; i < 2; i++ ) std::cout << std::exp(msg.getIdx(i)) << " ";
				pass = false;
	}

	passed = pass;
}

InferenceTestSpectrumMessageWithIsotopesAndNoisePeak::InferenceTestSpectrumMessageWithIsotopesAndNoisePeak(){
	description = "Test computing of spectrum message with isotopes and noise peak.";
}

void InferenceTestSpectrumMessageWithIsotopesAndNoisePeak::runTest(){

	bool pass = true;
	double tol = 0.00001;

	config_t cfg;
	initDefaultConfig(cfg);
	cfg.abs_mass_tol = 0.1;
	cfg.ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;
	cfg.fg_depth = 2;
	cfg.include_isotopes = true; cfg.isotope_thresh = 0.01;

	IsotopeSpectrumTestMolNoisePeak moldata(&cfg);

	Inference infer( &moldata, &cfg);
	Message msg, down_msg;
	down_msg.reset(3); for(int i=0; i<3; i++ ) down_msg.addToIdx(i, 0.0);
	infer.createSpectrumMessage( msg, 0, down_msg );

	if( fabs( std::exp(msg.getIdx(0)) - 0.567333 ) > tol ||
		fabs( std::exp(msg.getIdx(1)) - 0.432667 ) > tol ||
		fabs( std::exp(msg.getIdx(2)) - 0.0 ) > tol
		){
			std::cout << "Unexpected value for isotope spectrum message" << std::endl;
				for( int i = 0; i < 3; i++ ) std::cout << std::exp(msg.getIdx(i)) << " ";
				pass = false;
	}

	passed = pass;
}