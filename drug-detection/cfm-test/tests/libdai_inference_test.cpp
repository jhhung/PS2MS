/*#########################################################################
# Mass Spec Prediction of HMDB Molecules
#
# inference_test.cpp
#
# Description: Test code for inference.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "libdai_inference_test.h"
#include "Inference.h"
#include "MolData.h"
#include "Config.h"

#include "dai/factor.h"
#include "dai/factorgraph.h"
#include "dai/properties.h"
#include "dai/jtree.h"

class InferTestMol : public MolData {
public:
	InferTestMol() : MolData("Infer Test Mol", "", POSITIVE_ESI_IONIZATION_MODE){
	
		//Create a molecule based on what was previously in test_bn_transition_ipfp.txt
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		fg->addToGraph( FragmentTreeNode( createMolPtr("O=C(O)C[NH2+]C(=O)C(NC(=O)CN)Cc1ccccc1"), basic_nl, -1, -1, &fh), -1 ); //id = 0
		fg->addToGraph( FragmentTreeNode( createMolPtr("NCC(=O)[NH2+]C(=C=O)Cc1ccccc1"), basic_nl, -1, -1, &fh), 0 ); // id = 1, 0 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("N=CC(=O)[NH2+]CCc1ccccc1"), basic_nl, -1, -1, &fh), 0 ); //id = 2, 0 -> 2
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH2+]=CCc1ccccc1"), basic_nl, -1, -1, &fh), 2 ); // id = 3, 2 -> 3

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
			MolData::normalizeAndSortSpectrum(spectra[i]);

	}
};


class LibDaiTestMol : public MolData {
public:
	LibDaiTestMol() : MolData("Libdai Simple Test Mol", "", POSITIVE_ESI_IONIZATION_MODE){
	
		//Create a molecule based on what was previously in test_bn_transition.txt
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C") ;
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH3+]CCCN"), basic_nl, -1, -1, &fh), -1 ); //id = 0
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH4+]"), basic_nl, -1, -1, &fh), 0 ); // id = 1, 0 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("C=CC[NH3+]"), basic_nl, -1, -1, &fh), 0 ); //id = 2, 0 -> 2
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH4+]"), basic_nl, -1, -1, &fh), 2 ); // id = 1, 2 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("C=C=[NH2+]"), basic_nl, -1, -1, &fh), 2 ); // id = 3, 2 -> 3
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH3+]"), basic_nl, -1, -1, &fh), 2 ); // id = 3, 2 -> 4
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH3+]"), basic_nl, -1, -1, &fh), 3 ); // id = 3, 3 -> 4
		fg->addToGraph( FragmentTreeNode( createMolPtr("[C+]#N"), basic_nl, -1, -1, &fh), 3 ); // id = 5, 3 -> 5

		//Set a uniformly random theta value between -3.0 and 3.0 for each transition
		thetas.resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			thetas[energy].resize( fg->getNumTransitions() );
			for( int i = 0; i < fg->getNumTransitions(); i++ )
				thetas[energy][i] = (double(std::rand())/double(RAND_MAX) - 0.5)*6;
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(30.42, 2.82) );
		spectra[0].push_back( Peak(58.21 , 90.00 ) );
		spectra[0].push_back( Peak(75.13, 30.00) );
		spectra[1].push_back( Peak(74.89, 45.83) );
		spectra[2].push_back( Peak(42.36, 1234.0) );
		spectra[2].push_back( Peak(43.38, 12.79) );
		spectra[2].push_back( Peak(58.08, 1234.0) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);

	}
};

class LibDaiComplexTestMol : public MolData {
public:
	LibDaiComplexTestMol() : MolData("Libdai Complex Test Mol", "CCOP(O)(=S)OCC", POSITIVE_ESI_IONIZATION_MODE){
	
		//Create a fragmentation graph based on what was previously in HMDB01460_bn_transition.txt
		computeFragmentGraph(2);

		//Set a uniformly random theta value between -3.0 and 3.0 for each transition
		thetas.resize( 3 );
		for( int energy = 0; energy < 3; energy++ ){
			thetas[energy].resize( fg->getNumTransitions() );
			for( int i = 0; i < fg->getNumTransitions(); i++ )
				thetas[energy][i] = (double(std::rand())/double(RAND_MAX) - 0.5)*6;
		}
		computeTransitionProbabilities();

		//Set the spectra
		spectra.resize(3);
		spectra[0].push_back( Peak(115.00, 100.00) );
		spectra[0].push_back( Peak(143.00 , 62.53 ) );
		spectra[0].push_back( Peak(171.00, 15.11) );
		spectra[1].push_back( Peak(81.00, 31.04) );
		spectra[1].push_back( Peak(96.00, 5.54) );
		spectra[1].push_back( Peak(97.00, 100.00) );
		spectra[1].push_back( Peak(114.00, 6.68) );
		spectra[1].push_back( Peak(115.00, 56.63) );
		spectra[2].push_back( Peak(63.00, 10.29) );
		spectra[2].push_back( Peak(65.00 , 31.86) );
		spectra[2].push_back( Peak(79.00, 29.32) );
		spectra[2].push_back( Peak(81.00, 77.84) );
		spectra[2].push_back( Peak(82.00, 6.94) );
		spectra[2].push_back( Peak(96.00, 6.10) );
		spectra[2].push_back( Peak(97.00, 100.0) );
		spectra[2].push_back( Peak(115.00, 16.43) );
		for( int i = 0; i <=2; i++ ) 
			MolData::normalizeAndSortSpectrum(spectra[i]);

	}
};

bool compareWithLibdai( MolData &moldata, config_t &cfg );
void initialiseFactors( std::vector<dai::Factor> &factors, MolData &moldata, config_t &cfg );
void setFactorGraphEvidence( dai::FactorGraph &fg, std::vector<dai::Factor> &factors, MolData &moldata, config_t &cfg, int energy );

InferenceTestSimpleDepth3::InferenceTestSimpleDepth3(){
	description = "Test inference on simple input at depth 3";
}

void InferenceTestSimpleDepth3::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 9;
	int tmp_array[3] = {3,6,9};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	LibDaiTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;

}

InferenceTestComplexDepth3::InferenceTestComplexDepth3(){
	description = "Test inference on more complex input at depth 3";
}

void InferenceTestComplexDepth3::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 9;
	int tmp_array[3] = {3,6,9};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	LibDaiComplexTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;

}

InferenceTestvsIPFPDepth2::InferenceTestvsIPFPDepth2(){
	description = "Test inference on infer v ipfp input at depth 2";
}

void InferenceTestvsIPFPDepth2::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	InferTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;

}

InferenceTestSimpleDepth2::InferenceTestSimpleDepth2(){
	description = "Test inference on simple input at depth 2";
}

void InferenceTestSimpleDepth2::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	LibDaiTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;

}

InferenceTestComplexDepth2::InferenceTestComplexDepth2(){
	description = "Test inference on more complex input at depth 2";
}

void InferenceTestComplexDepth2::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 6;
	int tmp_array[3] = {2,4,6};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	LibDaiComplexTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;

}

InferenceTestSimpleDepth1::InferenceTestSimpleDepth1(){
	description = "Test inference on simple input at depth 1";
}

void InferenceTestSimpleDepth1::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 3;
	int tmp_array[3] = {1,2,3};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	LibDaiTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;

}

InferenceTestComplexDepth1::InferenceTestComplexDepth1(){
	description = "Test inference on more complex input at depth 1";
}

void InferenceTestComplexDepth1::runTest(){
	
	bool pass = true;
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.model_depth = 3;
	int tmp_array[3] = {1,2,3};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.use_single_energy_cfm = 1;
	initDerivedConfig(cfg);

	LibDaiComplexTestMol moldata;
	for( int energy = 0; energy < 3; energy ++ ){
		config_t se_cfg;
		initSingleEnergyConfig(se_cfg, cfg, energy);
		pass &= compareWithLibdai(moldata, se_cfg);
	}
	passed = pass;
}

bool compareWithLibdai( MolData &moldata, config_t &cfg ){

	bool pass = true;
	double tol = 1e-12;

	//Run our inference
	Inference infer( &moldata, &cfg ); 
    beliefs_t infer_beliefs;
	infer.calculateBeliefs(infer_beliefs);

	//Run LibDai Inference
	std::vector<dai::Factor> factors;
	initialiseFactors( factors, moldata, cfg );
	dai::FactorGraph fg( factors );	
	setFactorGraphEvidence( fg, factors, moldata, cfg, cfg.map_d_to_energy[0] );
	dai::PropertySet opts;
	dai::JTree bp(fg, opts("updates",std::string("HUGIN"))("logdomain",true));
	bp.init();
	bp.run();

	const FragmentGraph *frg = moldata.getFragmentGraph();

	//Compare the results
	for( unsigned int i = 0; i < frg->getNumTransitions(); i++ ){

		const Transition *t = frg->getTransitionAtIdx(i);

		//F1 Marginal Probabilities
		if( t->getFromId() == 0 ){	//main ion is always id = 0
			double my_belief = exp( infer_beliefs.tn[i][0] );
			double dai_belief = bp.belief(fg.var(1))[t->getToId()];
			if( fabs( my_belief - dai_belief ) > tol ){
				std::cout << "Mismatch for F1 marginal in transition " << i << " ";
				std::cout << my_belief << " vs " << dai_belief << std::endl;
				pass = false;
			} 
		}

		//F1,F2 and F2,F3 Marginal Joint Probabilities
		for( unsigned int j = 2; j <= cfg.model_depth; j++ ){
			double my_belief = exp( infer_beliefs.tn[i][j-1] );
			double dai_belief = bp.belief(fg.factor(j).vars())[t->getToId() * frg->getNumFragments() + t->getFromId()];
			if( fabs( my_belief - dai_belief ) > tol ){
				std::cout << "Mismatch for " << j << " joint marginal in transition " << i << " ";
				std::cout << my_belief << " vs " << dai_belief << std::endl;
				pass = false;
			} 
		}

	}

	//Persistence
	for( unsigned int i = 0; i < frg->getNumFragments(); i++ ){
		
		//F1 Marginal Probabilities
		if( i == 0 ){	//main ion is always id = 0
			double my_belief = exp( infer_beliefs.ps[i][0] );
			double dai_belief = bp.belief(fg.var(1))[i];
			if( fabs( my_belief - dai_belief ) > tol ){
				std::cout << "Mismatch for F1 persistence marginal in transition " << i << " ";
				std::cout << my_belief << " vs " << dai_belief << std::endl;
				pass = false;
			} 
		}

		//F1,F2 and F2,F3 Marginal Joint Probabilities
		for( unsigned int j = 2; j <= cfg.model_depth; j++ ){
			double my_belief = exp( infer_beliefs.ps[i][j-1] );
			double dai_belief = bp.belief(fg.factor(j).vars())[i * frg->getNumFragments() + i];
			if( fabs( my_belief - dai_belief ) > tol ){
				std::cout << "Mismatch for " << j << " persistence joint marginal in transition " << i << " ";
				std::cout << my_belief << " vs " << dai_belief << std::endl;
				pass = false;
			} 
		}
	}

	return pass;
}


void createFactorProbabilities( std::vector<double> &output, MolData &moldata, int energy);

//Create a vector of factors forming the tree
void initialiseFactors( std::vector<dai::Factor> &factors, MolData &moldata, config_t &cfg ){

	unsigned int i;
	std::vector<dai::Var> variables;
	int numfragments = moldata.getFragmentGraph()->getNumFragments();

	//Create the variable for the fragment nodes
	for( i = 0; i <= cfg.model_depth; i++ )
		variables.push_back(dai::Var( i, numfragments ));

	//Create the variable for the peak node
	variables.push_back(dai::Var( cfg.model_depth + 1,  2 ));

	//Initial Factor
	std::vector<double> probs;
	probs.assign( numfragments, 0 );
	probs[0] = 1;	//Set the Main Ion to 100%
	dai::VarSet vars; 
	vars |= variables[0];
	factors.push_back( dai::Factor( vars, probs ) );

	int energy = cfg.map_d_to_energy[0];

	//Create the remaining factors
	for( i = 1; i <= cfg.model_depth + 1; i++ ){
		
		//Create the Varset
		dai::VarSet vars; 
		vars |= variables[i-1];
		vars |= variables[i];

		//Add the factor
		std::vector<double> probs;
		createFactorProbabilities( probs, moldata, energy );
		
		if( i < cfg.model_depth + 1 ) factors.push_back( dai::Factor( vars, probs ) );
		else factors.push_back( dai::Factor( vars ) );	//Peak factor
	}

}

//Set the parameters of the factor graph according to the current probabilities
void createFactorProbabilities( std::vector<double> &output, MolData &moldata, int energy){


	unsigned int i;
	const FragmentGraph *fg =  moldata.getFragmentGraph();
	unsigned int numfragments = fg->getNumFragments();
	unsigned int numtransitions = fg->getNumTransitions();
	
	//Note: Would be better if we could somehow enter the log probabilities
	//directly into libdai without converting back to exp?

	//Inialise all transition probabilities to 0
	int num_states = numfragments * numfragments;
	output.assign( num_states, 0 );

	//Now set the non-zero transition (break) probabilities
	for( i = 0; i < numtransitions; i++ ){
		
		const Transition *t = fg->getTransitionAtIdx(i);
		int idx = t->getToId() * numfragments + t->getFromId();
		output[idx] = exp(moldata.getLogTransitionProbForIdx(energy, i));

	}

	//Now set the persistence probabilities (i.e. no break)
	for( i = 0; i < numfragments; i++ ){
		
		int idx = i * numfragments + i;
		output[idx] = exp(moldata.getLogPersistenceProbForIdx(energy,i));
	}

}

void setFactorGraphEvidence( dai::FactorGraph &fg, std::vector<dai::Factor> &factors, MolData &moldata, config_t &cfg, int energy ){

	double prob_sum = 0.0;
	std::vector<double> probs;
	const FragmentGraph *frg = moldata.getFragmentGraph();
	unsigned int numfragments = frg->getNumFragments();
	probs.resize(numfragments);

	//Note: better to set values directly ui

	//Set the main ion (always the first fragment)
	fg.clamp( 0, 0 );

	//Set the spectrum information
	Message msg;
	Inference infer( &moldata, &cfg );
	infer.createSpectrumMessage( msg, energy );
	for( unsigned int i = 0; i < numfragments; i++ ){
		double prob = std::exp(msg.getIdx(i));
		factors[cfg.model_depth + 1].set( i , prob );
		factors[cfg.model_depth + 1].set( i + numfragments, 1 - prob );
	}
	fg.setFactor( cfg.model_depth + 1, factors[cfg.model_depth + 1]);
	fg.clamp( cfg.model_depth + 1, 0 );
	
}