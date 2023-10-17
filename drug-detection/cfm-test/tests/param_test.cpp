/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param_test.cpp
#
# Description: Test code for Param.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "mpi.h"

#include "param_test.h"
#include "MolData.h"
#include "Features.h"
#include "Config.h"
#include "EM.h"

#include <boost/filesystem.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>


ParamsTestComputeTransitionThetas::ParamsTestComputeTransitionThetas(){
	description = "Test computing of transition thetas";
}

void ParamsTestComputeTransitionThetas::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	std::vector<std::string> fnames;
	fnames.push_back("BreakAtomPair");
	Param param(fnames, 3);

	//Set some parameter weights
	param.setWeightAtIdx(0.234, 0);
	param.setWeightAtIdx(-3.124, 1);
	param.setWeightAtIdx(0.89, 2);

	int med_offset = param.getNumWeightsPerEnergyLevel();
	param.setWeightAtIdx(-5.23432, 0 + med_offset);
	param.setWeightAtIdx(0.0, 1 + med_offset);
	param.setWeightAtIdx(1.0, 2 + med_offset);

	int high_offset = 2*med_offset;
	param.setWeightAtIdx(-35.0, 0 + high_offset);
	param.setWeightAtIdx(-10.0, 1 + high_offset);
	param.setWeightAtIdx( 0.0, 2 + high_offset);

	//Create some Feature Vectors
	FeatureVector fv0, fv1;
	fv0.addFeatureAtIdx(0.0, med_offset-1);	//Set the feature length
	fv1.addFeatureAtIdx(0.0, med_offset-1);
	fv0.addFeatureAtIdx(1.0, 1);
	fv0.addFeatureAtIdx(1.0, 2);
	fv1.addFeatureAtIdx(1.0, 0);

	//Compute the transition thetas
	double theta0_low = param.computeTheta(fv0, 0);
	double theta0_med = param.computeTheta(fv0, 1);
	double theta0_high = param.computeTheta(fv0, 2);
	double theta1_low = param.computeTheta(fv1, 0);
	double theta1_med = param.computeTheta(fv1, 1);
	double theta1_high = param.computeTheta(fv1, 2);

	if( fabs(theta0_low - -2.234 )  > tol ){
		std::cout << "Mismatch in Transition 0 Low Energy Theta" << std::endl;
		pass = false;
	}

	if( fabs(theta0_med - 1.0 )  > tol ){
		std::cout << "Mismatch in Transition 0 Med Energy Theta" << std::endl;
		pass = false;
	}

	if( fabs(theta0_high - -10.0 )  > tol ){
		std::cout << "Mismatch in Transition 0 High Energy Theta" << std::endl;
		pass = false;
	}

	if( fabs(theta1_low - 0.234 )  > tol ){
		std::cout << "Mismatch in Transition 1 Low Energy Theta" << std::endl;
		pass = false;
	}

	if( fabs(theta1_med - -5.23432 )  > tol ){
		std::cout << "Mismatch in Transition 1 Med Energy Theta" << std::endl;
		pass = false;
	}

	if( fabs(theta1_high - -35.0 )  > tol ){
		std::cout << "Mismatch in Transition 1 High Energy Theta" << std::endl;
		pass = false;
	}

	passed = pass;

}

std::vector<int> param_null_eloc;

class ParamTestMol : public MolData {
public:
	ParamTestMol(config_t *cfg) : MolData("Param Test Mol", "", cfg){

		//Create a molecule based on what was previously in test_bn_transition.txt
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C");
		fg->addToGraph( FragmentTreeNode( createMolPtr("NCCCN"), basic_nl,-1, -1, &fh, param_null_eloc), -1 ); //id = 0
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH4+]"), basic_nl,-1, -1, &fh, param_null_eloc), 0 ); // id = 1, 0 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("C=CC[NH3+]"), basic_nl,-1, -1, &fh, param_null_eloc), 0 ); //id = 2, 0 -> 2
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH4+]"), basic_nl,-1, -1, &fh, param_null_eloc), 2 ); // 2 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("C=C=[NH2+]"), basic_nl,-1, -1, &fh, param_null_eloc), 2 );  //id = 3, 2 -> 3
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH3+]"), basic_nl,-1, -1, &fh, param_null_eloc), 2 ); //id = 4, 2->4
		fg->addToGraph( FragmentTreeNode( createMolPtr("[CH3+]"), basic_nl, -1, -1, &fh, param_null_eloc),3 ); //3->4
		fg->addToGraph( FragmentTreeNode( createMolPtr("[C+]#N"), basic_nl,-1, -1, &fh, param_null_eloc), 3 ); //id = 5, 3->5
		graph_computed = true;

		//Add a simple feature vector to each transition
		unsigned int num_trans = fg->getNumTransitions();
		fvs.resize( num_trans ); 
		for( unsigned int i = 0; i < num_trans; i++ ){
			fvs[i] = new FeatureVector();
			fvs[i]->addFeature(1.0);			//Bias term
			fvs[i]->addFeatureAtIdx(0.0, 10);	//Resize to 11 (Using HydrogenMovement feature to initialise params)
		}
		fvs[1]->addFeatureAtIdx(1.0,1);	//Transition 0->2
		fvs[3]->addFeatureAtIdx(1.0,1); //Transition 2->3
		fvs[4]->addFeatureAtIdx(1.0,1); //Transition 2->4
	}
};


ParamsTestComputeAndAccumulateGradient::ParamsTestComputeAndAccumulateGradient(){
	description = "Test gradient computation and accumulation";
}

void ParamsTestComputeAndAccumulateGradient::runTest(){
	
	bool pass = true;
	double tol = 1e-3;
    
	std::vector<double> grads;

	suft_counts_t suft;

	//Model Config
	config_t cfg;
	cfg.model_depth = 3;
	int tmp_array[3] = {1,2,3};
	cfg.spectrum_depths.assign(tmp_array, tmp_array+3);
	cfg.lambda = 0.01;
	initDerivedConfig(cfg);

	//Create the molecule data
	ParamTestMol moldata(&cfg);

	//Set some arbitrary parameter weights (only used in the regularizer term and for sizing)
	std::vector<std::string> fnames;
	fnames.push_back("HydrogenMovement");	//Size = 10
	Param param(fnames, 3);
	for( unsigned int i = 0; i < param.getNumWeights(); i++ ) param.setWeightAtIdx(0.0, i);
	int med_offset = param.getNumWeightsPerEnergyLevel();
	int high_offset = 2*med_offset;
	param.setWeightAtIdx(0.2, 0);
	param.setWeightAtIdx(-0.5, 1);
	param.setWeightAtIdx(0.25, med_offset);
	param.setWeightAtIdx(0.1, high_offset);
	param.setWeightAtIdx(-0.1, 1+high_offset);
	
	std::string param_filename = "tmp_param_file.log";
	if( boost::filesystem::exists( param_filename ) )
		boost::filesystem::remove( param_filename );
	param.saveToFile( param_filename );

	//Set some arbitrary suft values
	suft.values.resize(1);

	const FragmentGraph *fg = moldata.getFragmentGraph();
	unsigned int N = fg->getNumTransitions() + fg->getNumFragments();
	suft.values[0].assign(3*N, 0.0);
	suft.values[0][1] = 0.5;
	suft.values[0][7] = 0.5;
	suft.values[0][1 + N] = 0.4;
	suft.values[0][7 + N] = 0.1;
	suft.values[0][3 + N] = 0.3;
	suft.values[0][9 + N] = 0.2;
	suft.values[0][1 + 2*N] = 0.2;
	suft.values[0][7 + 2*N] = 0.05;
	suft.values[0][3 + 2*N] = 0.7;
	suft.values[0][9 + 2*N] = 0.025;
	suft.values[0][10 + 2*N] = 0.025;

	//Initialise all gradients to 0
	grads.resize(33);
	for( unsigned int i = 0; i < grads.size(); i++ ) grads[i] = 0.0;

	std::set<unsigned int> used_idxs;
	FeatureCalculator fc_null(fnames); 
	std::string null_str = "null";
	EM em(&cfg, &fc_null, null_str, param_filename );
	double Q_only = em.computeQ( 0, moldata, suft );
	double Q = em.computeAndAccumulateGradient(&grads[0], 0, moldata, suft, true, used_idxs);
	
	//Check Q
	if( fabs(Q - -3.823 )  > tol ){
		std::cout << "Unexpected Q value resulting from gradient computation: " << Q << std::endl;
		pass = false;
	}
	if( fabs(Q_only - -3.823 )  > tol ){
		std::cout << "Unexpected Q value resulting from Q only computation: " << Q << std::endl;
		pass = false;
	}

	//Check the gradients
	double expected_vals[6] = {-0.1624,0.2499,-0.0568,0.2554,0.1649,0.4663};
	for( int energy = 0; energy < 3; energy++ ){
		for( unsigned int i = 0; i < 2; i++ ){
			if( fabs(grads[energy*11 + i] - expected_vals[energy*2+i])  > tol ){
				std::cout << "Unexpected grad[" << energy*11 + i << "] value resulting from gradient computation: " << grads[energy*11 + i] << std::endl;
				pass = false;
			}
		}
	}

	//Check the used flags
	for( int energy = 0; energy < 3; energy++ ){
		//Check the ones that should be on
		for( unsigned int i = 0; i < 2; i++ ){
			if( used_idxs.find(energy*11 + i) == used_idxs.end()){
				std::cout << "Expected used_flag to be set but found not set" << std::endl;
				pass = false;
			}
		}
		//Check one that shouldn't be on
		if( used_idxs.find(energy*11 + 5) != used_idxs.end()){
			std::cout << "Expected used_flag to be not set but found set" << std::endl;
			pass = false;
		}
	}
	passed = pass;

}
