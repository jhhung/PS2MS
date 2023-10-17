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

#include "nn_param_test.h"
#include "MolData.h"
#include "Features.h"
#include "Config.h"
#include "EM_NN.h"

#include <boost/filesystem.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>


class NNTestParam : public NNParam {
public:
	NNTestParam(std::vector<std::string> a_feature_list, int a_num_energy_levels, std::vector<int> &a_hlayer_num_nodes, std::vector<int> &a_act_func_ids) : 
	  NNParam(a_feature_list, a_num_energy_levels, a_hlayer_num_nodes, a_act_func_ids) {	

		randomInit();	//There will be a bunch of weights we don't set below

		//Set the parameter weights
		weights[0] = 1.0;  weights[2] = -2.0; weights[3] = -4.0; weights[5] = 3.0;   //H1
		weights[6] = -3.0; weights[8] = 5.0; weights[9] = 12.0;  weights[11] = -7.0; //H2
		weights[12]= -2.0; weights[14] = 1.0; weights[15] = 6.0; weights[17] = -3.0; //H3
		weights[18] = 5.0; weights[20] = 1.0; weights[21] = -2.0; weights[23]= 2.0;	  //H4

		weights[24] = 4.0; weights[25] = 2.0; weights[26] = 1.0; weights[27] = -3.0; //H5
		weights[29] = -5.0; weights[30] = 1.0; weights[31] = 2.0; weights[32] = -2.0; //H6

		weights[34] = 5.0; weights[35] = -2.0; weights[36] = -1.0;	//Theta
	}
};


NNParamsTestComputeTransitionThetas::NNParamsTestComputeTransitionThetas(){
	description = "Test computing of neural net transition thetas";
}

void NNParamsTestComputeTransitionThetas::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");	//5 features + bias =  6 features
	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 4; hlayer_numnodes[1] = 2; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear

	NNTestParam param(fnames, 1, hlayer_numnodes, act_ids);

	//Create the Feature Vector
	FeatureVector fv;
	fv.addFeatureAtIdx(1.0, 0); fv.addFeatureAtIdx(1.0, 2); fv.addFeatureAtIdx(1.0, 5);

	//Compute the theta
	double theta = param.computeTheta(fv, 0);

	if( fabs(theta - 12.0 )  > tol ){
		std::cout << "Expected theta of 12.0 but found theta of " << theta << std::endl;
		pass = false;
	}

	passed = pass;

}

NNParamsTestSaveAndLoadFromFile::NNParamsTestSaveAndLoadFromFile(){
	description = "Test saving and loading of neural net parameters to/from file";
}

void NNParamsTestSaveAndLoadFromFile::runTest(){
	
	bool pass = true;
    double tol = 1e-6;

	//Create some parameters
	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");
	fnames.push_back("RingFeatures");	//18 features in total

	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 12; hlayer_numnodes[1] = 4; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear?

	NNParam param(fnames, 1, hlayer_numnodes, act_ids);
	param.randomInit();

	//Save to file
	std::string filename = "tmp_param_file.log";
	if( boost::filesystem::exists( filename ) )
		boost::filesystem::remove( filename );
	param.saveToFile( filename );

	//Load from file
	NNParam param_load( filename );

	//Check the configuration
	if( param.getNumEnergyLevels() != param_load.getNumEnergyLevels() ){
		std::cout << "Mismatch in number of energy levels in saved and loaded parameters" << std::endl;
		pass = false;
	}

	//Check the feature names
	std::vector<std::string> *fnames1 = param.getFeatureNames();
	std::vector<std::string> *fnames2 = param_load.getFeatureNames();
	if( fnames1->size() != fnames2->size() ){
		std::cout << "Mismatch in feature names length in saved and loaded parameters" << std::endl;
		pass = false;		
	}
	else{
		for( int i = 0; i < fnames1->size(); i++ ){
			if( (*fnames1)[i] != (*fnames2)[i] ){
				std::cout << "Mismatch in feature names in saved and loaded parameters: " << (*fnames1)[i] << " vs " << (*fnames2)[i] << std::endl;
				pass = false;
			}
		}
	}

	//Check the weights
	if( param.getNumWeights() != param_load.getNumWeights() ){
		std::cout << "Mismatch in number of weights in saved and loaded parameters" << std::endl;
		pass = false;			
	}
	else{
		for( int i = 0; i < param.getNumWeights(); i++ ){
			if( fabs( param.getWeightAtIdx(i) - param_load.getWeightAtIdx(i) ) > tol ){
				std::cout << "Mismatch in weight in saved and loaded parameters at idx " << i << std::endl;
				pass = false;	
				break;
			}
		}
	}

	//Check that a computed theta value is the same
	FeatureVector fv;
	fv.addFeatureAtIdx(1.0, 0); fv.addFeatureAtIdx(1.0, 2); fv.addFeatureAtIdx(1.0, 5);
	fv.addFeatureAtIdx(1.0, 8); fv.addFeatureAtIdx(1.0, 11); fv.addFeatureAtIdx(1.0, 13);
	fv.addFeatureAtIdx(1.0, 14); fv.addFeatureAtIdx(1.0, 16); fv.addFeatureAtIdx(1.0, 17);

	double theta1 = param.computeTheta( fv, 0 );
	double theta2 = param_load.computeTheta( fv, 0 );
	if( fabs( theta1 - theta2 ) > tol ){
		std::cout << "Mismatch in thetas computed with saved and loaded parameters: " << theta1 << " vs " << theta2 << std::endl;
		pass = false;		
	}

	passed = pass;
}




NNParamsTestComputeDeltas::NNParamsTestComputeDeltas(){
	description = "Test computing of neural net delta values";
}

void NNParamsTestComputeDeltas::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");	//5 features + bias =  6 features
	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 4; hlayer_numnodes[1] = 2; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear

	NNTestParam testparam(fnames, 1, hlayer_numnodes, act_ids);

	//Create the Feature Vector
	FeatureVector fv1, fv2;
	fv1.addFeatureAtIdx(1.0, 0); fv1.addFeatureAtIdx(1.0, 2); fv1.addFeatureAtIdx(1.0, 5);
	fv2.addFeatureAtIdx(1.0, 0); fv2.addFeatureAtIdx(1.0, 3); fv2.addFeatureAtIdx(1.0, 5);

	for( unsigned int energy = 0; energy <= 2; energy++ ){

		NNParam *param;
		if( energy == 0 ) param = new NNTestParam(fnames, 1, hlayer_numnodes, act_ids);
		else{ 
			param = new NNParam(fnames, energy, hlayer_numnodes, act_ids);
			param->appendNextEnergyParams( testparam );
		}

		//Run the forwards paths
		std::vector<azd_vals_t> z_values(2), a_values(2);
		double theta1 = param->computeTheta(fv1, energy, z_values[0], a_values[0]);
		double theta2 = param->computeTheta(fv2, energy, z_values[1], a_values[1]);
		double rho_denom = 1.0 + exp(theta1) + exp(theta2);

		//Now run backwards to compute the delta values
		std::vector<azd_vals_t> deltasA, deltasB;
		param->computeDeltas(deltasA, deltasB, z_values, a_values, rho_denom, energy);

		//Check the values
		double x1 = exp(theta1)/rho_denom, x2 = exp(theta2)/rho_denom;
		double exp_deltaA_1[] = { -5.0, -4.0, 0.0, 0.0, -2.0, -1.0, 1.0 };
		double exp_deltaB_1[] = { -5*x1, -4*x1, 0.0, 0.0, -2*x1, -x1, x1 };
		double exp_deltaA_2[] = { 0.0, 0.0, 8.0, 0.0, -2.0, -1.0, 1.0 };
		double exp_deltaB_2[] = { 0.0, 0.0, 8*x2, 0.0, -2*x2, -x2, x2 };

		if( deltasA.size() != 2 || deltasB.size() != 2 || deltasA[0].size() != 7 || deltasB[0].size() != 7 || deltasA[1].size() != 7 || deltasB[1].size() != 7 ){
			std::cout << "Unexpected size of delta output" << std::endl;
			pass = false;
		}
		else{
			for( int i = 0; i < 7; i++ ){
				if( fabs(exp_deltaA_1[i] - deltasA[0][i]) > tol ||
					fabs(exp_deltaB_1[i] - deltasB[0][i]) > tol ||
					fabs(exp_deltaA_2[i] - deltasA[1][i]) > tol ||
					fabs(exp_deltaB_2[i] - deltasB[1][i]) > tol ){
						std::cout << "Unexpected delta values at index " << i << std::endl;
						pass = false;			
						break;
				}
			}
	
		}
		delete param;
	}
	passed = pass;

}

NNParamsTestComputeUnweightedGradients::NNParamsTestComputeUnweightedGradients(){
	description = "Test computing of neural net unweighted gradients";
}

void NNParamsTestComputeUnweightedGradients::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");	//5 features + bias =  6 features
	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 4; hlayer_numnodes[1] = 2; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear

	NNTestParam testparam(fnames, 1, hlayer_numnodes, act_ids);

	//Create the Feature Vector
	FeatureVector fv1, fv2;
	fv1.addFeatureAtIdx(1.0, 0); fv1.addFeatureAtIdx(1.0, 2); fv1.addFeatureAtIdx(1.0, 5);
	fv2.addFeatureAtIdx(1.0, 0); fv2.addFeatureAtIdx(1.0, 3); fv2.addFeatureAtIdx(1.0, 5);

	for( unsigned int energy = 0; energy <= 2; energy++ ){

		NNParam *param;
		if( energy == 0 ) param = new NNTestParam(fnames, 1, hlayer_numnodes, act_ids);
		else{ 
			param = new NNParam(fnames, energy, hlayer_numnodes, act_ids);
			param->appendNextEnergyParams( testparam );
		}

		//Run the forwards paths
		std::vector<azd_vals_t> z_values(2), a_values(2);
		double theta1 = param->computeTheta(fv1, energy, z_values[0], a_values[0]);
		double theta2 = param->computeTheta(fv2, energy, z_values[1], a_values[1]);
		double rho_denom = 1.0 + exp(theta1) + exp(theta2);

		//Now run backwards to compute the delta values
		std::vector<azd_vals_t> deltasA, deltasB;
		param->computeDeltas(deltasA, deltasB, z_values, a_values, rho_denom, energy);

		//Now compute the unweighted gradients
		std::vector<std::vector<double> > unweighted_grads;
		std::set<unsigned int> used_idxs;
		std::vector<const FeatureVector *> fvs(2); fvs[0] = &fv1; fvs[1] = &fv2;
		param->computeUnweightedGradients( unweighted_grads, used_idxs, fvs, deltasA, deltasB, a_values );

		//Check the values
		if( unweighted_grads.size() != 3 || unweighted_grads[0].size() != 37 || unweighted_grads[1].size() != 37 || unweighted_grads[2].size() != 37 ){
			std::cout << "Unexpected size of unweighted gradients" << std::endl;
			pass = false;
		}
		else{
		
			double x1 = exp(theta1)/rho_denom, x2 = exp(theta2)/rho_denom;
			double exp_unweighted_1[] = {-5+5*x1,0,-5+5*x1,0,0,-5+5*x1,-4+4*x1,0,-4+4*x1,0,0,-4+4*x1,-8*x2, 0, 0, -8*x2,0,-8*x2,0,0,0,0,0,0,-2+2*x1+2*x2,-2*2+2*2*x1, 5*2-5*2*x1,2*x2,0,-1+x1+x2,-2+2*x1,5-5*x1,x2,0,1-x1-x2,3-3*x1-x2, -13+13*x1+7*x2 };
			double exp_unweighted_2[] = {5*x1,0,5*x1,0,0,5*x1, 4*x1,0,4*x1,0,0,4*x1, 8-8*x2, 0,0,8-8*x2,0,8-8*x2, 0,0,0,0,0,0, -2+2*x1+2*x2,4*x1,-5*2*x1,-2+2*x2,0,-1+x1+x2, 2*x1,-5*x1, -1+x2, 0, 1-x1-x2,+1-3*x1-x2,-7+13*x1+7*x2 };
			double exp_unweighted_persist[] = {5*x1,0,5*x1,0,0,5*x1, 4*x1,0,4*x1,0,0,4*x1,-8*x2, 0,0,-8*x2,0,-8*x2,0,0,0,0,0,0,2*x1+2*x2,4*x1,-5*2*x1,2*x2,0,x1+x2, 2*x1,-5*x1, x2, 0, -x1-x2,-3*x1-x2,13*x1+7*x2 };

			for( int i = 0; i < 37; i++ ){
				if( fabs( exp_unweighted_1[i] - unweighted_grads[0][i]) > tol ||
					fabs( exp_unweighted_2[i] - unweighted_grads[1][i]) > tol ||
					fabs( exp_unweighted_persist[i] - unweighted_grads[2][i]) > tol ){
						std::cout << "Unexpected unweighted gradient values at index " << i << std::endl;
						pass = false;			
						break;	
				}
			}
		}

		//Check the used indexes
		if( used_idxs.size() != 16 || 
			used_idxs.find(0) == used_idxs.end() || used_idxs.find(2) == used_idxs.end() || used_idxs.find(3) == used_idxs.end() || used_idxs.find(5) == used_idxs.end() ||
			used_idxs.find(6) == used_idxs.end() || used_idxs.find(8) == used_idxs.end() || used_idxs.find(9) == used_idxs.end() || used_idxs.find(11) == used_idxs.end() ||
			used_idxs.find(12) == used_idxs.end() || used_idxs.find(14) == used_idxs.end() || used_idxs.find(15) == used_idxs.end() || used_idxs.find(17) == used_idxs.end() ||
			used_idxs.find(18) == used_idxs.end() || used_idxs.find(20) == used_idxs.end() || used_idxs.find(21) == used_idxs.end() || used_idxs.find(23) == used_idxs.end() ){
				std::cout << "Unexpected used indexes" << std::endl;
				pass = false;			
		}

		delete param;
	}

	passed = pass;
}

std::vector<int> nn_param_null_eloc;

class NNParamTestMol : public MolData {
public:
	NNParamTestMol(config_t *cfg) : MolData("NN Param Test Mol", "", cfg){


		//Create two transitions with the feature vectors as in previous tests
		fg = new FragmentGraph();
		FeatureHelper fh;
		romol_ptr_t basic_nl = createMolPtr("C");
		fg->addToGraph( FragmentTreeNode( createMolPtr("NCCCN"), basic_nl,-1, -1, &fh, nn_param_null_eloc), -1 ); //id = 0
		fg->addToGraph( FragmentTreeNode( createMolPtr("[NH4+]"), basic_nl,-1, -1, &fh, nn_param_null_eloc), 0 ); // id = 1, 0 -> 1
		fg->addToGraph( FragmentTreeNode( createMolPtr("C=CC[NH3+]"), basic_nl,-1, -1, &fh, nn_param_null_eloc), 0 ); //id = 2, 0 -> 2
		graph_computed = true;

		fvs.resize( 2 ); 
		fvs[0] = new FeatureVector(); fvs[0]->addFeatureAtIdx(1.0, 0); fvs[0]->addFeatureAtIdx(1.0, 2); fvs[0]->addFeatureAtIdx(1.0, 5);
		fvs[1] = new FeatureVector(); fvs[1]->addFeatureAtIdx(1.0, 0); fvs[1]->addFeatureAtIdx(1.0, 3); fvs[1]->addFeatureAtIdx(1.0, 5);

	}
};


NNParamsTestComputeAndAccumulateGradient::NNParamsTestComputeAndAccumulateGradient(){
	description = "Test neural net gradient computation and accumulation";
}

void NNParamsTestComputeAndAccumulateGradient::runTest(){
	
	bool pass = true;
	double tol = 1e-3;
    
	std::vector<double> grads;

	suft_counts_t suft;

	//Model Config
	config_t cfg; initDefaultConfig(cfg);
	cfg.model_depth = 6;
	cfg.spectrum_depths.push_back(2); cfg.spectrum_depths.push_back(4); cfg.spectrum_depths.push_back(6);
	cfg.spectrum_weights.push_back(1); cfg.spectrum_weights.push_back(1); cfg.spectrum_weights.push_back(1);
	cfg.lambda = 0.01;
	initDerivedConfig(cfg);

	//Create the molecule data
	NNParamTestMol moldata(&cfg);

	//Set some parameter weights
	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");	//5 features + bias =  6 features
	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 4; hlayer_numnodes[1] = 2; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear
	cfg.theta_nn_hlayer_num_nodes = hlayer_numnodes;
	cfg.theta_nn_layer_act_func_ids = act_ids;

	NNTestParam testparam(fnames, 1, hlayer_numnodes, act_ids);
	NNTestParam param(fnames, 1, hlayer_numnodes, act_ids);
	param.appendNextEnergyParams(testparam);
	param.appendNextEnergyParams(testparam);

	std::string param_filename = "tmp_param_file.log";
	if( boost::filesystem::exists( param_filename ) )
		boost::filesystem::remove( param_filename );
	param.saveToFile( param_filename );

	//Set some arbitrary suft values
	suft.values.resize(1);

	const FragmentGraph *fg = moldata.getFragmentGraph();
	unsigned int N = fg->getNumTransitions() + fg->getNumFragments();
	suft.values[0].assign(3*N, 0.0);
	suft.values[0][0] = 0.2; suft.values[0][1] = 0.3; suft.values[0][2] = 0.5; suft.values[0][3] = 0.0; suft.values[0][4] = 0.0;
	suft.values[0][N] = 0.9; suft.values[0][N+1] = 0.05; suft.values[0][N+2] = 0.05; suft.values[0][N+3] = 0.0; suft.values[0][N+4] = 0.0;
	suft.values[0][2*N] = 0.08; suft.values[0][2*N+1] = 0.9; suft.values[0][2*N+2] = 0.02; suft.values[0][2*N+3] = 0.0; suft.values[0][2*N+4] = 0.0;

	//Initialise all gradients to 0
	grads.resize(param.getNumWeights());
	for( unsigned int i = 0; i < grads.size(); i++ ) grads[i] = 0.0;

	std::set<unsigned int> used_idxs;
	FeatureCalculator fc_null(fnames); 
	std::string null_str = "null";
	EM_NN em(&cfg, &fc_null, null_str, param_filename );
	double Qonly  = em.computeQ(0, moldata, suft );
	double Q = em.computeAndAccumulateGradient(&grads[0], 0, moldata, suft, true, used_idxs);

	//Check Q
	double theta1 = 12.0, theta2 = 10.0;
	double rho_denom = 1.0 + exp(theta1) + exp(theta2);
	double x1 = exp(theta1)/rho_denom, x2 = exp(theta2)/rho_denom;	
	
	double expected_Q = (0.2+0.9+0.08)*(theta1-log(rho_denom)) + (0.3+0.05+0.9)*(theta2-log(rho_denom)) + (0.5+0.05+0.02)*(-log(rho_denom));
	if( fabs(Q - expected_Q )  > tol ){
		std::cout << "Unexpected Q value resulting from gradient computation: " << Q << std::endl;
		pass = false;
	}
	if( fabs(Qonly - expected_Q )  > tol ){
		std::cout << "Unexpected Q value resulting from Q only computation: " << Q << std::endl;
		pass = false;
	}

	//Check the gradients
	double exp_unweighted_1[] = {-5+5*x1,0,-5+5*x1,0,0,-5+5*x1,-4+4*x1,0,-4+4*x1,0,0,-4+4*x1,-8*x2, 0, 0, -8*x2,0,-8*x2,0,0,0,0,0,0,-2+2*x1+2*x2,-2*2+2*2*x1, 5*2-5*2*x1,2*x2,0,-1+x1+x2,-2+2*x1,5-5*x1,x2,0,1-x1-x2,3-3*x1-x2, -13+13*x1+7*x2 };
	double exp_unweighted_2[] = {5*x1,0,5*x1,0,0,5*x1, 4*x1,0,4*x1,0,0,4*x1, 8-8*x2, 0,0,8-8*x2,0,8-8*x2, 0,0,0,0,0,0, -2+2*x1+2*x2,4*x1,-5*2*x1,-2+2*x2,0,-1+x1+x2, 2*x1,-5*x1, -1+x2, 0, 1-x1-x2,+1-3*x1-x2,-7+13*x1+7*x2 };
	double exp_unweighted_persist[] = {5*x1,0,5*x1,0,0,5*x1, 4*x1,0,4*x1,0,0,4*x1,-8*x2, 0,0,-8*x2,0,-8*x2,0,0,0,0,0,0,2*x1+2*x2,4*x1,-5*2*x1,2*x2,0,x1+x2, 2*x1,-5*x1, x2, 0, -x1-x2,-3*x1-x2,13*x1+7*x2 };

	//Energy 0
	for( unsigned int i = 0; i < param.getNumWeightsPerEnergyLevel(); i++ ){
		if( fabs(grads[i] - (0.2*exp_unweighted_1[i] + 0.3*exp_unweighted_2[i] + 0.5*exp_unweighted_persist[i]))  > tol ){
			std::cout << "Unexpected gradient value at idx " << i;
			pass = false;
		}
	}
	//Energy 1
	int offset = param.getNumWeightsPerEnergyLevel();
	for( unsigned int i = 0; i < param.getNumWeightsPerEnergyLevel(); i++ ){
		if( fabs(grads[i+offset] - (0.9*exp_unweighted_1[i] + 0.05*exp_unweighted_2[i] + 0.05*exp_unweighted_persist[i]))  > tol ){
			std::cout << "Unexpected gradient value at idx " << i+offset;
			pass = false;
		}
	}
	//Energy 2
	offset = 2*param.getNumWeightsPerEnergyLevel();
	for( unsigned int i = 0; i < param.getNumWeightsPerEnergyLevel(); i++ ){
		if( fabs(grads[i+offset] - (0.08*exp_unweighted_1[i] + 0.9*exp_unweighted_2[i] + 0.02*exp_unweighted_persist[i]))  > tol ){
			std::cout << "Unexpected gradient value at idx " << i+offset;
			pass = false;
		}
	}


	//Check the used flags
	int expected_idxs[29] = {0,2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36};	
	if( used_idxs.size() != 29*3 ){
			std::cout << "Unexpected used index length " << used_idxs.size() << std::endl;
			pass = false;	
	}
	else{
		for( int energy = 0; energy <= 2; energy++ ){
			offset = energy*param.getNumWeightsPerEnergyLevel();
			for( unsigned int i = 0; i < 29; i++ ){
				if( used_idxs.find(expected_idxs[i]+offset) == used_idxs.end() ){
					std::cout << "Could not find expected used index " << expected_idxs[i]+offset << std::endl;
					pass = false;
				}
			}
		}
	}
	passed = pass;

}

NNParamsTestBiasIndexes::NNParamsTestBiasIndexes(){
	description = "Test return of bias indexes";
}

void NNParamsTestBiasIndexes::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	std::vector<std::string> fnames;
	fnames.push_back("IonicFeatures");	//5 features + bias =  6 features
	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 4; hlayer_numnodes[1] = 2; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear

	NNTestParam param(fnames, 1, hlayer_numnodes, act_ids);

	std::vector<unsigned int> bias_indexes;
	param.getBiasIndexes( bias_indexes );
	
	int expected_idxs[7] = {0,6,12,18,24,29,34};

	if( bias_indexes.size() != 7 ){
			std::cout << "Unexpected bias indexes length " << bias_indexes.size() << std::endl;
			pass = false;	
	}
	else{
		for( unsigned int i = 0; i < 7; i++ ){
			if( bias_indexes[i] != expected_idxs[i] ){
				std::cout << "Unexpected bias index " << bias_indexes[i] << std::endl;
				pass = false;
			}
		}
	}

	passed = pass;

}
