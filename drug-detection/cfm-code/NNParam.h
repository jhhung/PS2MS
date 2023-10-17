/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param.h
#
# Description: 	Class for parameterization of fragmentation probabilities
#				within the bayesian network fragmentation trees.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __NN_PARAM_H__
#define __NN_PARAM_H__


#include "Param.h"

#include <string>

typedef std::vector<double> azd_vals_t;


//Exception to throw when the activation function id is unknown
class NNParamActivationFunctionIdException: public std::exception{

	virtual const char* what() const throw(){
		return "Unknown activation function";
	}
};

//Exception to throw when a param file is expected to contain the neural net configuration
//information but doesn't
class NNParamFileReadException: public std::exception{

	virtual const char* what() const throw(){
		return "Couldn't find neural net configuration information in param file";
	}
};


class NNParam : public Param {
public:
	NNParam( std::vector<std::string> a_feature_list, int a_num_energy_levels, std::vector<int> &a_hlayer_num_nodes, std::vector<int> &a_act_func_ids );
		
	//Constructor for loading parameters from file
	NNParam( std::string &filename );

	//Initialisation options
	void randomInit();
	
	//Save parameters to file (non-sparse format, includes neural net configuration)
	void saveToFile( std::string &filename );

	//Compute the theta value for an input feature vector and energy
	//based on the current weight settings
	double computeTheta( const FeatureVector &fv, int energy );
	//Same, but keeps a record of the z and a values along the way (for forwards step in forwards-backwards algorithm)
	double computeTheta( const FeatureVector &fv, int energy, azd_vals_t &z_values, azd_vals_t &a_values, bool already_sized = false );

	void computeDeltas( std::vector<azd_vals_t> &deltasA, std::vector<azd_vals_t> &deltasB, std::vector<azd_vals_t> &z_values, std::vector<azd_vals_t> &a_values, double rho_denom, int energy );
	void computeUnweightedGradients(  std::vector<std::vector<double> > &unweighted_grads, std::set<unsigned int> &used_idxs, std::vector<const FeatureVector *> &fvs, std::vector<azd_vals_t> &deltasA, std::vector<azd_vals_t> &deltasB, std::vector<azd_vals_t> &a_values );

	unsigned int getTotalNumNodes() const { return total_nodes; };	//Only hidden and output nodes (not input nodes)
	unsigned int getSecondLayerWeightOffset() const { return hlayer_num_nodes[0]*expected_num_input_features;};
	void getBiasIndexes( std::vector<unsigned int> &bias_indexes );

private:
	std::vector<int> hlayer_num_nodes;
	unsigned int total_nodes;

	//Activation functions and their derivatives (kept general in case we want to try other activation functions...)
	std::vector<int>act_func_ids;
	std::vector<double (*)(double)> act_funcs;
	std::vector<double (*)(double)> deriv_funcs;

	//Linear Activation (usually used for the last layer)
	static double linear_activation( double input ){ return input; };
	static double linear_derivative( double input ){ return 1; };

	//ReLU Activation
	static double relu_activation( double input ){ if( input > 0 ) return input; else return 0;};
	static double relu_derivative( double input ){ if( input > 0 ) return 1; else return 0;};
	static double neg_relu_activation( double input ){ if( input < 0 ) return input; else return 0;};	
	static double neg_relu_derivative( double input ){ if( input < 0 ) return 1; else return 0;};

	//Function to configure activation functions used in each layer
	void setActivationFunctionsFromIds();

	azd_vals_t tmp_z_values, tmp_a_values;

};


#endif // __NN_PARAM_H__
