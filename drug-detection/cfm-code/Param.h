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

#ifndef __PARAM_H__
#define __PARAM_H__

#include "Features.h"

#include <string>

//Exception to throw when the input feature vector configuration doesn't match the parameters
class ParamFeatureMismatchException: public std::exception{

	virtual const char* what() const throw(){
		return "Mismatch between feature vector length and num parameters";
	}
};


class Param{
public:
	//Constructor to initialise parameter weight size from a feature list
	Param( std::vector<std::string> a_feature_list, int a_num_energy_levels );

	//Constructor to create skeleton parameter copy 
	//(doesn't actually copy weights, just resizes to same)
	Param( Param &param_instance );

	//Constructor for loading parameters from file
	Param( std::string &filename );

	//Append a set of parameters to the current parameters
	//Either all energy levels to the next slot (if energy < 0), 
	//or just the energy level specified.
	void appendNextEnergyParams( Param &next_param, int energy = -1 );

	//Copy the highest energy parameters to create another higher energy level
	void appendRepeatedPrevEnergyParams();

	//Initialisation options
	virtual void randomInit();
	void zeroInit();
	void fullZeroInit();

	//Compute the theta value for an input feature vector and energy based
	//on the current weight settings
	virtual double computeTheta( const FeatureVector &fv, int energy );

	void adjustWeightsByGrads( std::vector<double> &grads, std::set<unsigned int> &used_idxs, double learning_rate, double momentum, std::vector<double> &prev_v );

	//Set the value of a weight
	void setWeightAtIdx( double value, int index ){ weights[index] = value; };

	//Print the current parameters to the given output stream
	void reportParameters( std::ostream &out );

	//Save parameters to file
	virtual void saveToFile( std::string &filename );

	//Access functions
	double getWeightAtIdx(int index){ return weights[index];};
	std::vector<double> *getWeightsPtr(){ return &weights; };
	unsigned int const getNumWeights(){ return weights.size(); };
	unsigned int const getNumWeightsPerEnergyLevel(){ return weights.size()/num_energy_levels; };
	unsigned int const getNumEnergyLevels(){ return num_energy_levels; };
	std::vector<std::string> *const getFeatureNames(){ return &feature_list; };

protected:
	std::vector<double> weights;
	unsigned int num_energy_levels;
	std::vector<std::string> feature_list;
	int expected_num_input_features;

};

#endif // __PARAM_H__
