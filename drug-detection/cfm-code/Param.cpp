/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# param.cpp
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

#include "Param.h"
#include "Config.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

//Constructor to initialise parameter weight size from a feature list
Param::Param( std::vector<std::string> a_feature_list, int a_num_energy_levels ) :
	feature_list(a_feature_list), num_energy_levels(a_num_energy_levels) {
	
	FeatureCalculator fc( feature_list );
	unsigned int len = fc.getNumFeatures();
	unsigned int total_len = len * num_energy_levels;
	weights.resize( total_len );
	expected_num_input_features = len;
}

//Constructor to create skeleton parameter copy 
//(doesn't actually copy weights, just resizes to same)
Param::Param( Param &param_instance ){
	num_energy_levels = param_instance.num_energy_levels;
	weights.resize( param_instance.getNumWeights() );
}

//Append a set of parameters for the next energy level
void Param::appendNextEnergyParams( Param &next_param, int energy ){
	
	//Check that the features match
	if( feature_list.size() !=  next_param.getFeatureNames()->size() 
		|| getNumWeightsPerEnergyLevel() != next_param.getNumWeightsPerEnergyLevel() ){
		std::cout << "Mismatch in features for parameters to be appended" << std::endl;
		throw std::exception();
	}

	//Fetch the dimensions of the current weights
	int num_per_e_level = weights.size()/num_energy_levels;
	int start_offset = weights.size();
	
	//Append the new ones
	unsigned int new_start_offset = 0;
	unsigned int num_new_weights = next_param.getNumEnergyLevels()*num_per_e_level;
	if( energy < 0 ){ 
		num_energy_levels += next_param.getNumEnergyLevels();
	}else{ 
		new_start_offset = energy*num_per_e_level;
		num_new_weights = num_per_e_level;
		num_energy_levels++;
	}

	weights.resize( num_energy_levels*num_per_e_level );
	std::vector<double> *new_weights = next_param.getWeightsPtr();
	for( unsigned int i = 0; i < num_new_weights; i++ )
		weights[start_offset + i] = (*new_weights)[i + new_start_offset];
}

//Append a repeat of the highest energy's parameters (used to initialise high params with med etc).
void Param::appendRepeatedPrevEnergyParams(){
	
	//Fetch the dimensions of the current weights
	int num_per_e_level = weights.size()/num_energy_levels;
	
	//Append the repeated parameters
	unsigned int to_offset = weights.size();
	unsigned int from_offset = weights.size() - num_per_e_level;
	num_energy_levels++;

	weights.resize( num_energy_levels*num_per_e_level );
	for( unsigned int i = 0; i < num_per_e_level; i++ )
		weights[to_offset + i] = weights[from_offset + i];
}

//Randomly initialise all weights
void Param::randomInit(){

	// Non-Bias Terms: to uniform values between -0.5 and 0.5	
	for( unsigned int i = 1; i < weights.size(); i++ )
		weights[i] = (double(std::rand())/double(RAND_MAX) - 0.5);

	// Bias Terms: to uniform values between -3 and 0	
	unsigned int len = getNumWeightsPerEnergyLevel();
	for( unsigned int i = 0; i < num_energy_levels; i++ )
		weights[i*len] = (double(std::rand())/double(RAND_MAX) - 1.0)*3;

}

//Set all weights to zero except bias
void Param::zeroInit(){

	// Non-Bias Terms: All 0.0	
	for( unsigned int i = 1; i < weights.size(); i++ )
		weights[i] = 0.0;

	// Bias Terms: to uniform values between -3 and 0	
	unsigned int len = getNumWeightsPerEnergyLevel();
	double rand_bias = (double(std::rand())/double(RAND_MAX) - 1.0)*3;
	for( unsigned int i = 0; i < num_energy_levels; i++ )
		weights[i*len] = rand_bias;

}

//Set all weights to zero
void Param::fullZeroInit(){

	// All Terms 0.0	
	for( unsigned int i = 0; i < weights.size(); i++ )
		weights[i] = 0.0;

}

double Param::computeTheta( const FeatureVector &fv, int energy ){
	
	double theta = 0.0;
	
	//Check Feature Length
	int len = fv.getTotalLength();
	if( len != expected_num_input_features ){
		std::cout << "Expecting feature vector of length " << expected_num_input_features;
		std::cout << " but found " << len << std::endl;
		throw( ParamFeatureMismatchException() );
	}

	//Compute theta
	int energy_offset = len*energy;
	std::vector<feature_t>::const_iterator it = fv.getFeatureBegin();
	for( ; it != fv.getFeatureEnd(); ++it )
		theta += weights[*it + energy_offset];
	return theta;
}

void Param::adjustWeightsByGrads( std::vector<double> &grads, std::set<unsigned int> &used_idxs, double learning_rate, double momentum, std::vector<double> &prev_v ){

	std::set<unsigned int>::iterator it = used_idxs.begin();
	for( ; it != used_idxs.end(); ++it ){
		double v = momentum*prev_v[*it] + learning_rate*grads[*it];
		weights[*it] += v;
		prev_v[*it] = v;
	}
}


void Param::saveToFile( std::string &filename ){

	std::ofstream out;
	out.open(filename.c_str());
	if( !out.is_open() ){
		std::cout << "Warning: Trouble opening parameter file" << std::endl;
	}else{

		//Check the number of non-zero weights
		std::vector<double>::iterator itt = weights.begin();
		int num_used = 0;
		for( ; itt != weights.end(); ++itt )
			num_used += ( *itt != 0 );

		//Determine whether or not to use the sparse format
		bool use_sparse = false;
		if( (double)num_used/weights.size() < 0.25 ) use_sparse = true;

		//Use sparse format
		if(use_sparse) out << "SPARSE" << std::endl;

		//Print out the number of feature names
		out << feature_list.size() << std::endl;

		//Print out the feature names
		std::vector<std::string>::iterator it = feature_list.begin();
		for( ; it != feature_list.end(); ++it ) 
			out << *it << " ";
		out << std::endl;

		//Print out the number of energy levels
		out << num_energy_levels << std::endl;

		//Print out the total length of the weights
		out << weights.size() << std::endl;

		//Print out the number of used weights
		if(use_sparse) out << num_used << std::endl;

		out << std::setprecision(8);
		if(use_sparse){
			//Print out the used weights with indexes (in lines of 20)
			itt = weights.begin();
			for( int count = 0, idx = 0; itt != weights.end(); ++itt, idx++ ){ 
				if( *itt != 0 ){
					out << idx << " " << *itt;
					if( count % 20 == 19 ) out << std::endl;
					else out << " ";
					count++;
				}
			}
		}
		else{
			//Print out all the weights (in lines of 50) - no indexes
			itt = weights.begin();
			for( int idx = 0; itt != weights.end(); ++itt, idx++ ){ 
					out << *itt;
					if( idx % 50 == 49 ) out << std::endl;
					else out << " ";
			}		
		}
		out << std::endl;
		out.close();
	}
}

Param::Param( std::string &filename ){

	std::string line;
	std::ifstream ifs ( filename.c_str(), std::ifstream::in  );

	if(!ifs) std::cout << "Could not open file " << filename << std::endl;
	
	//Check for sparse representation
	getline( ifs, line );
	bool sparse = false;
	if( line[0] == 'S' ){ 
		sparse = true;
		getline( ifs, line );
	}

	//Get the number of feature names
	int num_features = atoi(line.c_str());

	//Get the feature names and compute the expected number of features
	std::string fname;
	getline( ifs, line );
	std::stringstream ss1(line);
	for( int i = 0; i < num_features; i++ ){
		ss1 >> fname;
		feature_list.push_back( fname );
	}
	FeatureCalculator fc( feature_list );
	expected_num_input_features = fc.getNumFeatures();

	//Get the number of energy levels
	getline( ifs, line );
	num_energy_levels = atoi(line.c_str());

	//Get the number of weights
	getline( ifs, line );
	int num_weights = atoi(line.c_str());
	weights.resize( num_weights, 0.0 );
	
	if(sparse){	//SPARSE FORMAT
		
		//Get the number of used weights
		getline( ifs, line );
		int num_used_weights = atoi(line.c_str());
		
		//Get the weights
		double weight;
		int count = 0, idx = -1;
		while( count < num_used_weights ){
			getline( ifs, line );
			std::stringstream ss2(line);
			int num_on_line = std::min( num_weights - count, 20);
			for( int i = 0; i < num_on_line; i++ ){
				ss2 >> idx >> weight;
				weights[idx] = weight;
				count++;
			}
		}
	}
	else{	//NON-SPARSE FORMAT
		
		//Get the weights
		double weight;
		int count = 0;
		while( count < num_weights ){
			getline( ifs, line );
			std::stringstream ss2(line);
			int num_on_line = std::min( num_weights - count, 50);
			for( int i = 0; i < num_on_line; i++ ){
				ss2 >> weight;
				weights[count++] = weight;
			}
		}
	}
	ifs.close();
}
