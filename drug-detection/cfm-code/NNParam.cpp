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

#include "NNParam.h"
#include "Config.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

NNParam::NNParam( std::vector<std::string> a_feature_list, int a_num_energy_levels, std::vector<int> &a_hlayer_num_nodes, std::vector<int> &a_act_func_ids ) :
		hlayer_num_nodes( a_hlayer_num_nodes ),	act_func_ids(a_act_func_ids ), Param( a_feature_list, a_num_energy_levels ){

	//The weight length was set in the Param constructor, now alter it to match the neural net layout
	unsigned int num_features = weights.size() / num_energy_levels;
	int total_len = 0, num_input = num_features;
	std::vector<int>::iterator it =  hlayer_num_nodes.begin();
	total_len += (*it)*num_features; //The first layer takes the fv as input (which has in-built bias).
	total_nodes = *it;	
	num_input = *it; ++it;	
	for( ; it != hlayer_num_nodes.end(); ++it ){
		total_len += (*it)*(num_input+1);
		total_nodes += *it;
		num_input = *it;
	}
	weights.resize( total_len * num_energy_levels );

	//Set pointers to the activation functions and their derivatives used in each layer
	setActivationFunctionsFromIds();

	//Resize some temporary a and z vectors so we don't have to allocate memory for a and z within the computeTheta function if we don't need them after.
	tmp_z_values.resize(total_nodes);
	tmp_a_values.resize(total_nodes);
}

//Randomly initialise all weights
void NNParam::randomInit(){

	// All Terms: to uniform values between -0.1 and 0.1 - Biases too
	for( unsigned int i = 0; i < weights.size(); i++ )
		weights[i] = (double(std::rand())/double(RAND_MAX) - 0.5)*0.2;
}


double NNParam::computeTheta( const FeatureVector &fv, int energy ){
	return computeTheta( fv, energy, tmp_z_values, tmp_a_values, true );
}

double NNParam::computeTheta( const FeatureVector &fv, int energy, azd_vals_t &z_values, azd_vals_t &a_values, bool already_sized ){

	//Check Feature Length
	int len = fv.getTotalLength();
	if( len != expected_num_input_features ){
		std::cout << "Expecting feature vector of length " << expected_num_input_features;
		std::cout << " but found " << len << std::endl;
		throw( ParamFeatureMismatchException() );
	}
	int energy_offset = getNumWeightsPerEnergyLevel()*energy;
	std::vector<double>::iterator wit = weights.begin() + energy_offset;

	//Resize the z and a vectors to the required sizes
	std::vector<int>::iterator itlayer = hlayer_num_nodes.begin();
	if( !already_sized){ z_values.resize( total_nodes ); a_values.resize( total_nodes ); }
	azd_vals_t::iterator zit = z_values.begin();
	azd_vals_t::iterator ait = a_values.begin();
	std::vector<double (*)(double)>::iterator itaf = act_funcs.begin();

	//The first hidden layer takes the fv as input (which already has a bias feature, so no addtional biases in this layer)
	itlayer = hlayer_num_nodes.begin();
	for( int hnode = 0; hnode < (*itlayer); hnode++ ){
		double z_val = 0.0;
		std::vector<feature_t>::const_iterator it = fv.getFeatureBegin();
		for( ; it != fv.getFeatureEnd(); ++it )
			z_val += *(wit + *it);
		wit += len;
		*zit++ = z_val;
		*ait++ = (*itaf++)(z_val);
	}

	//Subsequent layers take the previous layer as input
	azd_vals_t::iterator ait_input_tmp, ait_input = a_values.begin();
	int num_input = *itlayer++;
	for( ; itlayer != hlayer_num_nodes.end(); ++itlayer ){
		for( int hnode = 0; hnode < (*itlayer); hnode++ ){
			ait_input_tmp = ait_input;
			double z_val = *wit++; //Bias
			for( int i = 0; i < num_input; i++ ) z_val += (*ait_input_tmp++)*(*wit++);
			*zit++ = z_val;
			*ait++ = (*itaf++)(z_val);			
		}
		num_input = *itlayer;
		ait_input = ait_input_tmp;
	}	
	return *(--ait);	//The output of the last layer is theta

}

void NNParam::saveToFile( std::string &filename ){

	//Write all the weights etc using the parent function
	Param::saveToFile( filename );	

	//Re-open the file and add the neural net configuration parameters
	std::ofstream out;
	out.open(filename.c_str(), std::fstream::app);
	if( !out.is_open() ){
		std::cout << "Warning: Trouble opening parameter file" << std::endl;
	}else{
		out << "Neural Net Config" << std::endl;
		out << hlayer_num_nodes.size() << std::endl;
		std::vector<int>::iterator iit = hlayer_num_nodes.begin();
		for( ; iit != hlayer_num_nodes.end(); ++iit )
			out << *iit << " ";
		out << std::endl;
		out << act_func_ids.size() << std::endl;
		iit = act_func_ids.begin();
		for( ; iit != act_func_ids.end(); ++iit )
			out << *iit << " ";
		out << std::endl;
		out.close();
	}
}

NNParam::NNParam( std::string &filename ) : Param(filename) {

	//We've already read all the weights, here we want to read the neural net parameters,
	//which should be appended right at the end of the file
	std::string line;
	std::ifstream ifs ( filename.c_str(), std::ifstream::in  );
	
	bool found_nn_details = false;
	if(!ifs) std::cout << "Could not open file " << filename << std::endl;
	while( getline( ifs, line ) ){
		if( line.size() > 6 && line.substr(0,6) == "Neural" ){
			//Get the hidden node number configuration
			getline( ifs, line );
			total_nodes = 0;
			int num_layers = atoi(line.c_str());
			hlayer_num_nodes.resize(num_layers);
			getline( ifs, line );
			std::stringstream ss1(line);
			for( int i = 0; i < num_layers; i++ ){
				ss1 >> hlayer_num_nodes[i];
				total_nodes += hlayer_num_nodes[i];
			}
			//Get the activation function configuration
			getline( ifs, line );
			int act_func_size = atoi(line.c_str());
			act_func_ids.resize(act_func_size);
			getline( ifs, line );
			std::stringstream ss2(line);
			for( int i = 0; i < num_layers; i++ )
				ss2 >> act_func_ids[i];
			setActivationFunctionsFromIds();
			tmp_z_values.resize(total_nodes);
			tmp_a_values.resize(total_nodes);
			found_nn_details = true;
			break;
		}
	}
	ifs.close();
	if( !found_nn_details ) throw NNParamFileReadException();
}

void NNParam::setActivationFunctionsFromIds(){

	std::vector<int>::iterator it, itt;
	itt = hlayer_num_nodes.begin();
	for( it = act_func_ids.begin(); it != act_func_ids.end(); ++it, ++itt ){
		for ( int hnode = 0; hnode < *itt; hnode++ ){	
			if( *it == LINEAR_NN_ACTIVATION_FUNCTION ){
				act_funcs.push_back( linear_activation );
				deriv_funcs.push_back( linear_derivative );		
			}else if( *it == RELU_NN_ACTIVATION_FUNCTION ){ 
				if( hnode % 2 == 1 ){
					act_funcs.push_back( neg_relu_activation );
					deriv_funcs.push_back( neg_relu_derivative );
				}
				else{
					act_funcs.push_back( relu_activation );
					deriv_funcs.push_back( relu_derivative );			
				}
			}else throw NNParamActivationFunctionIdException();
		}
	}

}

void NNParam::computeDeltas( std::vector<azd_vals_t> &deltasA, std::vector<azd_vals_t> &deltasB, std::vector<azd_vals_t> &z_values, std::vector<azd_vals_t> &a_values, double rho_denom, int energy ){

	//Resize the output delta vectors
	unsigned int num_trans_from_id = z_values.size();
	deltasA.resize( num_trans_from_id ); deltasB.resize( num_trans_from_id );
	for( int idx = 0; idx < num_trans_from_id; idx++ ){
		deltasA[idx].resize( total_nodes ); deltasB[idx].resize( total_nodes );
		for( int i = 0; i < total_nodes; i++ ) deltasA[idx][i] = i+1;	//DEBUG!
	}

	//Set up the iterators (we are going backwards through the network, so set up iterators in reverse)
	std::vector<azd_vals_t::reverse_iterator> itAs( num_trans_from_id ), itBs( num_trans_from_id );
	std::vector<azd_vals_t::reverse_iterator> aits( num_trans_from_id ), zits( num_trans_from_id );
	for( int idx = 0; idx < num_trans_from_id; idx++ ){
		aits[idx] = a_values[idx].rbegin(); zits[idx] = z_values[idx].rbegin();
		itAs[idx] = deltasA[idx].rbegin();  itBs[idx] = deltasB[idx].rbegin();
	}
	std::vector<double (*)(double)>::reverse_iterator itdf = deriv_funcs.rbegin();
	std::vector<int>::reverse_iterator itlayer = hlayer_num_nodes.rbegin();
	std::vector<azd_vals_t::reverse_iterator> itAs_input = itAs, itBs_input = itBs;

	int energy_offset = (num_energy_levels - energy - 1)*getNumWeightsPerEnergyLevel();
	std::vector<double>::reverse_iterator wit = weights.rbegin() + energy_offset;

	//Compute the highest level deltas
	for( int idx = num_trans_from_id-1; idx >= 0; idx-- ){
		double deriv_val = (*itdf)(*(zits[idx])++);
		double theta = (*(aits[idx])++);
		double rho_val = exp(theta)/rho_denom;
		*(itAs[idx])++ = deriv_val;
		*(itBs[idx])++ = deriv_val * rho_val;
	}
	itlayer++;
	itdf++;
	int num_output = 1;

	//Compute the deltas for the subsequent layers
	std::vector<double>::reverse_iterator itA_input_tmp, itB_input_tmp, wit_tmp;
	for( ; itlayer != hlayer_num_nodes.rend(); ++itlayer ){
		std::vector<double (*)(double)>::reverse_iterator itdf_tmp;
		for( int idx = num_trans_from_id-1; idx >= 0; idx-- ){
			itdf_tmp = itdf;
			for( int hnode_r = 0; hnode_r < *itlayer; hnode_r++ ){
				wit_tmp = wit + hnode_r;
				itA_input_tmp = itAs_input[idx]; itB_input_tmp = itBs_input[idx];			
				double Aval = 0.0, Bval = 0.0;
				double deriv_val = (*itdf_tmp)(*zits[idx]++);
				for( int oidx = num_output-1; oidx >= 0; oidx-- ){
					double w_deriv_val = (*wit_tmp)*deriv_val;
					Aval += (*itA_input_tmp++)*w_deriv_val;
					Bval += (*itB_input_tmp++)*w_deriv_val;
					wit_tmp += (oidx > 0)*(*itlayer+1);
				}
				*(itAs[idx])++ = Aval;
				*(itBs[idx])++ = Bval;
				itdf_tmp++;
			}
			itAs_input[idx] = itA_input_tmp; itBs_input[idx] = itB_input_tmp;
			wit_tmp++;
			wit_tmp++;	//Skip the bias
		}
		num_output = *itlayer;
		itdf = itdf_tmp;
		wit = wit_tmp;
	}

}

void NNParam::computeUnweightedGradients( std::vector<std::vector<double> > &unweighted_grads, 	std::set<unsigned int> &used_idxs, std::vector<const FeatureVector *> &fvs, std::vector<azd_vals_t> &deltasA, std::vector<azd_vals_t> &deltasB, std::vector<azd_vals_t> &a_values ){

	unsigned int num_trans_from_id = a_values.size();
	unweighted_grads.resize(num_trans_from_id + 1); //+1 for the persistence transition (stored last)
	std::vector<std::vector<double>::iterator> itgrads(num_trans_from_id), itAs(num_trans_from_id), itBs(num_trans_from_id), aits(num_trans_from_id);
	for( int idx = 0; idx < num_trans_from_id; idx++ ){
		unweighted_grads[idx].clear(); unweighted_grads[idx].resize( getNumWeightsPerEnergyLevel(), 0.0 );
		itgrads[idx] = unweighted_grads[idx].begin(); aits[idx] = a_values[idx].begin();
		itAs[idx] = deltasA[idx].begin(); itBs[idx] = deltasB[idx].begin();
	}
	//Persistence (and normalization) terms
	unweighted_grads[num_trans_from_id].clear(); unweighted_grads[num_trans_from_id].resize( getNumWeightsPerEnergyLevel(), 0.0 );
	
	//If there's no transitions from this id except persistence, the unweighted gradients for that persistence transition will all be 0
	if( num_trans_from_id == 0 ) return;
	
	std::vector<double>::iterator normit = unweighted_grads[num_trans_from_id].begin();
	unsigned int feature_len = fvs[0]->getTotalLength();


	//Collect the used feature idxs (we'll need these for the first layer normalization)
	std::set<unsigned int> tmp_used_idxs;
	for( int idx = 0; idx < num_trans_from_id; idx++ ){
		std::vector<feature_t>::const_iterator fit = fvs[idx]->getFeatureBegin();
		for( ; fit != fvs[idx]->getFeatureEnd(); ++fit ){
			tmp_used_idxs.insert( *fit );
			for( int hnode = 0; hnode < hlayer_num_nodes[0]; hnode++ )
				used_idxs.insert( hnode*feature_len + *fit );
		}
	}

	//First layer (only compute gradients for indices relevant to used features)	
	std::vector<int>::iterator itlayer = hlayer_num_nodes.begin();
	for( int hnode = 0; hnode < *itlayer; hnode++ ){

		//Compute the unnormalized terms, tracking the normalizers (persistence terms) as we go
		for( int idx = 0; idx < num_trans_from_id; idx++ ){
			std::vector<feature_t>::const_iterator fit = fvs[idx]->getFeatureBegin();
			double deltaA = *itAs[idx]++;
			double deltaB = *itBs[idx]++;
			for( ; fit != fvs[idx]->getFeatureEnd(); ++fit ){
				*(itgrads[idx] + *fit) = deltaA;
				*(normit + *fit) -= deltaB;
			}	
		}

		//Apply the normalizers 
		//(Note: we need to normalize for all transitions if any transition used that feature)
		for( int idx = 0; idx < num_trans_from_id; idx++ ){
			std::set<unsigned int>::iterator sit = tmp_used_idxs.begin();
			for( ; sit != tmp_used_idxs.end(); ++sit )
				*(itgrads[idx] + *sit) += *(normit + *sit);
			itgrads[idx] += feature_len;
		}
		normit += feature_len;
	}

	//Subsequent layers
	int num_input = *itlayer++;
	for( ; itlayer != hlayer_num_nodes.end(); ++itlayer ){
		for( int hnode = 0; hnode < *itlayer; hnode++ ){
			
			//Compute the unnormalized terms, tracking the normalizers (persistence terms) as we go
			for( int idx = 0; idx < num_trans_from_id; idx++ ){
				double deltaA = *itAs[idx]++;
				double deltaB = *itBs[idx]++;
				*itgrads[idx] = deltaA; //Bias term
				*normit -= deltaB;
				for( int input_node = 0; input_node < num_input; ++input_node ){
					double a_val = *(aits[idx] + input_node);
					*(itgrads[idx] + input_node + 1) = deltaA*a_val;
					*(normit + input_node + 1) -= deltaB*a_val;
				}	
			}

			//Apply the normalizers
			for( int idx = 0; idx < num_trans_from_id; idx++ ){
				*itgrads[idx] += *normit; //Bias term
				for( int input_node = 0; input_node < num_input; ++input_node )
					*(itgrads[idx] + input_node + 1) += *(normit + input_node + 1);
				itgrads[idx] += (num_input+1);
			}
			normit += (num_input+1);
		} 
		for( int idx = 0; idx < num_trans_from_id; idx++ ) aits[idx] += num_input;
		num_input = *itlayer;
	}
	
}


void NNParam::getBiasIndexes( std::vector<unsigned int> &bias_indexes ){

	bias_indexes.resize( total_nodes );
	int idx = 0, count = 0;
	std::vector<int>::iterator itlayer = hlayer_num_nodes.begin();
	for( int hnode = 0; hnode < *itlayer; hnode++, count++ ){
		bias_indexes[count] = idx;
		idx += expected_num_input_features;
	}
	int num_input = *itlayer++;
	for( ; itlayer != hlayer_num_nodes.end(); ++itlayer ){
		for( int hnode = 0; hnode < *itlayer; hnode++, count++ ){
			bias_indexes[count] = idx;
			idx += (num_input + 1);
		}
	}
}