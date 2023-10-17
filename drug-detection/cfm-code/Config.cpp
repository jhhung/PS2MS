/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Config.cpp
#
# Description: 	Structs, functions, defaults for setting general 
#			    configuration data.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Config.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

void initDefaultConfig( config_t &cfg ){

	cfg.lambda = DEFAULT_LAMBDA;
	cfg.use_lbfgs_for_ga = DEFAULT_USE_LBFGS_FOR_GA;
	cfg.converge_count_thresh = DEFAULT_CONVERGE_COUNT_THRESH;
	cfg.em_converge_thresh = DEFAULT_EM_CONVERGE_THRESH;
	cfg.ga_converge_thresh = DEFAULT_GA_CONVERGE_THRESH;
	cfg.model_depth = DEFAULT_MODEL_DEPTH;
	cfg.abs_mass_tol = DEFAULT_ABS_MASS_TOL;
	cfg.ppm_mass_tol = DEFAULT_PPM_MASS_TOL;
	cfg.num_em_restarts = DEFAULT_NUM_EM_RESTARTS;
	cfg.line_search_alpha = DEFAULT_LINE_SEARCH_ALPHA;
	cfg.line_search_beta = DEFAULT_LINE_SEARCH_BETA;
	cfg.starting_step_size = 1.0;
	cfg.max_search_count = DEFAULT_MAX_SEARCH_COUNT;
	cfg.spectrum_depths.clear();
	cfg.spectrum_weights.clear();
	cfg.dv_spectrum_indexes.clear();
	cfg.intermediate_weights = 0.0;
	cfg.fg_depth = DEFAULT_FRAGGRAPH_DEPTH;
	cfg.ipfp_algorithm = DEFAULT_IPFP_ALGORITHM;
	cfg.ipfp_converge_thresh = DEFAULT_IPFP_CONVERGE_THRESH;
	cfg.osc_ipfp_converge_thresh = DEFAULT_IPFP_OSC_CONVERGE_THRESH;
	cfg.use_single_energy_cfm = 0;
	cfg.ionization_mode = DEFAULT_IONIZATION_MODE;
	cfg.update_bias_first = 0;
	cfg.em_init_type = PARAM_DEFAULT_INIT;
	cfg.use_lower_energy_params_for_init = 0;
	cfg.include_isotopes = DEFAULT_INCLUDE_ISOTOPES;
	cfg.isotope_thresh = DEFAULT_ISOTOPE_THRESH;
	cfg.allow_frag_detours = DEFAULT_ALLOW_FRAG_DETOURS;
	cfg.do_prelim_bfs = DEFAULT_DO_PRELIM_BFS;
	cfg.max_ring_breaks = DEFAULT_MAX_RING_BREAKS;
	cfg.theta_function = DEFAULT_THETA_FUNCTION;
	cfg.ga_minibatch_nth_size = DEFAULT_GA_MINIBATCH_NTH_SIZE;
	cfg.ga_max_iterations = DEFAULT_GA_MAX_ITERATIONS;
	cfg.ga_momentum = DEFAULT_GA_MOMENTUM;
	cfg.obs_function = DEFAULT_OBS_FUNCTION;
	cfg.include_h_losses = DEFAULT_INCLUDE_H_LOSSES;
	cfg.include_precursor_h_losses_only = DEFAULT_INCLUDE_PRECURSOR_H_LOSSES_ONLY;
	cfg.fragraph_compute_timeout_in_secs = DEFAULT_FRAGGRAPH_COMPUTE_TIMEOUT_IN_SECS;
}


void initConfig( config_t &cfg, std::string &filename, bool report_all ){

	std::string line, name;
	double value;
	std::ifstream ifs ( filename.c_str(), std::ifstream::in  );

	initDefaultConfig( cfg );

	//Read the config file into the paramater update config structure
	if(!ifs) std::cout << "Could not open file " << filename << std::endl;
	while( ifs.good() ){

		getline( ifs, line );
		if( line.size() < 3 ) continue;	//in case of empty line

		std::stringstream ss1(line);
		ss1 >> name >> value;

		if( name == "lambda" ) cfg.lambda = value;
		else if( name == "ionization_mode" ) cfg.ionization_mode = (int)value;
		else if( name == "converge_count_thresh" ) cfg.converge_count_thresh = (int)value;
		else if( name == "em_converge_thresh" ) cfg.em_converge_thresh = value;
		else if( name == "ga_converge_thresh" ) cfg.ga_converge_thresh = value;
		else if( name == "update_bias_first" ) cfg.update_bias_first = (int)value;
		else if( name == "model_depth" ) cfg.model_depth = (unsigned int)value;
		else if( name == "spectrum_depth" ) cfg.spectrum_depths.push_back( (unsigned int)value );
		else if( name == "spectrum_weight" ) cfg.spectrum_weights.push_back( (double)value );
		else if( name == "abs_mass_tol" ) cfg.abs_mass_tol = (double)value;
		else if( name == "ppm_mass_tol" ) cfg.ppm_mass_tol = (double)value;
		else if( name == "num_em_restarts" ) cfg.num_em_restarts = (int)value;
		else if( name == "line_search_alpha" ) cfg.line_search_alpha = (double)value;
		else if( name == "line_search_beta" ) cfg.line_search_beta = (double)value;
		else if( name == "starting_step_size" ) cfg.starting_step_size = (double)value;
		else if( name == "max_search_count" ) cfg.max_search_count = (int)value;
		else if( name == "fg_depth" ) cfg.fg_depth = (int)value;
		else if( name == "allow_frag_detours" ) cfg.allow_frag_detours = (int)value;
		else if( name == "do_prelim_bfs" ) cfg.do_prelim_bfs = (int)value;
		else if( name == "max_ring_breaks" ) cfg.max_ring_breaks = (int)value;
		else if( name == "ipfp_algorithm" ) cfg.ipfp_algorithm = (int)value;
		else if( name == "ipfp_converge_thresh" ) cfg.ipfp_converge_thresh = (double)value;
		else if( name == "osc_ipfp_converge_thresh" ) cfg.osc_ipfp_converge_thresh = (double)value;
		else if( name == "use_single_energy_cfm" ) cfg.use_single_energy_cfm = (int)value;
		else if( name == "include_isotopes" ) cfg.include_isotopes = (int)value;
		else if( name == "isotope_thresh" ) cfg.isotope_thresh = (double)value;
		else if( name == "use_lbfgs_for_ga" ) cfg.use_lbfgs_for_ga = (int)value;
		else if( name == "em_init_type" ) cfg.em_init_type = (int)value;
		else if( name == "use_lower_energy_params_for_init" ) cfg.use_lower_energy_params_for_init = (int)value;
		else if( name == "theta_function" ) cfg.theta_function = (int)value;
		else if( name == "theta_nn_hlayer_num_nodes" ) cfg.theta_nn_hlayer_num_nodes.push_back( (int)value );
		else if( name == "theta_nn_layer_act_func_ids" ) cfg.theta_nn_layer_act_func_ids.push_back( (int)value );
		else if( name == "ga_minibatch_nth_size") cfg.ga_minibatch_nth_size = (int)value;
		else if( name == "ga_max_iterations" ) cfg.ga_max_iterations = (int)value;
		else if( name == "ga_momentum" ) cfg.ga_momentum = (double)value;
		else if( name == "obs_function" ) cfg.obs_function = (int)value;
		else if( name == "include_h_losses" ) cfg.include_h_losses = (int)value;
		else if( name == "include_precursor_h_losses_only" ) cfg.include_precursor_h_losses_only = (int)value;
		else if( name == "fragraph_compute_timeout_in_secs" ) cfg.fragraph_compute_timeout_in_secs = (int)value;
		else std::cout << "Warning: Unknown paramater configuration identifier " << name << std::endl;
	}
	ifs.close();

	if( cfg.spectrum_depths.size() != cfg.spectrum_weights.size() )
		std::cout << "Warning: Mismatch between size of spectrum depths and weights" << std::endl;

	if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION ){
		if( cfg.theta_nn_layer_act_func_ids.size() < cfg.theta_nn_hlayer_num_nodes.size() + 1  ){
			std::cout << "Warning: Activations function types not specified for all neural net layers, using default activations for unspecified layers." << std::endl;
			while( cfg.theta_nn_layer_act_func_ids.size() < cfg.theta_nn_hlayer_num_nodes.size() + 1 )
				cfg.theta_nn_layer_act_func_ids.push_back( DEFAULT_NN_ACTIVATION_FUNCTION );
		}
		else if( cfg.theta_nn_layer_act_func_ids.size() != cfg.theta_nn_hlayer_num_nodes.size() + 1 ){
			std::cout << "Warning: More activations function types than neural net layers, ignoring some activations." << std::endl;
			cfg.theta_nn_layer_act_func_ids.resize( cfg.theta_nn_hlayer_num_nodes.size() + 1 );
		}
		for( int i = 0; i < cfg.theta_nn_hlayer_num_nodes.size(); i++ ){
			if( cfg.theta_nn_layer_act_func_ids[i] == RELU_NN_ACTIVATION_FUNCTION && 
				cfg.theta_nn_hlayer_num_nodes[i] % 2 > 0 ){
				std::cout << "Warning: Invalid to have odd number of nodes for ReLU hidden layer. Adding extra node." << std::endl;
				cfg.theta_nn_hlayer_num_nodes[i] += 1;
			}
		}
		cfg.theta_nn_hlayer_num_nodes.push_back(1);	//Last layer has one node
	}

	initDerivedConfig(cfg);

	//Report config parameters
	if( report_all ){
		if( cfg.use_single_energy_cfm ) std::cout << "Using Single Energy CFM" << std::endl;
		else std::cout << "Using Combined Energy CFM" << std::endl;
		if( cfg.ionization_mode == POSITIVE_ESI_IONIZATION_MODE ) std::cout << "Positive ESI Ionization Mode" << std::endl;
		else if( cfg.ionization_mode == NEGATIVE_ESI_IONIZATION_MODE ) std::cout << "Negative ESI Ionization Mode" << std::endl;
		else if( cfg.ionization_mode == POSITIVE_EI_IONIZATION_MODE ) std::cout << "Positive EI Ionization Mode" << std::endl;
		else{ 
			std::cout << "Warning: Unknown Ionization Mode, reverting to default mode (positive)!" << std::endl;
			cfg.ionization_mode = DEFAULT_IONIZATION_MODE;
		}
		if( cfg.include_isotopes ) std::cout << "Including fragment isotopes above intensity " << cfg.isotope_thresh << std::endl;
		else std::cout << "Not including fragment isotopes" << std::endl;
		if( cfg.include_h_losses || cfg.include_precursor_h_losses_only ){ 
			std::cout << "Including Hydrogen losses";
			if( cfg.include_precursor_h_losses_only ) std::cout << " from precursor only";
			std::cout << std::endl;
		}if( cfg.em_init_type == PARAM_RANDOM_INIT ) std::cout << "Using Random Parameter Initialisation" << std::endl;
		else if( cfg.em_init_type == PARAM_FULL_ZERO_INIT ) std::cout << "Using Full Zero Initialisation" << std::endl;
		else if( cfg.em_init_type == PARAM_ZERO_INIT ) std::cout << "Using Zero Initialisation (non-zero Bias)" << std::endl;
		else std::cout << "Warning: Unknown parameter initialization, revering to default mode (full random)!" << std::endl;
	
		std::cout << "Using EM Convergence Threshold " << cfg.em_converge_thresh << std::endl;
		std::cout << "Using Lambda " << cfg.lambda << std::endl;
		if( cfg.use_lbfgs_for_ga ) std::cout << "Using LBFGS package for gradient ascent" << std::endl;
		else{
			std::cout << "Using simple gradient ascent implementation" << std::endl;
			std::cout << "Using Starting Step Size " << cfg.starting_step_size << " and momentum " << cfg.ga_momentum << std::endl;
		}
		std::cout << "Using GA max iterations " << cfg.ga_max_iterations << std::endl;
		std::cout << "Using GA Convergence Threshold " << cfg.ga_converge_thresh << std::endl;
		std::cout << "Using GA mini batch taking 1 in " << cfg.ga_minibatch_nth_size << " of processor data" << std::endl;
		std::cout << "Using Fragmentation Graph Depth " << cfg.fg_depth << std::endl;
		if( cfg.allow_frag_detours ) std::cout << "Allowing fragmentation detours " << std::endl;
		else{ 
			std::cout << "Disallowing fragmentation detours ";
			if( cfg.do_prelim_bfs ) std::cout << "with preliminary breadth-first search" << std::endl;
			else std::cout << "without preliminary breadth-first search" << std::endl;
		}
		std::cout << "Maximum Ring Breaks " << cfg.max_ring_breaks << std::endl;
		std::cout << "Using Model Depth " << cfg.model_depth << std::endl;
		std::cout << "Using Spectrum Depths and Weights: ";
		for( unsigned int i =0; i < cfg.spectrum_depths.size(); i++ )
			std::cout << "(" << cfg.spectrum_depths[i] << "," << cfg.spectrum_weights[i] << ") ";
		std::cout << std::endl;
		std::cout << "Using Absolute mass tolerance " << cfg.abs_mass_tol << std::endl;
		std::cout << "Using PPM mass tolerance " << cfg.ppm_mass_tol << std::endl;
		if( !cfg.use_single_energy_cfm ){
			if( cfg.ipfp_algorithm == 0 ) std::cout << "Using standard IPFP" << std::endl;
			else if( cfg.ipfp_algorithm == 1 ) std::cout << "Using GEMA" << std::endl;
			else if( cfg.ipfp_algorithm == 2 ) std::cout << "Using IPFP with Oscillatory Adjustment" << std::endl;
			else std::cout << "Warning: Unknown IPFP algorithm id" << std::endl;
			std::cout << "Using IPFP Converge Thresh " << cfg.ipfp_converge_thresh << std::endl;
			std::cout << "Using IPFP Oscillatory Converge Thresh " << cfg.osc_ipfp_converge_thresh << std::endl;
		}
		if( cfg.use_lower_energy_params_for_init ) std::cout << "Initialising higher energy params with those of one level lower" << std::endl;
		if( cfg.theta_function == LINEAR_THETA_FUNCTION ) std::cout << "Using linear function for theta" << std::endl;
		else if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION ){ 
			std::cout << "Using neural net for theta with " << cfg.theta_nn_hlayer_num_nodes.size() << " hidden layers: ";
			for( int i = 0; i < cfg.theta_nn_hlayer_num_nodes.size(); i++ )
				std::cout << cfg.theta_nn_hlayer_num_nodes[i] << " ";
			std::cout << "and activation functions: ";
			for( int i = 0; i < cfg.theta_nn_layer_act_func_ids.size(); i++ )
				std::cout << cfg.theta_nn_layer_act_func_ids[i] << " ";
			std::cout << std::endl;
		}
		if( cfg.obs_function == NORMAL_OBS_FUNCTION ) std::cout << "Using normally distributed observation function" << std::endl;
		else if( cfg.obs_function == UNIFORM_OBS_FUNCTION ) std::cout << "Using windowed uniform distribution observation function" << std::endl;
		else{
			std::cout << "Warning: Unrecognised observation function (" << cfg.obs_function << "). Using default" << std::endl;
			cfg.obs_function = DEFAULT_OBS_FUNCTION;
		}
		if( cfg.fragraph_compute_timeout_in_secs > 0) std::cout << "Timeout set on fragment graph computation to " << cfg.fragraph_compute_timeout_in_secs << " mins" << std::endl;
	}
}

void initDerivedConfig( config_t &cfg, int se_energy ){

	//Force single energy where only one spectrum level
	if( cfg.spectrum_depths.size() == 1 )
		cfg.use_single_energy_cfm = 1;

	//Derived Parameters
	cfg.map_d_to_energy.resize( cfg.model_depth );
	int energy = 0;
	if( se_energy > 0 ) energy = se_energy;
	for( unsigned int d = 0; d < cfg.model_depth; d++ ){
		std::vector<int>::iterator it = cfg.spectrum_depths.begin();
		for( ; it != cfg.spectrum_depths.end(); ++it ){
			if( d == *it ) energy++;
		}
		cfg.map_d_to_energy[d] = energy;
	}

	cfg.dv_spectrum_depths = cfg.spectrum_depths;
	cfg.dv_spectrum_weights = cfg.spectrum_weights;

	//Re-normalise weights
	double sum = 0.0;
	std::vector<double>::iterator it = cfg.dv_spectrum_weights.begin();
	for( ; it != cfg.dv_spectrum_weights.end(); ++it ) sum += *it;
	it = cfg.dv_spectrum_weights.begin();
	for( ; it != cfg.dv_spectrum_weights.end(); ++it ) *it = *it/sum;

	//Set spectrum indexes
	cfg.dv_spectrum_indexes.clear();
	if( se_energy >= 0 ) cfg.dv_spectrum_indexes.push_back(se_energy);
	else{
		for( int i = 0; i < cfg.dv_spectrum_depths.size(); i++ )
			cfg.dv_spectrum_indexes.push_back(i);
	}

}

void initSingleEnergyConfig( config_t &se_cfg, config_t &cfg, int energy ){

	se_cfg = cfg;

	//Adjust the depth parameters to include only one spectrum
	se_cfg.model_depth = cfg.spectrum_depths[energy];	
	se_cfg.spectrum_depths.resize(1);
	se_cfg.spectrum_depths[0] = cfg.spectrum_depths[energy];
	se_cfg.spectrum_weights.resize(1);
	se_cfg.spectrum_weights[0] = 1.0;

	//Re-derive the derived parameters
	initDerivedConfig(se_cfg, energy);

}