/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.h
#
# Description: 	Class to apply Expectation Maximization algorithm to derive 
#				model parameters when using a neural net for theta.
#				- all identical to linear model, except params are NNParam and gradient 
#				  computation is different.
#					E-step: IPFP or equivalent.
#					M-step: Gradient Ascent
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __EM_NN_TRAIN_H__
#define __EM_NN_TRAIN_H__

#include "EM.h"
#include "NNParam.h"

class EM_NN : public EM {
public:
	//Constructor
	//Note: To include the group in the status filename, include _GRP_ in the name
	EM_NN(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename, std::string initial_params_filename = "");

	//This is public so the test can access it....there must be a better way?
	double computeAndAccumulateGradient(double *grads, int molidx, MolData &moldata, suft_counts_t &suft, bool record_used_idxs, std::set<unsigned int> &used_idxs );
	double computeQ(int molidx, MolData &moldata, suft_counts_t &suft);	
	double addRegularizers( double *grads );

	void writeParamsToFile( std::string &filename );

private:
	//The current parameters
	boost::shared_ptr<NNParam> nn_param; //EM::param will also be set to point to this location.
	void initParams();
	void computeThetas( MolData *moldata );

};

#endif // __EM_NN_TRAIN_H__
