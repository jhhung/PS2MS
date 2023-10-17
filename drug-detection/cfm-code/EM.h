/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.h
#
# Description: 	Class to apply Expectation Maximization algorithm to derive 
#				model parameters.
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

#ifndef __EM_TRAIN_H__
#define __EM_TRAIN_H__

static const int MAX_EM_ITERATIONS = 100;

#include "Config.h"
#include "MolData.h"
#include "Param.h"
#include "NNParam.h"
#include "Comms.h"
#include "IPFP.h"
#include <boost/shared_ptr.hpp>

#include <lbfgs.h>

//Sufficient stats, access by transition index for a given molecule
typedef std::vector<double> suft_t;

struct suft_counts_t{
	//Access each by molecule index
	std::vector<suft_t> values;
};

class EMComputationException: public std::exception{

	virtual const char* what() const throw(){
		return "Error during EM, unable to proceed.";
	}
};


class EM{
public:
	//Constructor
	//Note: To include the group in the status filename, include _GRP_ in the name
	EM(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename, std::string initial_params_filename = "");
	~EM();

	//Run the EM algorithm on the supplied data (except the specified group), 
	//return the final likelihood value.
	double run( std::vector<MolData> &data, int group, std::string &out_param_filename  );
	
	//After running EM, the final params can be written out to file
	virtual void writeParamsToFile( std::string &filename );

	//This is public so the test can access it....there must be a better way?
	virtual double computeAndAccumulateGradient(double *grads, int molidx, MolData &moldata, suft_counts_t &suft, bool record_used_idxs, std::set<unsigned int> &used_idxs);
	virtual double computeQ(int molidx, MolData &moldata, suft_counts_t &suft);

	//For use by LBFGS only (public for this reason).
	lbfgsfloatval_t evaluateLBFGS( const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);
	void progressLBFGS( const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);

	//Select mini batch (exposed publicly for testing...)
	void selectMiniBatch( std::vector<int> &initialized_minibatch_flags );

protected:
	//The feature calculator to use - preconfigured with feature spec
	FeatureCalculator *fc;

	//The current parameters
	boost::shared_ptr<Param> param;
	bool initial_params_provided;
	virtual void initParams();

	//Further virtual functions
	virtual void computeThetas( MolData *moldata );

	//Configuration data
	config_t *cfg;

	//Communicator (for exchanging data between threads)
	Comms *comm;
	int unused_zeroed;	//Use to note when unused parameters have been zeroed (so we don't
						//do it more times than we need to)
	void initComms();

	//For writing status messages to a log file
	std::string status_filename;
	void writeStatus( const char *msg );

	//Initialise sufficient statistics
	void initSuft( suft_counts_t &suft, std::vector<MolData> &data );

	//Update sufficient statistics based on beliefs
	void recordSufficientStatistics( suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs );

	//M-step: Run gradient ascent (using LBFGS package)
	double updateParametersLBFGS( std::vector<MolData> &data, suft_counts_t &suft );
	lbfgsfloatval_t *convertCurrentParamsToLBFGS(int N);
	void copyLBFGSToParams( const lbfgsfloatval_t *x );
	void copyGradsToLBFGS(lbfgsfloatval_t *g, std::vector<double> &grads, int n);

	//Simple gradient ascent
	double updateParametersSimpleGradientDescent( std::vector<MolData> &data, suft_counts_t &suft );

	//Helper functions
	virtual double addRegularizers( double *grads );
	void zeroUnusedParams();

	//tmp pointers for use during LBFGS
	suft_counts_t *tmp_suft_ptr_lbfgs;
	std::vector<MolData> *tmp_moldata_ptr_lbfgs;
	std::vector<int> tmp_minibatch_flags;

	int validation_group;
	bool sparse_params;
};

#endif // __EM_TRAIN_H__
