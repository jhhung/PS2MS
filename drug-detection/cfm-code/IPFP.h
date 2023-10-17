/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# IPFP.h
#
# Description: 	Custom implementation of iterative proportional fitting
#				procedure. This is intended to replace the jtree inference
#				method for direct incorporation of marginal spectrum
#				evidence.
#				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __IPFP_H__
#define __IPFP_H__

#include "MolData.h"
#include "Message.h"
#include "Inference.h"

static const int MAX_IPFP_ITERATIONS = 2000;

//Status variables for the IPFP state of convergence
static const int NON_CONVERGE = 0;
static const int COMPLETE_CONVERGE = 1;
static const int OSC_CONVERGE = 2;
static const int CONVERGE_AFTER_MOD = 3;

//Algorithm identifiers
static const int IPFP_ALGORITHM = 0;
static const int GEMA_ALGORITHM = 1;
static const int IPFP_WITH_MOD_ALGORITHM = 2;

struct config_t;

//Exception to throw when the input feature configuration file is invalid 
class UnknownIPFPAlgorithmException: public std::exception{

	virtual const char* what() const throw(){
		return "Unknown algorithm id in IPFP belief calculation";
	}
};

class IPFP{

public:
	IPFP( MolData *data, config_t *config);
	IPFP( IPFP &ipfp_instance ){ *this = ipfp_instance; }
	beliefs_t *calculateBeliefs();
	int status;
	int algorithm_to_use;	//As defined above
	double converge_thresh;
	double osc_converge_thresh;
	IPFP &operator=(const IPFP &rhs);
	void applyIPFPstep( std::vector<Message> &actual_marginals, std::vector<Message> &desired_marginals, int idx );

private:
	MolData *moldata;
	config_t *cfg;
	
	beliefs_t beliefs;
	factor_probs_t fprobs;	
	
	std::vector<Message> down_msgs;
	std::vector<Message> up_msgs;
	
	std::vector<double> diff_buf;
	std::vector<double> prev_diff_buf;
	bool targets_modified;
	
	std::vector<Message> prev_marginals;
	
	int checkConvergence(bool allow_osc_converge);
	void applyEvidenceAndPropagate( unsigned int spec_depth, Message &desired_marginals, Message &actual_marginals);
	double computeBeliefs();
	void computeMarginals(std::vector<Message> &marginals);
	void computeMarginal( Message &out, unsigned int spec_depth );
	void initialiseFactorProbsAndBeliefs();
	double copyBeliefsAndComputeDiff( beliefs_t *new_beliefs );
	int runIPFP( );
	int runGEMA( );
	int runIPFPwithOscAdjust( );
	void recalculateTargets( int last_iter );
	void initLogSpecFactor( spec_factor_t &log_spec_factor, const Spectrum *spectrum );
	void setDesiredMarginals( std::vector<Message> &desired_margs, std::vector<Message> &marginals );
	void initPeakTargets( Message &peak_targets, const Spectrum *spectrum );
	void computePredictedPeakHeights( Message &out, spec_factor_t &log_spec_factor, Message &marginals);
	void printMarginals(std::vector<Message> &marginals);
};

#endif // __IPFP_H__
