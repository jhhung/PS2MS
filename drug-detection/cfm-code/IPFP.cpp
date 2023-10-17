/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# ipfp.cpp
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

#include "IPFP.h"
#include "Config.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions.hpp>

//Big log negative value that will effectively result in zero probability
static const double NULL_PROB = -1000000000000.0;

//Maximum number of times to retry
static const int MAX_IPFP_OSC_RETRIES = 10;

double computeDifferenceRatio( double diff, double prev_diff ){
	double ratio = 0.0;
	if( diff == 0.0 && prev_diff == 0.0 ) ratio = 1.0;
	else if( prev_diff > 0 ) ratio = diff/prev_diff;
	return ratio;
}

IPFP::IPFP( MolData *data, config_t *config ){
	
	//Record pointers
	moldata = data;
	cfg = config;
	targets_modified = false;

	//Initialise difference buffer (for convergence check)
	unsigned int num_spec = config->spectrum_depths.size();
	diff_buf.resize( num_spec );
	prev_diff_buf.resize( num_spec );
	for( unsigned int i = 0; i < num_spec; i++){
		diff_buf[i] = 100.0;
		prev_diff_buf[i] = 1000.0;
	}

	//Initialise messages
	down_msgs.resize(cfg->model_depth - 1);
	up_msgs.resize(cfg->model_depth - 1);
	unsigned int num_fragments = moldata->getFragmentGraph()->getNumFragments();
	for( unsigned int i = 0; i < cfg->model_depth-1; i++ ){
		
		//Initialise all up messages to uniform (since there
		//is no evidence from below to indicate otherwise)
		up_msgs[i].reset(num_fragments);
		for( unsigned int j = 0; j < num_fragments; j++ )
			up_msgs[i].addToIdx(j, 0.0);
	}

}

void IPFP::initialiseFactorProbsAndBeliefs(){

	unsigned int num_fragments = moldata->getFragmentGraph()->getNumFragments();
	unsigned int num_transitions = moldata->getFragmentGraph()->getNumTransitions();

	//Note: copy them all across, regardless whether they are used
	fprobs.tn.resize(num_transitions);
	beliefs.tn.resize(num_transitions);
	for( unsigned int i = 0; i < num_transitions; i++ ){
		fprobs.tn[i].resize( cfg->model_depth );
		beliefs.tn[i].resize( cfg->model_depth );
		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			int energy = cfg->map_d_to_energy[d];
			fprobs.tn[i][d] = moldata->getLogTransitionProbForIdx(energy, i);
		}
	}
	fprobs.ps.resize(num_fragments);
	beliefs.ps.resize(num_fragments);
	for( unsigned int i = 0; i < num_fragments; i++ ){
		fprobs.ps[i].resize( cfg->model_depth );
		beliefs.ps[i].resize( cfg->model_depth );
		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			int energy = cfg->map_d_to_energy[d];
			fprobs.ps[i][d] = moldata->getLogPersistenceProbForIdx(energy, i);
		}	
	}
}

beliefs_t *IPFP::calculateBeliefs(){

	prev_marginals.resize(0);

	//Intialise all messages using simple junction tree inference
	Inference infer( moldata, cfg );
	infer.runInferenceDownwardPass( down_msgs, cfg->model_depth );

	//Initialise factor transition/persistence probabilities 
	//(these will change on evidence so we need this copy)
	initialiseFactorProbsAndBeliefs();

	//Apply Iterative Proportional Fitting Procedure or similar
	int num_iter = 0;
	if( cfg->ipfp_algorithm == IPFP_ALGORITHM ) 
		num_iter = runIPFP();
	else if( cfg->ipfp_algorithm == GEMA_ALGORITHM ) 
		num_iter = runGEMA();
	else if( cfg->ipfp_algorithm == IPFP_WITH_MOD_ALGORITHM ) 
		num_iter = runIPFPwithOscAdjust();
	else throw UnknownIPFPAlgorithmException();
	
	if( status == COMPLETE_CONVERGE ){ 
		// std::cout << "IPFP(GEMA) converged in " << num_iter << " iterations" << std::endl;
		return &beliefs;
	}

	if( status == CONVERGE_AFTER_MOD ){ 
		// std::cout << "IPFP converged after oscillatory modification in " << num_iter << " iterations" << std::endl;
		return &beliefs;
	}

	if( status == NON_CONVERGE ){
		std::cout << "Warning: IPFP did not converge: ";
		for( int i = 0; i < diff_buf.size(); i++ ) std::cout << diff_buf[i] << " ";
		return &beliefs;
	}
	std::cout << "Warning: Unknown return status from IPFP - " << status << std::endl;
	return &beliefs;
}


int IPFP::runIPFP( ){
	
	//Apply Iterative Proportional Fitting Procedure
	int iter = 0;
	targets_modified = false;
	double diff = computeBeliefs();
	std::vector<Message> desired_marginals, actual_marginals;
    computeMarginals( actual_marginals );
	setDesiredMarginals( desired_marginals, actual_marginals );
	while( !checkConvergence(false) && iter < MAX_IPFP_ITERATIONS ){
		iter++;
		for( unsigned int idx = 0; idx < cfg->dv_spectrum_depths.size(); idx++ )
			applyIPFPstep( actual_marginals, desired_marginals, idx );
	}
	return iter;
}

int IPFP::runIPFPwithOscAdjust( ){
	
	//Apply Iterative Proportional Fitting Procedure - detecting and adjusting for oscillatory convergence
	int iter = 0;
	int reverse = 0;
	targets_modified = false;

	//Store a copy of the initial msgs
	std::vector<Message> init_down_msgs = down_msgs;
	std::vector<Message> init_up_msgs = up_msgs;
	
	//Initialise the target marginals
	std::vector<Message> desired_marginals, actual_marginals, reverse_marginals;	
	initialiseFactorProbsAndBeliefs();
	double diff = computeBeliefs();
	computeMarginals( actual_marginals );
	setDesiredMarginals( desired_marginals, actual_marginals );	

	int retries = 0;
	while( (retries == 0 || status == OSC_CONVERGE) && retries <= MAX_IPFP_OSC_RETRIES ){

		//Run IPFP to convergence (oscillatory or otherwise) in the forward direction
		iter = 0;
		while( !checkConvergence(true) && iter < MAX_IPFP_ITERATIONS ){
			iter++;
			//for( int idx = cfg->dv_spectrum_depths.size()-1; idx >= 0; idx-- ){
			for( unsigned int idx = 0; idx < cfg->dv_spectrum_depths.size(); idx++ )
				applyIPFPstep( actual_marginals, desired_marginals, idx );
		}
		if( status != OSC_CONVERGE ) break;

		//Oscillatory Convergence detected:
		// std::cout << "Oscillatory convergence - modifying targets for retry after " << iter << " iterations ("  << retries << " retry)" << std::endl;

		//Re-Initialise msgs and factor transition/persistence probabilities
		down_msgs = init_down_msgs;
		up_msgs = init_up_msgs;
		initialiseFactorProbsAndBeliefs();
		diff = computeBeliefs();
		computeMarginals( reverse_marginals );
		for( unsigned int i = 0; i < diff_buf.size(); i++){
			diff_buf[i] = 100.0;
			prev_diff_buf[i] = 1000.0;
		}

		//Run IPFP to convergence (oscillatory or otherwise) in the reverse direction
		iter = 0;
		while( !checkConvergence(true) && iter < MAX_IPFP_ITERATIONS ){
			iter++;
			for( int idx = cfg->dv_spectrum_depths.size()-1; idx >= 0; idx-- )
				applyIPFPstep( reverse_marginals, desired_marginals, idx );
		}
		
		//Set the target marginals to the average of the forward and reverse marginals
		desired_marginals = actual_marginals;
		for( unsigned int j = 0; j < cfg->dv_spectrum_depths.size(); j++ ) 
			desired_marginals[j].addWeightedMessage(reverse_marginals[j], 1.0);	
		targets_modified = true;
		
		//Re-Initialise everything ready to go again
		down_msgs = init_down_msgs;
		up_msgs = init_up_msgs;
		initialiseFactorProbsAndBeliefs();
		diff = computeBeliefs();
		computeMarginals( actual_marginals );
		for( unsigned int i = 0; i < diff_buf.size(); i++){
			diff_buf[i] = 100.0;
			prev_diff_buf[i] = 1000.0;
		}

		retries++;

	}
	return iter;
}


int IPFP::runGEMA( ){
	
	//Apply GEMA extensions of Iterative Proportional Fitting Procedure
	int iter = 0;
	targets_modified = false;
	IPFP ipfp_copy(*this);
	double diff = computeBeliefs();
	std::vector<Message> desired_marginals, actual_marginals;
	std::vector<Message> gema_marginals(cfg->dv_spectrum_depths.size());
    computeMarginals( actual_marginals );
	setDesiredMarginals( desired_marginals, actual_marginals );
	while( !checkConvergence(false) && iter < MAX_IPFP_ITERATIONS ){
		iter++;

		//Re-initialise the GEMA target marginals	
		for( unsigned int idx = 0; idx < cfg->dv_spectrum_depths.size(); idx++ ) 
			gema_marginals[idx].reset(actual_marginals[idx].size());
		
		//Update the GEMA target marginals		
		for( unsigned int idx = 0; idx < cfg->dv_spectrum_depths.size(); idx++ ){
			ipfp_copy = *this;	//Create a copy of the IPFP instance so we don't change this one
			std::vector<Message> tmp_marginals( actual_marginals );
			ipfp_copy.applyIPFPstep( tmp_marginals, desired_marginals, idx );
			for( unsigned int j = 0; j < cfg->dv_spectrum_depths.size(); j++ ) 
				gema_marginals[j].addWeightedMessage(tmp_marginals[j], cfg->dv_spectrum_weights[idx]);
		}

		//Run standard IPFP for number of spectra iterations
		for( unsigned int idx = 0; idx < cfg->dv_spectrum_depths.size(); idx++ )
			applyIPFPstep( actual_marginals, gema_marginals, idx );

	}
	return iter;
}

void IPFP::applyIPFPstep( std::vector<Message> &actual_marginals, std::vector<Message> &desired_marginals, int idx ){
	
	applyEvidenceAndPropagate( cfg->dv_spectrum_depths[idx] - 1, desired_marginals[idx], actual_marginals[idx] );
	double diff = computeBeliefs();
	computeMarginals( actual_marginals );
	diff_buf[idx] = diff;
}

int IPFP::checkConvergence(bool allow_osc_converge){

	int prev_status = status;

	//Check for complete convergence - all diff values must be below tol
	if( targets_modified ) status = CONVERGE_AFTER_MOD;
	else status = COMPLETE_CONVERGE;

	std::vector<double>::iterator it = diff_buf.begin();
	for( ; it != diff_buf.end(); ++it ){
		if( *it >= cfg->ipfp_converge_thresh ){
			status = NON_CONVERGE;
			break;
		}
	}

	if( status == NON_CONVERGE && allow_osc_converge){
		//Check for oscillatory convergence - diff values must be cycling
		status = OSC_CONVERGE;
		for( unsigned int i = 0; i < diff_buf.size(); i++ ){
			if( prev_diff_buf[i] > 0.0 && (diff_buf[i] / prev_diff_buf[i]) <= cfg->osc_ipfp_converge_thresh ){
				status = NON_CONVERGE;
				break;
			} 
		}
	}
	
	//Copy across the diff_buf to the prev_diff_buf
	for( unsigned int i = 0; i < diff_buf.size(); i++ )
		prev_diff_buf[i] = diff_buf[i];

	return (status != NON_CONVERGE);
}

void IPFP::applyEvidenceAndPropagate( unsigned int spec_depth, Message &desired_marginals, Message &actual_marginals){

	const FragmentGraph *fg = moldata->getFragmentGraph(); 

	//Set the iterator for updating the factor
	Message *msg, tmp_msg(1);
	if( spec_depth > 0 ) msg = &(down_msgs[spec_depth-1]);
	else{	//Spectrum at depth 1 means there is no down message,
			//so only iterate over the 0 index
		tmp_msg.addToIdx(0, 0.0);
		msg = &tmp_msg;
	}
	
	//Update the affected factor	
	Message::const_iterator it = msg->begin();
	for( ; it != msg->end(); ++it ){
		
		unsigned int idx = it.index();

		//Persistence
		double tmp = beliefs.ps[idx][spec_depth];
		tmp += desired_marginals.getIdx(idx);
		tmp -= actual_marginals.getIdx(idx);
		tmp -= *it;	//down msg
		if( spec_depth < cfg->model_depth-1 )
			tmp -= up_msgs[spec_depth].getIdx(idx);
		fprobs.ps[idx][spec_depth] = tmp;

		//Transitions
		std::vector<int>::const_iterator itt = (*fg->getFromIdTMap())[idx].begin();
		for( ; itt != (*fg->getFromIdTMap())[idx].end(); ++itt ){

			const Transition *t = fg->getTransitionAtIdx(*itt);
			double tmp = beliefs.tn[*itt][spec_depth];
			tmp += desired_marginals.getIdx(t->getToId());
			tmp -= actual_marginals.getIdx(t->getToId());
			if( spec_depth > 0 )		
				tmp -= down_msgs[spec_depth-1].getIdx(t->getFromId());
			if( spec_depth < cfg->model_depth-1 )
				tmp -= up_msgs[spec_depth].getIdx(t->getToId());
			fprobs.tn[*itt][spec_depth] = tmp;
		}

	}

	//Propagate up to the top (update outgoing messages)
	for( int d = spec_depth-1; d >= 0; d-- ){
		it = down_msgs[d].begin();
		up_msgs[d].reset( fg->getNumFragments() );
		for( ; it != down_msgs[d].end(); ++it ){

			unsigned int idx = it.index();

			//Persistence
			if( (unsigned int)d < cfg->model_depth-2 ){
				double up_val = up_msgs[d+1].getIdx(idx);
				up_msgs[d].addToIdx(idx, fprobs.ps[idx][d+1] + up_val);
			}else up_msgs[d].addToIdx(idx, fprobs.ps[idx][d+1]);

			//Transitions
			std::vector<int>::const_iterator itt = (*fg->getFromIdTMap())[idx].begin();
			for( ; itt != (*fg->getFromIdTMap())[idx].end(); ++itt ){
				const Transition *t = fg->getTransitionAtIdx( *itt );
				double tmp;
				if( (unsigned int)d < cfg->model_depth-2 ){
					double up_val = up_msgs[d+1].getIdx(t->getToId());
					tmp = fprobs.tn[*itt][d+1] + up_val;
				}else tmp = fprobs.tn[*itt][d+1];
				up_msgs[d].addToIdx(idx, tmp);
			}
		}
	}

	//Propagate down to the bottom (update outgoing messages)
	for( unsigned int d = spec_depth; d < cfg->model_depth-1; d++ ){
		
		//Use a copy of the old down msg to get the used_idxs before resetting it
		Message old_down_msg = down_msgs[d];
		down_msgs[d].reset( fg->getNumFragments() );
		for( it = old_down_msg.begin(); it != old_down_msg.end(); ++it ){

			unsigned int idx = it.index();

			//Persistence
			if( d > 0 ){
				double down_val = down_msgs[d-1].getIdx(idx);
				down_msgs[d].addToIdx(idx, fprobs.ps[idx][d] + down_val);
			}else if( idx == 0 ) 
				down_msgs[d].addToIdx(idx, fprobs.ps[idx][d]);

			//Transitions
			std::vector<int>::const_iterator itt = (*fg->getToIdTMap())[idx].begin();
			for( ; itt != (*fg->getToIdTMap())[idx].end(); ++itt ){
				const Transition *t = fg->getTransitionAtIdx( *itt );
				if( d > 0 ){
					double down_val = down_msgs[d-1].getIdx(t->getFromId());
					down_msgs[d].addToIdx(idx, fprobs.tn[*itt][d] + down_val);
				}else if( t->getFromId() == 0 ) 
					down_msgs[d].addToIdx(idx, fprobs.tn[*itt][d]);
			} 
		}
	}
}

double IPFP::computeBeliefs( ){

	double diff = 0.0;
	beliefs_t tmp_beliefs;

	std::vector<double> norms;
	norms.resize(cfg->model_depth);

	const FragmentGraph *fg = moldata->getFragmentGraph();

	//Compute Persistence Beliefs (and track norms)
	tmp_beliefs.ps.resize(fg->getNumFragments());
	for( unsigned int i = 0; i < fg->getNumFragments(); i++ ){
		tmp_beliefs.ps[i].resize(cfg->model_depth);

		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			
			double tmp;
			if( (d == 0 && i == 0) || (d > 0 && down_msgs[d-1].getIdx(i) > -DBL_MAXIMUM) ){
				tmp = fprobs.ps[i][d];
				if( d < cfg->model_depth-1 ) 
					tmp += up_msgs[d].getIdx(i);
				if( d > 0 ){
					double msg_val = down_msgs[d-1].getIdx(i);
					if( msg_val > -DBL_MAXIMUM )
						tmp += msg_val;
					else tmp = NULL_PROB;
				}
				if ( i == 0 ) norms[d] = tmp;
				else norms[d] = logAdd(norms[d], tmp);
			}else tmp = NULL_PROB;
			tmp_beliefs.ps[i][d] = tmp;

		}	
	}

	//Compute Transition Beliefs (and track norms)
	tmp_beliefs.tn.resize(fg->getNumTransitions());
	for( unsigned int i = 0; i < fg->getNumTransitions(); i++ ){
		const Transition *t = fg->getTransitionAtIdx(i);
		tmp_beliefs.tn[i].resize(cfg->model_depth);

		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			
			double tmp;
			if( (d == 0 && t->getFromId() == 0) || (d > 0 && down_msgs[d-1].getIdx(t->getFromId()) > -DBL_MAXIMUM) ){
				tmp= fprobs.tn[i][d];
				if( d < cfg->model_depth-1 ) 
					tmp += up_msgs[d].getIdx(t->getToId());
				if( d > 0 ){
					double msg_val = down_msgs[d-1].getIdx(t->getFromId());
					if( msg_val > -DBL_MAXIMUM )
						tmp += msg_val;
					else tmp = NULL_PROB;
				}
				norms[d] = logAdd(norms[d], tmp);
			}else tmp = NULL_PROB;
			tmp_beliefs.tn[i][d] = tmp;
		}

	}

	//Normalise
	for( unsigned int i = 0; i < tmp_beliefs.tn.size(); i++ ){
		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			tmp_beliefs.tn[i][d] -= norms[d];
		}
	}
	for( unsigned int i = 0; i < tmp_beliefs.ps.size(); i++ ){
		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			tmp_beliefs.ps[i][d] -= norms[d];
		}
	}

	//Replace the old beliefs with the new ones, tracking the total difference
	diff = copyBeliefsAndComputeDiff( &tmp_beliefs );
	return diff;

}

double IPFP::copyBeliefsAndComputeDiff( beliefs_t *new_beliefs ){

	double diff = 0.0;
	for( unsigned int i = 0; i < new_beliefs->tn.size(); i++ ){
		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			diff += fabs( exp(beliefs.tn[i][d]) - exp(new_beliefs->tn[i][d]) );
			beliefs.tn[i][d] = new_beliefs->tn[i][d];
		}
	}
	for( unsigned int i = 0; i < new_beliefs->ps.size(); i++ ){
		for( unsigned int d = 0; d < cfg->model_depth; d++ ){
			diff += fabs( exp(beliefs.ps[i][d]) - exp(new_beliefs->ps[i][d]) );
			beliefs.ps[i][d] = new_beliefs->ps[i][d];
		}
	}
	return diff;
}

void IPFP::initLogSpecFactor( spec_factor_t &log_spec_factor, const Spectrum *spectrum ){

	//Store normpdf( pk mass, ion mass, sigma*sqrt2 )
	static const double pi = boost::math::constants::pi<double>();
	log_spec_factor.resize( spectrum->size() );
	const FragmentGraph *fg = moldata->getFragmentGraph();
	unsigned int num_fragments = fg->getNumFragments();
	Spectrum::const_iterator itp = spectrum->begin();
	for( unsigned int i = 0; itp != spectrum->end(); ++itp, i++ ){

		log_spec_factor[i].resize( num_fragments );
		
		//Set peak sigma based on the specified mass tolerances and the peak of interest
		//(assume tolerances cut things off at approx 3 x std deviation)
		double peak_sigma = 0.33*getMassTol(cfg->abs_mass_tol, cfg->ppm_mass_tol, itp->mass);
		
		double norm = -0.5*std::log(4*pi*peak_sigma*peak_sigma);
		double denom = 0.25/(peak_sigma*peak_sigma);
		for( unsigned int j = 0; j < num_fragments; j++ ){
			const Fragment *fgt = fg->getFragmentAtIdx(j);
			double tmp_mass_diff = fgt->getMass() - itp->mass;
			if( fabs(tmp_mass_diff) <= 3*peak_sigma ){
				if( cfg->obs_function == UNIFORM_OBS_FUNCTION ) log_spec_factor[i][j] = norm;
				else log_spec_factor[i][j] = norm - denom*tmp_mass_diff*tmp_mass_diff;
			}
			else
				log_spec_factor[i][j] = NULL_PROB;
		}

	}
}

void IPFP::initPeakTargets( Message &peak_targets, const Spectrum *spectrum ){

	peak_targets.reset(spectrum->size());
	Spectrum::const_iterator itp = spectrum->begin();
	for( unsigned int i = 0; itp != spectrum->end(); ++itp, i++ )
		peak_targets.addToIdx(i, std::log(itp->intensity * 0.01));
}

void IPFP::computeMarginal( Message &out, unsigned int spec_depth ){
	
	const FragmentGraph *fg = moldata->getFragmentGraph();
	out.reset(fg->getNumFragments());

	Message *msg, tmp_msg(1);
	if( spec_depth > 0 ) msg = &(down_msgs[spec_depth-1]);
	else{	//Spectrum at depth 1 means there is no down message,
			//so only iterate over the 0 index
		tmp_msg.addToIdx(0, 0.0);
		msg = &tmp_msg;
	}
	
	Message::const_iterator it = msg->begin();
	for( ; it != msg->end(); ++it ){
		unsigned int idx = it.index();

		//Persistence term
		out.addToIdx(idx, beliefs.ps[idx][spec_depth]);

		//Transition terms
		std::vector<int>::const_iterator itt = (*fg->getFromIdTMap())[idx].begin();
		for( ; itt != (*fg->getFromIdTMap())[idx].end(); ++itt ){
			const Transition *t = fg->getTransitionAtIdx( *itt );
			out.addToIdx( t->getToId(), beliefs.tn[*itt][spec_depth] );
		}
	}
}

void IPFP::computeMarginals(std::vector<Message> &marginals){
	marginals.resize( cfg->dv_spectrum_depths.size() );
	for( unsigned int i = 0; i < cfg->dv_spectrum_depths.size(); i++ )
		computeMarginal( marginals[i], cfg->dv_spectrum_depths[i]-1 );
}


void IPFP::setDesiredMarginals( std::vector<Message> &desired_margs, std::vector<Message> &marginals ){

	unsigned int num_spectra = cfg->dv_spectrum_depths.size();

	//Initialise factor probabilities between peaks and fragments
	//and peak targets
	std::vector<Message> peak_targets( num_spectra );
	std::vector<spec_factor_t> log_spec_factors( num_spectra );
	for( unsigned int energy = 0; energy < num_spectra; energy++ ){
		int spec_idx = cfg->dv_spectrum_indexes[energy];
		initLogSpecFactor( log_spec_factors[energy], moldata->getSpectrum(spec_idx) );
		initPeakTargets( peak_targets[energy], moldata->getSpectrum(spec_idx) );
	}

	//Set desired marginals
	Message pk_heights;
	desired_margs.resize(num_spectra);
	for( unsigned int energy = 0; energy < num_spectra; energy++ ){

		//Adjust the fragment marginals so that they sum to the correct target for each
		//peak yet maintain the correct ratio between fragments of similar mass
		//given prior likelihoods
		computePredictedPeakHeights( pk_heights, log_spec_factors[energy], marginals[energy] );

		//Use these to find the desired marginals for each fragment
		desired_margs[energy].reset(marginals[energy].size());
		Message::const_iterator it = marginals[energy].begin();
		for( ; it != marginals[energy].end(); ++it ){
			double log_sum = NULL_PROB;
			for( unsigned int i = 0; i < pk_heights.size(); i++ ){
				if( pk_heights.getIdx(i) > -10000.0 ){	//Ignore peaks that can't be realistically achieved 	
					double tmp = peak_targets[energy].getIdx(i);
					tmp -= pk_heights.getIdx(i);
					log_sum = logAdd( log_sum, tmp + log_spec_factors[energy][i][it.index()]);
				}
			}
			log_sum += *it;
			desired_margs[energy].addToIdx(it.index(), log_sum);
		}
	}
	
}

void IPFP::computePredictedPeakHeights( Message &out, spec_factor_t &log_spec_factor, Message &marginals){

	out.reset( log_spec_factor.size() );
	for( unsigned int i = 0; i < log_spec_factor.size(); i++ ){
		
		Message::const_iterator it = marginals.begin();
		double log_sum = NULL_PROB;
		for( ; it != marginals.end(); ++it ){
			double tmp = *it;
			tmp += log_spec_factor[i][it.index()];
			log_sum= logAdd( log_sum, tmp );
		}
		out.addToIdx(i,log_sum);
	}
}

void IPFP::printMarginals( std::vector<Message> &marginals ){

	std::cout << std::endl;
	for( unsigned int energy = 0; energy < cfg->dv_spectrum_depths.size(); energy++ ){
		Message::const_iterator it = marginals[energy].begin();
		for( ; it != marginals[energy].end(); ++it ){
			double prob = exp(*it);
			if( prob > 1e-2 ) 
				std::cout << it.index() << ":" << prob << " ";
		}
		std::cout << std::endl;
	}

}

IPFP &IPFP::operator=(const IPFP &rhs){

	//The moldata and config can point to the same structure without copying
	//since they are never modified anyway
	moldata = rhs.moldata;
	cfg = rhs.cfg;

	//Everything else will be vector-copied
	down_msgs = rhs.down_msgs;
	up_msgs = rhs.up_msgs;
	beliefs.tn = rhs.beliefs.tn;
	beliefs.ps = rhs.beliefs.ps;
	fprobs.tn = rhs.fprobs.tn;
	fprobs.ps = rhs.fprobs.ps;
	diff_buf = rhs.diff_buf;
	converge_thresh = rhs.converge_thresh;

	return *this;
}
