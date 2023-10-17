/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# EM.cpp
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

#include "mpi.h"
#include "EM.h"
#include "IPFP.h"
#include "Comms.h"
#include "Config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include "time.h"

EM::EM(config_t *a_cfg, FeatureCalculator *an_fc, std::string &a_status_filename, std::string initial_params_filename){
	cfg = a_cfg; 
	fc = an_fc;
	status_filename = a_status_filename;
	initComms();
	int num_energies_to_include = cfg->map_d_to_energy.back() + 1;
	if( initial_params_filename == "" ){
		param = boost::shared_ptr<Param>( new Param( fc->getFeatureNames(), num_energies_to_include ));
		initial_params_provided = false;
		comm->printToMasterOnly( "EM: No initial params provided" );
	}
	else{
		param = boost::shared_ptr<Param>( new Param( initial_params_filename ) );
		while( param->getNumEnergyLevels() < num_energies_to_include )
			param->appendRepeatedPrevEnergyParams();
		initial_params_provided = true; 
		std::string msg = "EM: Initial params provided from " + initial_params_filename;
		comm->printToMasterOnly( msg.c_str() );
	}
	sparse_params = true;
}

void EM::initComms(){
	//Initialise the communicator
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );
	if( mpi_rank == MASTER ) comm = new MasterComms();
	else comm = new WorkerComms();
}

EM::~EM(){
  delete comm;
}

void EM::writeStatus( const char *msg ){

	std::ofstream out;
	out.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
	out << msg << std::endl;
	out.close();
}

void EM::writeParamsToFile( std::string &filename ){
	param->saveToFile( filename );
}

void EM::initParams(){	
	if( cfg->em_init_type == PARAM_FULL_ZERO_INIT ) param->fullZeroInit();
	else if(cfg->em_init_type == PARAM_ZERO_INIT ) param->zeroInit();
	else param->randomInit();
}

void EM::computeThetas( MolData *moldata ){
	moldata->computeTransitionThetas( *param );
}

double EM::run( std::vector<MolData> &data, int group, std::string &out_param_filename ){

	unused_zeroed = 0;
	int iter = 0;

	if( !initial_params_provided ) initParams();
	comm->broadcastInitialParams( param.get() );
	validation_group = group;

	//Write the initialised params to file (we may get want to reload and use with saved suft state, even before updating)
	std::string init_out_param_filename = out_param_filename + "_init";
	if( comm->isMaster() ) writeParamsToFile( init_out_param_filename );

	// EM
	iter = 0;
	double Q, prevQ = -10000000.0;
	int count_no_progress = 0;
	while( iter < MAX_EM_ITERATIONS ){

		std::string iter_out_param_filename = out_param_filename + "_" + boost::lexical_cast<std::string>(iter);

		std::string msg = "EM Iteration " + boost::lexical_cast<std::string>(iter);
		if( comm->isMaster() ) writeStatus( msg.c_str() );
		comm->printToMasterOnly(msg.c_str());

		time_t before, after;
		std::vector<MolData>::iterator itdata;

		//Reset sufficient counts
		suft_counts_t suft;
		initSuft(suft, data);

		int num_converged = 0, num_nonconverged =0;
		int tot_numc = 0, total_numnonc = 0;
		before = time( NULL );

		//Do the inference part (E-step)
		itdata = data.begin();
		for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ ){
		
			if( !itdata->hasComputedGraph() ) continue;	//If we couldn't compute it's graph for some reason..
			if( itdata->hasEmptySpectrum() ){ 
				std::cout << "Warning: No peaks with explanatory fragment found for " << itdata->getId() << ", ignoring this input molecule." << std::endl;
				continue; //Ignore any molecule with poor (no peaks matched a fragment) or missing spectra.
			}

			if( itdata->getGroup() == validation_group ) continue;

			MolData *moldata = &(*itdata);

			//Compute the transition probabilities
			computeThetas(moldata);
			moldata->computeTransitionProbabilities();
			
			//Apply the peak evidencee, compute the beliefs and record the sufficient statistics
			if( cfg->use_single_energy_cfm ){
				beliefs_t beliefs;
				Inference infer( moldata, cfg );
				infer.calculateBeliefs( beliefs );
				recordSufficientStatistics( suft, molidx, moldata, &beliefs);
			}
			else{
				IPFP ipfp( moldata, cfg); 
				beliefs_t *beliefs = ipfp.calculateBeliefs();
				int status = ipfp.status;
				if( status == NON_CONVERGE || status == OSC_CONVERGE ) num_nonconverged++;
				else if( status == COMPLETE_CONVERGE || status == CONVERGE_AFTER_MOD ) num_converged++;
				recordSufficientStatistics( suft, molidx, moldata, beliefs);
			}

		}

		MPI_Barrier(MPI_COMM_WORLD);   	//All threads wait
		after = time( NULL );
		if( !cfg->use_single_energy_cfm ){
			total_numnonc = comm->collectSumInMaster( num_nonconverged );
			tot_numc = comm->collectSumInMaster( num_converged );
			std::string cvg_msg = "Num Converged: " + boost::lexical_cast<std::string>(tot_numc);
			std::string noncvg_msg = "Num Non-Converged: " + boost::lexical_cast<std::string>(total_numnonc);
			if( comm->isMaster() ){ writeStatus( cvg_msg.c_str() ); writeStatus( noncvg_msg.c_str() ); }
			comm->printToMasterOnly(cvg_msg.c_str()); comm->printToMasterOnly(noncvg_msg.c_str());
		}
		std::string estep_time_msg = "Completed E-step processing: Time Elapsed = " + boost::lexical_cast<std::string>(after - before) + " seconds";
		if( comm->isMaster() ) writeStatus(estep_time_msg.c_str());
		comm->printToMasterOnly(estep_time_msg.c_str());


		MPI_Barrier(MPI_COMM_WORLD);   	//All threads wait for master
		//Find a new set of parameters to maximize the expected log likelihood (M-step)
		before = time( NULL );
		if( cfg->use_lbfgs_for_ga ) 
			Q = updateParametersLBFGS(data, suft);
		else 
			Q = updateParametersSimpleGradientDescent(data, suft);
		after = time( NULL );
		std::string param_update_time_msg = "Completed M-step param update: Time Elapsed = " + boost::lexical_cast<std::string>(after - before) + " seconds";
		if( comm->isMaster() ) writeStatus(param_update_time_msg.c_str());
		comm->printToMasterOnly(param_update_time_msg.c_str());

		//Write the params
		if( comm->isMaster() ){ 
			writeParamsToFile( iter_out_param_filename );
			writeParamsToFile( out_param_filename );
		}
		MPI_Barrier(MPI_COMM_WORLD);   	//All threads wait for master
		after = time( NULL );
		std::string debug_time_msg = "Starting Q compute: Time Elapsed = " + boost::lexical_cast<std::string>(after - before) + " seconds";
		if( comm->isMaster() ) writeStatus(debug_time_msg.c_str());

		//Compute the final Q (with all molecules, in case only some were used in the mini-batch)
		if( cfg->ga_minibatch_nth_size > 1 ) Q = 0.0;
		double valQ = 0.0; int molidx = 0, numvalmols = 0, numnonvalmols = 0;
		for( itdata = data.begin(); itdata != data.end(); ++itdata, molidx++ ){
			if( itdata->getGroup() == validation_group ){ 
				//valQ += computeQ( molidx, *itdata, suft );
				numvalmols++;
			}
			else if(cfg->ga_minibatch_nth_size > 1 ){
				Q += computeQ( molidx, *itdata, suft );
				numnonvalmols++;
			}else	
				numnonvalmols++;
		}

		MPI_Barrier(MPI_COMM_WORLD);   	//All threads wait for master
		after = time( NULL );
		std::string debug_time_msg2 = "Finished Q compute: Time Elapsed = " + boost::lexical_cast<std::string>(after - before) + " seconds";
		if( comm->isMaster() ) writeStatus(debug_time_msg2.c_str());
		
		valQ = comm->collectQInMaster(valQ);
		if( cfg->ga_minibatch_nth_size > 1 ){
			Q = comm->collectQInMaster(Q);
			Q = comm->broadcastQ(Q);
		}
		numvalmols = comm->collectSumInMaster(numvalmols);
		numnonvalmols = comm->collectSumInMaster(numnonvalmols);

		//Check for convergence	
		double Qratio = fabs((Q-prevQ)/Q);
		if( comm->isMaster() ){
			std::string qdif_str = boost::lexical_cast<std::string>(Qratio) + " " + boost::lexical_cast<std::string>(prevQ) + " ";
			qdif_str += "Q=" + boost::lexical_cast<std::string>(Q) + " ValQ=" + boost::lexical_cast<std::string>(valQ) + " ";
			qdif_str += "Qn=" + boost::lexical_cast<std::string>(Q/numnonvalmols) + " ValQn=" + boost::lexical_cast<std::string>(valQ/numvalmols) + " ";
			writeStatus(qdif_str.c_str());
			comm->printToMasterOnly(qdif_str.c_str());
		}

		if( Qratio < 1e-15 ) count_no_progress += 1;
		else count_no_progress = 0;

		prevQ = Q;
		if( Qratio < cfg->em_converge_thresh || count_no_progress >= 3 ){ 
			comm->printToMasterOnly(("EM Converged after " + boost::lexical_cast<std::string>(iter) + " iterations").c_str());
			break;
		}
		iter++;
	}

	if( iter >= MAX_EM_ITERATIONS )
		comm->printToMasterOnly(("Warning: EM did not converge after " + boost::lexical_cast<std::string>(iter) + " iterations.").c_str());

	return Q;		
}

void EM::initSuft(suft_counts_t &suft, std::vector<MolData> &data ){

	//Resize the suft structure for each molecule
	unsigned int num_mols = data.size();
	suft.values.resize(num_mols);
	for( unsigned int i = 0; i < num_mols; i++){
		const FragmentGraph *fg = data[i].getFragmentGraph();
		int len =  fg->getNumTransitions() + fg->getNumFragments();
		int num_spectra = data[i].getNumSpectra();
		suft.values[i].resize( len*num_spectra );
	}
}

void EM::recordSufficientStatistics( suft_counts_t &suft, int molidx, MolData *moldata, beliefs_t *beliefs ){

	const FragmentGraph *fg = moldata->getFragmentGraph();

	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();

	int len_offset = num_transitions + num_fragments;

	//Accumulate the Sufficient Statistics
	for( unsigned int i = 0; i < num_transitions; i++ ){
		
		const Transition *t = fg->getTransitionAtIdx(i);

		double belief = 0.0;
		int energy = cfg->map_d_to_energy[0];
		if( t->getFromId() == 0 )	//main ion is always id = 0
			belief += exp(beliefs->tn[i][0]);
		for( unsigned int d = 1; d < cfg->model_depth; d++ ){
			energy = cfg->map_d_to_energy[d];
			if( energy != cfg->map_d_to_energy[d-1] ){
				suft.values[molidx][i + cfg->map_d_to_energy[d-1]*len_offset] = belief;
				belief = 0.0;
			}
			belief += exp(beliefs->tn[i][d]);
		}
		suft.values[molidx][i + energy*len_offset] = belief;
	}

	//Accumulate the persistence terms
	int offset = num_transitions;
	for( unsigned int i = 0; i < num_fragments; i++ ){
		
		double belief = 0.0;
		int energy = cfg->map_d_to_energy[0];
		if( i == 0 )	//main ion is always id = 0
			belief += exp(beliefs->ps[i][0]);
		for( unsigned int d = 1; d < cfg->model_depth; d++ ){
			energy = cfg->map_d_to_energy[d];
			if( energy != cfg->map_d_to_energy[d-1] ){
				suft.values[molidx][i + offset + cfg->map_d_to_energy[d-1]*len_offset] = belief;
				belief = 0.0;
			}
			belief += exp(beliefs->ps[i][d]);
		}
		suft.values[molidx][i + offset + energy*len_offset] = belief;
	}
}


static lbfgsfloatval_t lbfgs_evaluate( void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
{
	lbfgsfloatval_t fx = ((EM *)instance)->evaluateLBFGS(x, g, n, step);
    return fx;
}

static int lbfgs_progress( void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls)
{
	((EM *)instance)->progressLBFGS(x, g, fx, xnorm, gnorm, step, n, k,ls);
    return 0;
}


double EM::updateParametersLBFGS( std::vector<MolData> &data, suft_counts_t &suft ){

	int ret = 0;
	double Q = 0.0;
	lbfgs_parameter_t lparam;
    lbfgs_parameter_init(&lparam);
	lparam.delta = cfg->ga_converge_thresh;
	lparam.past = 1;
	lparam.max_iterations = cfg->ga_max_iterations;

	//Initial Q and gradient calculation (to determine used indexes - if we already have them don't bother)
	if( comm->used_idxs.size() == 0 ){
		std::vector<double> grads(param->getNumWeights(), 0.0);
		std::vector<MolData>::iterator itdata = data.begin();
		for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ ){
			if( itdata->getGroup() != validation_group )
				Q += computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft, true, comm->used_idxs);
		}

		//Collect the used_idxs from all the processors into the MASTER
		comm->setMasterUsedIdxs();

		//Copy the used parameters into the LBFGS array 
		if( comm->isMaster() ) zeroUnusedParams();
	}

	int N = 0;
	if( comm->isMaster() ) N = ((MasterComms *)comm)->master_used_idxs.size();
	N = comm->broadcastNumUsed( N );
	if( N > 0.1*param->getNumWeights() ) sparse_params = false;
	lbfgsfloatval_t *x = convertCurrentParamsToLBFGS(N);
	if( !sparse_params ) N = param->getNumWeights();

	//Select molecules to include in gradient mini-batch (will select a new set at each LBFGS iteration, but for every evaluation).
	tmp_minibatch_flags.resize(data.size());
	std::vector<MolData>::iterator itdata = data.begin();
	for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ )
		tmp_minibatch_flags[molidx] = ( itdata->getGroup() != validation_group ); //Don't include validation molecules
	if( cfg->ga_minibatch_nth_size > 1 ) selectMiniBatch(tmp_minibatch_flags);

	//Run LBFGS
	tmp_moldata_ptr_lbfgs = &data;
	tmp_suft_ptr_lbfgs = &suft;
	lbfgsfloatval_t fx;
	ret = lbfgs(N, x, &fx, lbfgs_evaluate, lbfgs_progress, this, &lparam);

	//Master converts and broadcasts final param weights and Q to all
	copyLBFGSToParams(x);
	if( sparse_params ) lbfgs_free(x);
	Q = -fx;
	comm->broadcastParams( param.get() );
	Q = comm->broadcastQ( Q );
	
	return( Q );

}

lbfgsfloatval_t *EM::convertCurrentParamsToLBFGS(int N){
    
	lbfgsfloatval_t *x;
	if( sparse_params){
		x = lbfgs_malloc(N);
		if (x == NULL) {
			std::cout << "ERROR: Failed to allocate a memory block for variables." << std::cout;
			throw EMComputationException();
		}

		//Only fill the actual parameters in the master (the others we can update at the start of
		//each evaluate call).
		if( comm->isMaster() ){
			std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
			for( unsigned int i = 0; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it )
				x[i++] = param->getWeightAtIdx(*it);	
		}
	}
	else x = &((*param->getWeightsPtr())[0]);	//If not sparse, just use the existing param array
	return x;
}

void EM::copyGradsToLBFGS(lbfgsfloatval_t *g, std::vector<double> &grads, int n){
    
	if( comm->isMaster() ){
		if( sparse_params ){
			std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
			for( unsigned int i = 0; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it )
				g[i++] = -grads[*it];	
		}
		else{
			for( int i = 0; i < n; i++ ) g[i] *= -1;
		}
	}
	comm->broadcastGorX( g, n );
}

void EM::copyLBFGSToParams( const lbfgsfloatval_t *x ){
	if( sparse_params && comm->isMaster() ){
		std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
		for( unsigned int i = 0; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it )
			param->setWeightAtIdx(x[i++], *it);	
	}
}

lbfgsfloatval_t EM::evaluateLBFGS( const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step){
	
	//Retrieve params from LBFGS and sync between processors
	copyLBFGSToParams(x);
	comm->broadcastParams( param.get() );

	//Set the location for the computed gradients (new array if sparse, else the provided g) and initialise
	double *grad_ptr;
	std::vector<double> grads;
	if( sparse_params ){ 
		grads.resize(param->getNumWeights(), 0.0);
		grad_ptr = &grads[0];
	}
	else{ 
		grad_ptr = g;
		std::set<unsigned int>::iterator sit = comm->used_idxs.begin();
		for( ; sit != comm->used_idxs.end(); ++sit ) *(grad_ptr + *sit) = 0.0;
	}

	//Compute Q and the gradient
	double Q = 0.0;
	std::vector<MolData>::iterator itdata = tmp_moldata_ptr_lbfgs->begin();
	for( int molidx = 0; itdata != tmp_moldata_ptr_lbfgs->end(); ++itdata, molidx++ ){
		if( tmp_minibatch_flags[molidx] )
			Q += computeAndAccumulateGradient(grad_ptr, molidx, *itdata, *tmp_suft_ptr_lbfgs, false, comm->used_idxs);
	}

	if( comm->isMaster() ) Q += addRegularizers( grad_ptr );
	comm->collectGradsInMaster( grad_ptr );
	Q = comm->collectQInMaster(Q);
	Q = comm->broadcastQ(Q);

	//Move the computed gradients into the lbfgs structure (note: only used idxs are included)
	copyGradsToLBFGS(g, grads, n);
	return -Q;
}

void EM::progressLBFGS( const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls){

	if( comm->isMaster() ){ 
		writeStatus(("LBFGS Iteration " + boost::lexical_cast<std::string>(k) + ": fx = " + boost::lexical_cast<std::string>(fx)).c_str());
		std::cout << "LBFGS Iteration " << k << ": fx = " << fx << std::endl;
	}

	//Select molecules to include in next gradient mini-batch (will select a new set at each LBFGS iteration, but not for every evaluation).
	std::vector<MolData>::iterator itdata = tmp_moldata_ptr_lbfgs->begin();
	for( int molidx = 0; itdata != tmp_moldata_ptr_lbfgs->end(); ++itdata, molidx++ )
		tmp_minibatch_flags[molidx] = ( itdata->getGroup() != validation_group ); //Don't include validation molecules
	if( cfg->ga_minibatch_nth_size > 1 ) selectMiniBatch(tmp_minibatch_flags);

}

double EM::updateParametersSimpleGradientDescent( std::vector<MolData> &data, suft_counts_t &suft ){
	
	double Q = 0.0, prev_Q = 100000.0;

	std::vector<double> grads(param->getNumWeights(), 0.0);
	std::vector<double> prev_v(param->getNumWeights(), 0.0);

	//Initial Q and gradient calculation (to determine used indexes)
	if( comm->used_idxs.size() == 0 ){
		std::vector<MolData>::iterator itdata = data.begin();
		for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ ){
			if( itdata->getGroup() != validation_group )
				Q += computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft, true, comm->used_idxs);
		}
		if( comm->isMaster() ) Q += addRegularizers( &grads[0] );	
		Q = comm->collectQInMaster(Q);
		Q = comm->broadcastQ(Q);
		comm->setMasterUsedIdxs();
		if( comm->isMaster() ) zeroUnusedParams();
	}
	int N = 0;
	if( comm->isMaster() ) N = ((MasterComms *)comm)->master_used_idxs.size();
	N = comm->broadcastNumUsed( N );

	int iter = 0;
	double learn_mult = 1.0;
	while( iter++ < cfg->ga_max_iterations && fabs((Q-prev_Q)/Q) >= cfg->ga_converge_thresh){

		if( Q < prev_Q && iter > 1) learn_mult = learn_mult*0.5;
		double learn_rate = cfg->starting_step_size*learn_mult/(1+ 0.01*iter);

		if(iter > 1) prev_Q = Q;

		//Select molecules to include in gradient mini-batch.
		std::vector<int> minibatch_flags(data.size());
		std::vector<MolData>::iterator itdata = data.begin();
		for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ )
			minibatch_flags[molidx] = ( itdata->getGroup() != validation_group ); //Don't include validation molecules
		if( cfg->ga_minibatch_nth_size > 1 ) selectMiniBatch(minibatch_flags);

		//Compute Q and the gradient
		std::vector<double>::iterator git =  grads.begin();
		for( ; git != grads.end(); ++git ) *git = 0.0;
		Q = 0.0; itdata = data.begin();
		for( int molidx = 0; itdata != data.end(); ++itdata, molidx++ ){
			if( minibatch_flags[molidx] ) 
				Q += computeAndAccumulateGradient(&grads[0], molidx, *itdata, suft, false, comm->used_idxs);
		}
		if( comm->isMaster() ) Q += addRegularizers( &grads[0] );
		comm->collectGradsInMaster( &grads[0] );
		Q = comm->collectQInMaster(Q);
		Q = comm->broadcastQ(Q);

		if( comm->isMaster() ) std::cout << iter << ": " << Q << " " << prev_Q << " " << learn_rate << std::endl;

		//Step the parameters
		if( comm->isMaster() )
			param->adjustWeightsByGrads( grads, ((MasterComms *)comm)->master_used_idxs, learn_rate, cfg->ga_momentum, prev_v );
		comm->broadcastParams( param.get() );

	}

	if( comm->isMaster() ){
		if( iter == cfg->ga_max_iterations )
			std::cout << "Gradient ascent did not converge" << std::endl;
		else
			std::cout << "Gradient ascent converged after " << iter << " iterations" << std::endl;
	}
	return Q;
}

double EM::computeAndAccumulateGradient(double *grads, int molidx, MolData &moldata, suft_counts_t &suft, bool record_used_idxs, std::set<unsigned int> &used_idxs){

	double Q = 0.0;
	const FragmentGraph *fg = moldata.getFragmentGraph();
	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();
	
	int offset = num_transitions;

	if( !moldata.hasComputedGraph() ) return Q;

	//Compute the latest transition thetas
	moldata.computeTransitionThetas( *param );
	suft_t *suft_values = &(suft.values[molidx]);
	
	//Collect energies to compute
	std::vector<unsigned int> energies;
	unsigned int energy;
	int prev_energy = -1;
	for( unsigned int d = 0; d < cfg->model_depth; d++ ){
		energy = cfg->map_d_to_energy[d];
		if( energy != prev_energy ) energies.push_back( energy );
		prev_energy = energy;
	}

	//Compute the gradients
	std::vector<unsigned int>::iterator eit = energies.begin();
	for( ; eit != energies.end(); ++eit ){
		energy = *eit;

		unsigned int grad_offset = energy * param->getNumWeightsPerEnergyLevel();
		unsigned int suft_offset = energy * (num_transitions + num_fragments);

		//Iterate over from_id (i)
		tmap_t::const_iterator it = fg->getFromIdTMap()->begin();
		for( int from_idx=0; it != fg->getFromIdTMap()->end(); ++it, from_idx++ ){

			//Calculate the denominator of the sum terms
			double denom = 1.0;
			std::vector<int>::const_iterator itt = it->begin();
			for( ; itt != it->end(); ++itt )
				denom += exp( moldata.getThetaForIdx(energy, *itt) );

			//Complete the innermost sum terms	(sum over j')	
			std::map<unsigned int, double> sum_terms;
			for( itt = it->begin(); itt != it->end(); ++itt ){
				const FeatureVector *fv = moldata.getFeatureVectorForIdx(*itt);
				std::vector<feature_t>::const_iterator fvit = fv->getFeatureBegin();				
				for( ; fvit != fv->getFeatureEnd(); ++fvit ){
					double val = exp( moldata.getThetaForIdx(energy, *itt))/denom;
					if( sum_terms.find(*fvit) != sum_terms.end() )
						sum_terms[*fvit] += val;
					else
						sum_terms[*fvit] = val;
				}
			}

			//Accumulate the transition (i \neq j) terms of the gradient (sum over j)
			double nu_sum = 0.0;
			for( itt = it->begin(); itt != it->end(); ++itt ){
				double nu = (*suft_values)[*itt + suft_offset];
				nu_sum += nu;
				const FeatureVector *fv = moldata.getFeatureVectorForIdx(*itt);
				std::vector<feature_t>::const_iterator fvit = fv->getFeatureBegin();				
				for( ; fvit != fv->getFeatureEnd(); ++fvit ){
					*( grads + *fvit + grad_offset ) += nu;
					if( record_used_idxs) used_idxs.insert(*fvit + grad_offset);
				}
				Q += nu*(moldata.getThetaForIdx(energy, *itt) - log(denom));
			}

			//Accumulate the last term of each transition and the 
			//persistence (i = j) terms of the gradient and Q
			std::map<unsigned int, double>::iterator sit = sum_terms.begin();
			double nu = (*suft_values)[offset + from_idx + suft_offset];	//persistence (i=j)
			for( ; sit != sum_terms.end(); ++sit ){ 
				*( grads + sit->first + grad_offset) -= (nu_sum + nu)*sit->second;
				if( record_used_idxs) used_idxs.insert(sit->first + grad_offset);
			}
			Q -= nu*log(denom);
		}

	}
	return Q;
}

double EM::computeQ(int molidx, MolData &moldata, suft_counts_t &suft){

	double Q = 0.0;
	const FragmentGraph *fg = moldata.getFragmentGraph();
	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();
	
	int offset = num_transitions;

	if( !moldata.hasComputedGraph() ) return Q;

	//Compute the latest transition thetas
	moldata.computeTransitionThetas( *param );
	suft_t *suft_values = &(suft.values[molidx]);
	
	//Collect energies to compute
	std::vector<unsigned int> energies;
	unsigned int energy;
	int prev_energy = -1;
	for( unsigned int d = 0; d < cfg->model_depth; d++ ){
		energy = cfg->map_d_to_energy[d];
		if( energy != prev_energy ) energies.push_back( energy );
		prev_energy = energy;
	}

	std::vector<unsigned int>::iterator eit = energies.begin();
	for( ; eit != energies.end(); ++eit ){
		energy = *eit;

		unsigned int suft_offset = energy * (num_transitions + num_fragments);

		//Iterate over from_id (i)
		tmap_t::const_iterator it = fg->getFromIdTMap()->begin();
		for( int from_idx=0; it != fg->getFromIdTMap()->end(); ++it, from_idx++ ){

			//Calculate the denominator of the sum terms
			double denom = 1.0;
			std::vector<int>::const_iterator itt = it->begin();
			for( ; itt != it->end(); ++itt )
				denom += exp( moldata.getThetaForIdx(energy, *itt) );

			//Accumulate the transition (i \neq j) terms of the gradient (sum over j)
			double nu_sum = 0.0;
			for( itt = it->begin(); itt != it->end(); ++itt ){
				double nu = (*suft_values)[*itt + suft_offset];
				Q += nu*(moldata.getThetaForIdx(energy, *itt) - log(denom));
			}

			//Accumulate the last term of each transition and the 
			//persistence (i = j) terms of the gradient and Q
			double nu = (*suft_values)[offset + from_idx + suft_offset];	//persistence (i=j)
			Q -= nu*log(denom);
		}

	}
	return Q;
}


double EM::addRegularizers( double *grads ){

	double Q = 0.0;
	std::set<unsigned int>::iterator it = ((MasterComms *)comm)->master_used_idxs.begin();
	for( ; it != ((MasterComms *)comm)->master_used_idxs.end(); ++it ){
		
		double weight = param->getWeightAtIdx(*it);
		Q -= 0.5*cfg->lambda*weight*weight;
		*(grads+*it) -= cfg->lambda*weight;		
	}

	//Remove the Bias terms (don't regularize the bias terms!)
	unsigned int weights_per_energy = param->getNumWeightsPerEnergyLevel();
	for( unsigned int energy = 0; energy < param->getNumEnergyLevels(); energy++ ){
		double bias = param->getWeightAtIdx(energy * weights_per_energy);
		Q += 0.5*cfg->lambda*bias*bias; 
		*(grads + energy * weights_per_energy) += cfg->lambda*bias;
	}
	return Q;
}

void EM::zeroUnusedParams(){

	unsigned int i;
	for( i = 0; i < param->getNumWeights(); i++ ){
		if( ((MasterComms *)comm)->master_used_idxs.find(i) == ((MasterComms *)comm)->master_used_idxs.end() )
			param->setWeightAtIdx(0.0, i);
	}
	
}


void EM::selectMiniBatch( std::vector<int> &initialized_minibatch_flags ){

	//The flags are initialized to 1's for selectable molecules, so set unwanted molecule flags to 0
	int num_mols = initialized_minibatch_flags.size();
	std::vector<int> idxs(num_mols);
	int count = 0;
	for( int i = 0; i < num_mols; i++ )
		if( initialized_minibatch_flags[i] ) idxs[count++] = i;
	idxs.resize(count);
	std::random_shuffle( idxs.begin(), idxs.end() );
	int num_minibatch_mols = (num_mols + cfg->ga_minibatch_nth_size - 1)/cfg->ga_minibatch_nth_size;
	for( int i = num_minibatch_mols; i < idxs.size(); i++ )
		initialized_minibatch_flags[idxs[i]] = 0;

}
