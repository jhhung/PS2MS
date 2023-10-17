/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# comms.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see param.cpp.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"
#include "Comms.h"

#include <iostream>
#include <sstream>
#include <exception>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/set.hpp>

Comms::Comms( ){
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );
}

void Comms::printWithWorkerId(const char *msg){
	std::cout << mpi_rank << ": " << msg << std::endl;
}

int Comms::collectSumInMaster( int partial){
	int tmp = partial;
	int total;
	MPI_Reduce( &tmp, &total, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
	return total;
}

void WorkerComms::setMasterUsedIdxs(){

	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	num_used = used_idxs.size();

	//Serialize the data for this worker's used idxs
	std::ostringstream str;
    boost::archive::text_oarchive ar(str);
	ar & used_idxs;
	std::string serialized_used_idxs = str.str();

	//Let the master know how many characters will be sent for this worker
	unsigned int data_size = serialized_used_idxs.size();
	MPI_Send( &data_size, 1, MPI_UNSIGNED, MASTER, 0, MPI_COMM_WORLD );

	//Send the master the (serialized) used idxs for this worker
	if( data_size > 0 ){
		//std::cout << mpi_rank << "start_sending" << "..." << std::endl;
		MPI_Send( &(serialized_used_idxs[0]), data_size, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD );
		//std::cout << mpi_rank << "end_sending" << "..." << std::endl;
	}
}

void MasterComms::printToMasterOnly(const char *msg){
	std::cout << msg << std::endl;
}

void MasterComms::setMasterUsedIdxs(){

	MPI_Status status;
	worker_used_idxs.resize( mpi_nump );
	worker_num_used.resize( mpi_nump );

	//Add the master's own used idxs first
	std::set<unsigned int>::iterator it = used_idxs.begin();
	for( ; it != used_idxs.end(); ++it )
		master_used_idxs.insert( *it );

	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait

	//Fetch each of the worker used idxs in turn
	for( int i = 1; i < mpi_nump; i++ ){
	
		//Find out the length of the incoming used idxs
		unsigned int data_size = 0;
		MPI_Recv( &data_size, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD, &status );
		if( data_size == 0 ){ 
			worker_num_used[i] = 0;
			worker_used_idxs[i].clear();
			continue;
		}
		
		//Receive the (serialized) used idxs for this worker
		std::string serialized_used_idxs;
		serialized_used_idxs.resize(data_size);
		//std::cout << i << "start_recv" << "..." << std::endl;
		MPI_Recv( &(serialized_used_idxs[0]), data_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status );
		//std::cout << i << "end_recv" << "..." << std::endl;

		//Unserialize the data and write to worker_used_idxs[i]
		std::istringstream ifs(serialized_used_idxs);
		boost::archive::text_iarchive ar(ifs);		
		ar & worker_used_idxs[i];
		worker_num_used[i] = worker_used_idxs[i].size();

		//Populate the master_used_idxs
		std::set<unsigned int>::iterator it = worker_used_idxs[i].begin();
		for( ; it != worker_used_idxs[i].end(); ++it )
			master_used_idxs.insert( *it );
		
	}

}

double Comms::collectQInMaster( double Q ){
	
	double Qsum;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	MPI_Reduce( &Q, &Qsum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	return Qsum;	//Note: Only the master has the real Qsum.
}

void Comms::broadcastInitialParams( Param *param ){
	std::vector<double> *weights = param->getWeightsPtr();
	MPI_Bcast(&((*weights)[0]), weights->size(), MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
}

void Comms::broadcastGorX( lbfgsfloatval_t *g, int n ){
	MPI_Bcast(g, n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
}

void WorkerComms::collectGradsInMaster( double *grads ){

	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	if( num_used == 0 ) return;

	std::vector<double> used_grads( num_used );
	std::set<unsigned int>::iterator it = used_idxs.begin();
	for( int i = 0; it != used_idxs.end(); ++it, i++ )
		used_grads[i] = *(grads + *it);
	MPI_Send( &(used_grads[0]), num_used, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD );
	

}

void MasterComms::collectGradsInMaster( double *grads ){

	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait

	//Receive and accumulate Gradients
	std::vector<double> used_grads;
	for( int i = 1; i < mpi_nump; i++ ){

		if( worker_num_used[i] > 0 ){
			used_grads.resize( worker_num_used[i] );
			MPI_Recv( &(used_grads[0]), worker_num_used[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );

			std::set<unsigned int>::iterator it = worker_used_idxs[i].begin();
			for( int j = 0; it != worker_used_idxs[i].end(); ++it, j++ )
				*(grads + *it) += used_grads[j];
		}
	}

}

void WorkerComms::broadcastParams( Param *param ){

	//Receive updated params from master
	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	if( num_used == 0 ) return;
	std::vector<double> used_params( num_used );
	MPI_Recv( &(used_params[0]), num_used, MPI_DOUBLE, MASTER, 0, MPI_COMM_WORLD, &status );
	
	//Update the params
	std::set<unsigned int>::iterator it = used_idxs.begin();
	for( int i = 0; it != used_idxs.end(); ++it, i++ )
		param->setWeightAtIdx( used_params[i], *it );
	
}

void MasterComms::broadcastParams( Param *param ){

	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait

	//Send each processor the changes of interest to them
	std::vector<double> used_params;
	for( int i = 0; i < mpi_nump; i++ ){
		
		if( worker_num_used[i] == 0 ) continue;
		used_params.resize( worker_num_used[i] );

		std::set<unsigned int>::iterator it = worker_used_idxs[i].begin();
		for( int j = 0; it != worker_used_idxs[i].end(); ++it, j++ )
			used_params[j] = param->getWeightAtIdx(*it);

		MPI_Send( &(used_params[0]), worker_num_used[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
		
	}
}

int Comms::broadcastConverged( int converged ){
	
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	MPI_Bcast(&converged, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	return converged;
	
} 

int Comms::broadcastNumUsed( int num_used ){
	
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	MPI_Bcast(&num_used, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	return num_used;
	
} 

double Comms::broadcastQ( double Q ){
	
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
	MPI_Bcast(&Q, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	return Q;
}
