/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# comms_test.cpp
#
# Description: Test code for Comms.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "mpi.h"

#include "comms_test.h"

CommsTestSetMasterUsedIdxs::CommsTestSetMasterUsedIdxs(){
	description = "Test communication of used idxs";
}

void CommsTestSetMasterUsedIdxs::runTest(){
	
	bool pass = true;
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump < 3 ){
		std::cout << "Test should be run with at least 3 processors" << std::endl;
		passed = false;
		return;
	}

	Comms *comm;
	if( mpi_rank == MASTER ) comm = new MasterComms();
	else comm = new WorkerComms();

	//Leave the last processor with no used idxs
	if( mpi_rank < mpi_nump - 1 ){
		comm->used_idxs.insert( mpi_rank ); //Set a different used idx for each
		comm->used_idxs.insert( mpi_nump ); //Set a common used idx to all
	}

	//Communicate the indexes
	comm->setMasterUsedIdxs();

	//Check the resulting indexes in the master
	if( mpi_rank == MASTER ){
		
		MasterComms *mcomm = (MasterComms *)comm;

		for( int i = 1; i < mpi_nump - 1; i++ ){

			if( mcomm->worker_num_used[i] != 2 ){
				std::cout << "Unexpected number of used indexes for worker " << i << std::endl;
				pass = false;
			}

			if( mcomm->worker_used_idxs[i].find( mpi_nump ) == mcomm->worker_used_idxs[i].end() ){
				std::cout << "Missing index " << mpi_nump << " for worker used idxs " << i << std::endl;
				pass = false;			
			}

			if( mcomm->worker_used_idxs[i].find( i ) == mcomm->worker_used_idxs[i].end() ){
				std::cout << "Missing index " << i << " for worker used idxs " << i << std::endl;
				pass = false;			
			}
		}

		if( mcomm->worker_num_used[mpi_nump-1] != 0 ){
			std::cout << "Unexpected number of used indexes for worker " << mpi_nump-1 << std::endl;
			pass = false;
		}

		if( mcomm->master_used_idxs.size() != mpi_nump ){
			std::cout << "Unexpected number of master_used_idxs " << mcomm->master_used_idxs.size() << std::endl;
			pass = false;
		}

		for( int i = 0; i <= mpi_nump; i++ ){ 
			if( i != mpi_nump-1 ){
				if( mcomm->master_used_idxs.find(i) == mcomm->master_used_idxs.end() ){
					std::cout << "Couldn't find expected master used idx " << i << std::endl;
					pass = false;			
				}
			}
		}
	}

	delete comm;
	passed = pass;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait

}

CommsTestCollectQInMaster::CommsTestCollectQInMaster(){
	description = "Test communication and accumulation of Q";
}

void CommsTestCollectQInMaster::runTest(){
	
	bool pass = true;
	
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	Comms *comm;
	if( mpi_rank == MASTER ) comm = new MasterComms();
	else comm = new WorkerComms();

	double Q = (double)mpi_rank + 1.0;
	Q = comm->collectQInMaster(Q);

	double expected_Q = (double)(mpi_nump*(mpi_nump+1)/2);

	if( mpi_rank == MASTER && Q != expected_Q ){
		std::cout << "Q mismatch: " << Q << " vs " << expected_Q << std::endl;
		pass = false;		
	}
	delete comm;
	passed = pass;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait

}

CommsTestCollectGradsInMaster::CommsTestCollectGradsInMaster(){
	description = "Test communication and accumulation of gradients";
}

void CommsTestCollectGradsInMaster::runTest(){
	
	bool pass = true;
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump < 3 ){
		std::cout << "Test should be run with at least 3 processors" << std::endl;
		passed = false;
		return;
	}

	Comms *comm;
	if( mpi_rank == MASTER ) comm = new MasterComms();
	else comm = new WorkerComms();

	//Set used indexes as in Used Idx test (last has none)
	if( mpi_rank < mpi_nump - 1 ){
		comm->used_idxs.insert( mpi_rank ); //Set a different used idx for each
		comm->used_idxs.insert( mpi_nump ); //Set a common used idx to all
	}
	comm->setMasterUsedIdxs();

	//Set some gradients
	std::vector<double> grads( mpi_nump + 1 );
	for( unsigned int i = 0; i < grads.size(); i++ ){
		if( mpi_rank == MASTER ) grads[i] = 0.0;
		else grads[i] = rand();	//Make sure other unused grads are ignored
	}
	if( mpi_rank != mpi_nump - 1 ){
		grads[mpi_nump] = (double)mpi_rank + 1.0;
		grads[mpi_rank] = (double)mpi_rank + 1.0;
	}

	//Run the code
	comm->collectGradsInMaster(&grads[0]);

	//Check results
	if( mpi_rank == MASTER ){
		for( int i = 0; i < mpi_nump - 1; i++ ){
			if( grads[i] != (double(i) + 1.0) ){
				std::cout << "Grad mismatch at idx " << i << ": " << grads[i] << std::endl;
				pass = false;			
			}
		}
		if( grads[mpi_nump-1] != 0.0 ){
			std::cout << "Grad mismatch at idx " << mpi_nump-1 << ": " << grads[mpi_nump-1] << std::endl;
			pass = false;		
		}
		if( grads[mpi_nump] != (double)mpi_nump*(mpi_nump-1)/2 ){
			std::cout << "Grad mismatch at idx " << mpi_nump << ": " << grads[mpi_nump] << std::endl;
			pass = false;		
		}
	}

	delete comm;
	passed = pass;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait

}

CommsTestBroadcastParams::CommsTestBroadcastParams(){
	description = "Test broadcasting of parameters";
}

void CommsTestBroadcastParams::runTest(){
	
	bool pass = true;
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	Comms *comm;
	if( mpi_rank == MASTER ) comm = new MasterComms();
	else comm = new WorkerComms();

	std::vector<std::string> feature_list;
	feature_list.push_back( "BreakAtomPair" );
	Param param(feature_list, 3 );	//3 energy levels, 73 features.

	//Set used indexes
	unsigned int weights_per_energy = param.getNumWeightsPerEnergyLevel();
	comm->used_idxs.insert( mpi_rank ); //Set a different used idx for each
	comm->used_idxs.insert( mpi_rank + weights_per_energy ); 
	comm->used_idxs.insert( mpi_rank + 2*weights_per_energy); 
	comm->used_idxs.insert( mpi_nump ); //Set a common used idx to all
	comm->used_idxs.insert( mpi_nump + weights_per_energy); 
	comm->used_idxs.insert( mpi_nump + 2*weights_per_energy); 
	comm->setMasterUsedIdxs();

	//Set some parameters
	if( mpi_rank == MASTER ){
		for( unsigned int i = 0; i < weights_per_energy; i++){
			param.setWeightAtIdx((double)i, i);
			param.setWeightAtIdx((double)i*2, i + weights_per_energy);
			param.setWeightAtIdx((double)i*3, i + 2*weights_per_energy);
		}
	}

	//Run the code
	comm->broadcastParams(&param);

	//Check the results
	if( param.getWeightAtIdx(mpi_rank) != mpi_rank
		|| param.getWeightAtIdx(mpi_rank + weights_per_energy) != mpi_rank*2
		|| param.getWeightAtIdx(mpi_rank + 2*weights_per_energy) != mpi_rank*3 ){
			std::cout << "Rank parameters incorrectly set for processor " << mpi_rank << ": " << param.getWeightAtIdx(mpi_rank) << ": " << param.getWeightAtIdx(mpi_rank + weights_per_energy) << ": " << param.getWeightAtIdx(mpi_rank + 2*weights_per_energy) << std::endl;
		pass = false;				
	}

	if( param.getWeightAtIdx(mpi_nump) != mpi_nump
		|| param.getWeightAtIdx(mpi_nump + weights_per_energy) != mpi_nump*2
		|| param.getWeightAtIdx(mpi_nump + 2*weights_per_energy)  != mpi_nump*3 ){
		std::cout << "Common parameters incorrectly set for processor " << mpi_rank << std::endl;
		pass = false;				
	}

	delete comm;
	passed = pass;
	MPI_Barrier(MPI_COMM_WORLD);  	//All threads wait
}
