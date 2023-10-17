/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Learn the parameters of a model for the mass spec 
#				fragmentation process using EM, then use it to predict
#				the mass spectra.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "mpi.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <sstream>
#include <iostream>

#include "EM.h"
#include "EM_NN.h"
#include "Features.h"
#include "Config.h"
#include "MolData.h"
#include "FragmentGraphGenerator.h"

void parseInputFile(std::vector<MolData> &data, std::string &input_filename, int mpi_rank, int mpi_nump, config_t *cfg );
void trainCombinedEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data, int start_repeat);
void trainSingleEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data, int start_energy, int no_train, int start_repeat);

int main(int argc, char *argv[])
{
	int mpi_rank, mpi_nump;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if (argc < 4 || argc > 10)
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-train.exe <input_filename> <feature_filename> <config_filename> <peakfile_dir> <group> <status_filename> <no_train> <start_energy>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "input_filename:" << std::endl << "Text file with number of mols on first line, then " << std::endl << "id smiles_or_inchi cross_validation_group" << std::endl << "on each line after that." << std::endl;
		std::cout << std::endl << "feature_filename:" << std::endl << "Text file with list of feature names to include, line separated:" << std::endl << "BreakAtomPair" << std::endl << "IonRootPairs...etc" << std::endl;
		std::cout << std::endl<< "config_filename:" << std::endl << "Text file listing configuration parameters. Line separated 'name value'." << std::endl;
		std::cout << std::endl << "peakfile_dir_or_msp:" << std::endl << "Input MSP file, with ID fields corresponding to id fields in input_file (the MSP filename not including the .msp extension) OR Directory containing files with spectra. Each file should be called <id>.txt, where <id> is the id specified in the input file, and contains a list of peaks 'mass intensity' on each line, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. ";
		std::cout << "e.g." << std::endl << "energy0" << std::endl << "65.02 40.0" << std::endl << "86.11 60.0" << std::endl << "energy1" << std::endl << "65.02 100.0 ... etc" << std::endl;
		std::cout << std::endl << "group (opt):" << std::endl << "Cross validation group to run. Otherwise will assume 10 groups and run all of them." << std::endl;
		std::cout << std::endl << "status_filename (opt):" << std::endl << "Name of file to write logging information as the program runs. If not specified will write to status.log<group>, or status.log if no group is specified" << std::endl;
		std::cout << std::endl << "no_train (opt):" << std::endl << "Set to 1 if the training part should be skipped (useful in debugging - default 0)" << std::endl;
		std::cout << std::endl << "start_energy (opt - se only)"  << std::endl << "Set to starting energy if want to start training part way through (single energy only -default 0)" << std::endl;
		std::cout << std::endl << "start_repeat (opt)"  << std::endl << "Set to starting repeat if want to start training part way through (default 0)" << std::endl;
		std::cout << std::endl;
		exit(1);
	}

    std::string input_filename = argv[1];	//List (one per line): id, smiles_or_inchi, group
	std::string feature_filename = argv[2]; //List of features, line-spaced
	std::string config_filename = argv[3];	//Parameter configuration
	std::string peakfile_dir_or_msp = argv[4];	//MSP file or Directory containing the peak files for each molecule (in format <id>.txt)
	
	//Cross validation groups to process
	int min_group = 0, max_group = 9;
	if( argc >= 6 ){ 
		min_group = atoi(argv[5]);
		max_group = min_group;
	}
	std::string status_filename("status.log");	//Status file to write to
	if( argc >= 7 ) status_filename = argv[6];					
	else if( argc >= 6 ) status_filename = "status.log" + std::string(argv[5]);	//status.log<group>

	int no_train = 0;
	int start_energy = 0, start_repeat = 0;
	if( argc >= 8 ) no_train = atoi(argv[7]);
	if( argc >= 9 ) start_energy = atoi(argv[8]);
	if( argc >= 10 ) start_repeat = atoi(argv[9]);

	if( mpi_rank == MASTER ){
		//Create the tmp_data directory if it doesn't exist
		if( !boost::filesystem::exists("tmp_data") )
			boost::filesystem::create_directory("tmp_data");
		if( !boost::filesystem::exists("tmp_data/enumerated_output") )
			boost::filesystem::create_directory("tmp_data/enumerated_output");
		if( !boost::filesystem::exists("tmp_data/predicted_output") )
			boost::filesystem::create_directory("tmp_data/predicted_output");
		if( !boost::filesystem::exists("tmp_data/fv_fragment_graphs") )
			boost::filesystem::create_directory("tmp_data/fv_fragment_graphs");
		//Delete the status file if it already exists
		if( boost::filesystem::exists( status_filename ) )
			boost::filesystem::remove_all( status_filename );
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if( mpi_rank == MASTER ) std::cout << "Initialising Feature Calculator..";
	FeatureCalculator fc( feature_filename );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	if( mpi_rank == MASTER ) std::cout << "Initialising Parameter Configuration..";
	config_t cfg;
	initConfig( cfg, config_filename, mpi_rank == MASTER );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	if( mpi_rank == MASTER ) std::cout << "Parsing input file...";
	std::vector<MolData> data;
	parseInputFile( data, input_filename, mpi_rank, mpi_nump, &cfg );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	//Configure fragment graph state files
	std::string fv_filename_out = "tmp_data/fv_fragment_graphs/P" + boost::lexical_cast<std::string>(mpi_rank) + "_graphs.fg";
	std::string fv_filename_in = "tmp_data/fv_fragment_graphs/P" + boost::lexical_cast<std::string>(mpi_rank) + "_graphs.fg_tmp";
	if( min_group != 0 ) fv_filename_in = fv_filename_out;
	if( min_group == 0 && boost::filesystem::exists(fv_filename_in) ) 
		boost::filesystem::remove( fv_filename_in );
	if( min_group == 0 && boost::filesystem::exists(fv_filename_out) ) 
		boost::filesystem::copy_file( fv_filename_out, fv_filename_in );
	std::ifstream *fv_ifs; std::string next_id = "";
	if( boost::filesystem::exists(fv_filename_in) ){
		fv_ifs = new std::ifstream( fv_filename_in.c_str(), std::ifstream::in | std::ios::binary );
		if(!(*fv_ifs)) std::cout << "Could not open file " << fv_filename_in << std::endl;
		else{
			unsigned int id_size; fv_ifs->read(reinterpret_cast<char *>(&id_size), sizeof(id_size));
			next_id.resize(id_size); fv_ifs->read(&next_id[0], id_size);
		}
	}
	std::ofstream fv_out; 
	if( min_group == 0 ) fv_out.open(fv_filename_out.c_str(), std::ios::out | std::ios::binary );

	//Fragment Graph Computation (or load from file)
	time_t before_fg, after_fg;
	before_fg = time( NULL );
	if( mpi_rank == MASTER ) std::cout << "Computing fragmentation graphs and features..";
	std::vector<MolData>::iterator mit = data.begin();
	int success_count = 0, except_count = 0;
	for( ; mit != data.end(); ++mit ){
		try{
			//If we're not training, only load the ones we'll be testing
			if( (mit->getGroup() >= min_group && mit->getGroup() <= max_group) || !no_train ){
				if( mit->getId() == next_id ){
					mit->readInFVFragmentGraphFromStream( *fv_ifs );
					unsigned int id_size; fv_ifs->read(reinterpret_cast<char *>(&id_size), sizeof(id_size));
					if(!fv_ifs->eof()){ next_id.resize(id_size); fv_ifs->read(&next_id[0], id_size); }
					else next_id = "NULL_ID";
					
				}
				else{
					time_t before, after;
					before = time( NULL );
					mit->computeFragmentGraphAndReplaceMolsWithFVs(&fc);
					std::ofstream eout;
					eout.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
					after = time( NULL );
					eout << mit->getId() << "Done. Time Elaspsed = " << (after - before) << " Seconds";
					eout << " :Num Frag = " << mit->getFragmentGraph()->getNumFragments();
					eout << " :Num Trans = " << mit->getFragmentGraph()->getNumTransitions() << std::endl; 
					eout.close();
				}
				unsigned int id_size = mit->getId().size();
				if( min_group == 0 ){
					fv_out.write(reinterpret_cast<const char *>(&id_size), sizeof(id_size));
					fv_out.write(&(mit->getId()[0]), id_size);
					mit->writeFVFragmentGraphToStream( fv_out );	//We always write it, in case we haven't already computed all of them
				}
				success_count++;
			}
		}
		catch( std::exception e ){
			std::ofstream eout;
			eout.open(status_filename.c_str(), std::fstream::out | std::fstream::app);
			eout << "Exception occurred computing fragment graph for " << mit->getId() << std::endl;
			eout << mit->getSmilesOrInchi() << std::endl;
			eout << e.what() << std::endl << std::endl;
			except_count++;
			eout << except_count << " exceptions, from " << except_count + success_count << " total" << std::endl;
			eout.close();
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	after_fg = time( NULL );
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;
	std::cout << mpi_rank << ": " << success_count << " successfully computed. " << except_count << " exceptions." << std::endl;
	if( mpi_rank == MASTER ) std::cout << "Total Fragmentation Graph Computation Time Elaspsed = " << (after_fg - before_fg) << " Seconds";
	if( min_group == 0 ) fv_out.close(); 
	if( boost::filesystem::exists( fv_filename_in ) ){ 
		if( *fv_ifs ) fv_ifs->close(); 
		delete fv_ifs;
		if( min_group == 0 ) boost::filesystem::remove(fv_filename_in);
	}

	//Loading input spectra
	MPI_Barrier(MPI_COMM_WORLD);
	if( mpi_rank == MASTER ) std::cout << "Loading spectra..";
	bool spectra_in_msp = false; std::string pre_id = "";
	if( peakfile_dir_or_msp.substr(peakfile_dir_or_msp.size()-4, 4) == ".msp" ){ 
		spectra_in_msp = true;
		//Allow splitting of input data into multiple msps by naming them
		//some_file_name<P>.msp, where <P> gets replaced with the processor rank that uses that file
		pre_id = peakfile_dir_or_msp.substr(0, peakfile_dir_or_msp.size()-4);
		if( pre_id.substr(pre_id.size()-3, 3) == "<P>" ){
			peakfile_dir_or_msp.replace( peakfile_dir_or_msp.size()-7, 3, boost::lexical_cast<std::string>(mpi_rank) );
		}
	}
	
	//MSP Setup
	MspReader *msp; 
	std::ostream *out_enum_msp, *out_pred_msp; 
	std::ofstream of_emsp, of_pmsp;
	//Create the MSP lookup
	if( spectra_in_msp ) msp = new MspReader( peakfile_dir_or_msp.c_str(), "" );
	for( mit = data.begin(); mit != data.end(); ++mit ){ 
		if( (mit->getGroup() >= min_group && mit->getGroup() <= max_group) || !no_train ){
			if( spectra_in_msp )
				mit->readInSpectraFromMSP( *msp );
			else{
				std::string spec_file = peakfile_dir_or_msp + "/" + mit->getId() + ".txt";
				mit->readInSpectraFromFile( spec_file );
			}
			mit->removePeaksWithNoFragment( cfg.abs_mass_tol, cfg.ppm_mass_tol );
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);  
	if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

	//Training
	for( int group = min_group; group <= max_group; group++ ){
		if( mpi_rank == MASTER ) std::cout << "Running EM to train parameters for Group " << group << std::endl;	

		time_t before, after;
		before = time( NULL );
		std::string param_filename = "tmp_data/param_output";
		param_filename += boost::lexical_cast<std::string>(group);
		param_filename += ".log";

		if( cfg.use_single_energy_cfm )
			trainSingleEnergyCFM( param_filename, cfg, fc, status_filename, group, data, start_energy, no_train, start_repeat);
		else if(!no_train)
			trainCombinedEnergyCFM( param_filename, cfg, fc, status_filename, group, data, start_repeat);
		MPI_Barrier(MPI_COMM_WORLD);   	//Wait for all threads
		if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;
	
		if( mpi_rank == MASTER ){ 
			after = time( NULL );
			std::cout << "EM: Time Elaspsed = " << (after - before) << " Seconds";
		}

		//Open the output MSP files
		if( spectra_in_msp ){
			std::string proc_group_str = "P" +  boost::lexical_cast<std::string>(mpi_rank) + "G" +  boost::lexical_cast<std::string>(group);
			std::string enum_msp_filename = "tmp_data/enumerated_output/especs_" + proc_group_str + ".msp";
			std::string pred_msp_filename = "tmp_data/predicted_output/pspecs_" + proc_group_str + ".msp";
			of_emsp.open( enum_msp_filename.c_str()); of_pmsp.open( pred_msp_filename.c_str() );
			if( !of_emsp.is_open() || !of_pmsp.is_open() ){
				std::cout << "Warning: Trouble opening msp output files" << std::endl;
				exit(1);
			}
			std::streambuf *ebuf = of_emsp.rdbuf(), *pbuf = of_pmsp.rdbuf();	
			out_enum_msp = new std::ostream(ebuf);
			out_pred_msp = new std::ostream(pbuf);
		}

		if( mpi_rank == MASTER ) std::cout << "Generating Peak Predictions for Group " << group << "..." << std::endl;
		Param *param;
		if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
			param = new NNParam( param_filename );
		else
			param = new Param(param_filename);
		for( mit = data.begin(); mit != data.end(); ++mit ){ 
			if( mit->getGroup() != group ) continue;
			if( !mit->hasComputedGraph() ) continue;	//If we couldn't compute it's graph for some reason..

			//Predicted spectrum
			mit->computePredictedSpectra( *param, false );
			
			if( spectra_in_msp ) mit->writePredictedSpectraToMspFileStream( *out_pred_msp );
			else{
				std::string spectra_filename = "tmp_data/predicted_output/" + mit->getId() + ".txt";
				mit->writePredictedSpectraToFile(spectra_filename);
			}
		}
		delete param;
		if( mpi_rank == MASTER ) std::cout << "Done" << std::endl;

		if( spectra_in_msp ){ 
			of_emsp.close(); of_pmsp.close();
			delete out_pred_msp; delete out_enum_msp; 
		}

		MPI_Barrier(MPI_COMM_WORLD);   	//Wait for all threads
	}

	if( spectra_in_msp ) delete msp;


	MPI_Barrier(MPI_COMM_WORLD);   	//Wait for all threads
	MPI_Finalize();

	return(0);    
}

void trainCombinedEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data, int start_repeat){

	int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

	//Run EM multiple times with random restarts, taking the final one with the best Q
	double prev_Q = -1000000000;
	for( int repeat = start_repeat; repeat < cfg.num_em_restarts; repeat++ ){
		EM *em;
		std::string repeat_filename = param_filename + boost::lexical_cast<std::string>(repeat);
		std::string out_filename = repeat_filename;
		if( !boost::filesystem::exists( repeat_filename ) ) repeat_filename = "";
		if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
			em = new EM_NN( &cfg, &fc, status_filename, repeat_filename );
		else
			em = new EM( &cfg, &fc, status_filename, repeat_filename );

		double Q = em->run( data, group, out_filename );
		if( Q > prev_Q ){ 
			if( mpi_rank == MASTER ){ 
				std::cout << "Found better Q!" << std::endl;
				em->writeParamsToFile(param_filename);
			}
			prev_Q = Q;
		}
		delete em;
	}
}

void trainSingleEnergyCFM( std::string &param_filename, config_t &cfg, FeatureCalculator &fc, std::string &status_filename, int group, std::vector<MolData> &data, int start_energy, int no_train, int start_repeat){
	
	int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	Param *final_params;

	for( int energy = 0; energy < cfg.spectrum_depths.size(); energy++ ){

		std::string eparam_filename = param_filename + boost::lexical_cast<std::string>(energy) + "e";
		std::string prev_eparam_filename = param_filename + boost::lexical_cast<std::string>(energy-1) + "e";

		config_t se_cfg;
		initSingleEnergyConfig( se_cfg, cfg, energy );

		if( energy >= start_energy && !no_train ){
			//Run EM multiple times with random restarts, taking the final one with the best Q
			double prev_Q = -1000000000;
			if( energy > start_energy ) start_repeat = 0;
			if( start_repeat > 0 ){ //Too messy to compute previous best Q...just print a warning
				std::cout << "Warning: best Q for previous repeats unknown, may not pick up best params - check manually!" << std::endl;
			}
			for( int repeat = start_repeat; repeat < se_cfg.num_em_restarts; repeat++ ){

				EM *em;
				std::string repeat_filename = eparam_filename + boost::lexical_cast<std::string>(repeat);
				std::string out_filename = repeat_filename;
				if( !boost::filesystem::exists( repeat_filename ) ) repeat_filename = "";
				else if( energy > 0 && cfg.use_lower_energy_params_for_init ) 
					repeat_filename = prev_eparam_filename + boost::lexical_cast<std::string>(repeat);
				if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
					em = new EM_NN( &se_cfg, &fc, status_filename, repeat_filename );
				else
					em = new EM( &se_cfg, &fc, status_filename, repeat_filename );

				double Q = em->run( data, group, out_filename );
				if( Q > prev_Q ){ 
					if( mpi_rank == MASTER ){ 
						std::cout << "Found better Q!" << std::endl;
						em->writeParamsToFile(eparam_filename);
					}
					prev_Q = Q;
				}
				delete em;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);  
		if(energy == 0 ){ 
			if(cfg.theta_function == NEURAL_NET_THETA_FUNCTION)
				final_params = new NNParam(eparam_filename);
			else
				final_params = new Param(eparam_filename);
		}else{
			Param *eparam;
			if(cfg.theta_function == NEURAL_NET_THETA_FUNCTION) 
				eparam = new NNParam(eparam_filename);
			else
				eparam = new Param(eparam_filename);
			final_params->appendNextEnergyParams( *eparam, energy );
			delete eparam;
		}
	}
	if( mpi_rank == MASTER ) final_params->saveToFile(param_filename);
	delete final_params;
}



void parseInputFile(std::vector<MolData> &data, std::string &input_filename, int mpi_rank, int mpi_nump, config_t *cfg ){

	std::string line, smiles_or_inchi, id;
	std::ifstream ifs ( input_filename.c_str() , std::ifstream::in );
	int group, num_mols = 0;

	//Get the first line - the number of input molecules
	if( ifs.good() ){ 
		getline( ifs, line );
		num_mols = atoi(line.c_str());
	}
	else{
		std::cout << "Could not open input file " << input_filename << std::endl;
	}

	//Now get all the molecules
	int i = 0;
	while( ifs.good() && i < num_mols){
		i++;

		getline( ifs, line );
		if( line.size() < 3 ) continue;

		std::stringstream ss(line);
		ss >> id >> smiles_or_inchi >> group;

		//Split the data between processors. Only load in data for this
		//processor
		if( (i % mpi_nump) == mpi_rank )
			data.push_back( MolData( id, smiles_or_inchi, group, cfg) );
	}

}
