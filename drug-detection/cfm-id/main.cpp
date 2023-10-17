/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Identify the most likely candidate structure for a given
#				spectrum from a list of candidates.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

int main(int argc, char *argv[]);

#include "FragmentGraphGenerator.h"
#include "MolData.h"
#include "Identifier.h"
#include "Comparators.h"

#include <iostream>
#include <fstream>
#include <string>

void readInCandidates( std::vector<Candidate> &candidates, std::string &candidate_file);
void reportResults( std::vector<Candidate> &candidates, std::ostream &out );

void printUsage(){

		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-id.exe <spectrum_file> <id> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol> <prob_thresh_for_prune> <param_filename> <config_filename> <score_type> <apply_postprocessing> <output_filename> <output_msp_or_mgf>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "spectrum_file:" << std::endl << "The filename where the input spectra can be found. This can be a .msp file in which the desired spectrum is listed under a corresponding id (next arg). Or it could be a single file with a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. ";
		std::cout << "e.g." << std::endl << "energy0" << std::endl << "65.02 40.0" << std::endl << "86.11 60.0" << std::endl << "energy1" << std::endl << "65.02 100.0 ... etc" << std::endl;
		std::cout << std::endl << "id:" << std::endl << "An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Otherwise not used but printed to output, in case of multiple concatenated results)" << std::endl;
		std::cout << std::endl << "candidate_file:" << std::endl << "The filename where the input list of candidate structures can be found - line separated 'id smiles_or_inchi' pairs." << std::endl;
		std::cout << std::endl << "num_highest (opt):" << std::endl << "The number of (ranked) candidates to return or -1 for all (if not given, returns all in ranked order)" << std::endl;
		std::cout << std::endl << "ppm_mass_tol (opt):" << std::endl << "The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm)" << std::endl;
		std::cout << std::endl << "abs_mass_tol (opt):" << std::endl << "The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs ( if not given defaults to 0.01Da)" << std::endl;
		std::cout << std::endl << "prob_thresh_for_prune (opt):" << std::endl << "The probabiltiy threshold at which to prune unlikely fragmnetations (default 0.001)" << std::endl;
		std::cout << std::endl << "param_filename (opt):" << std::endl << "The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory)" << std::endl;
		std::cout << std::endl << "config_filename (opt):" << std::endl << "The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory)" << std::endl;
		std::cout << std::endl << "score_type (opt):" << std::endl << "The type of scoring function to use when comparing spectra. Options: Jaccard (default for ESI-MS/MS), DotProduct (default for EI-MS)" << std::endl;
		std::cout << std::endl << "apply_postprocessing (opt):" << std::endl << "Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF (default for EI-MS), 1 = ON (default for ESI-MS/MS))." << std::endl; 
		std::cout << std::endl << "output_filename (opt):" << std::endl << "The filename of the output file to write to (if not given, prints to stdout)" << std::endl;
		std::cout << std::endl << "output_msp_or_mgf (opt):" << std::endl << "The filename for an output msp or mgf file to record predicted candidate spectra (if not given, doesn't save predicted spectra)" << std::endl;
};

int main(int argc, char *argv[])
{
	bool to_stdout = true;
	std::string output_filename;
	std::string param_filename = "param_output.log";
	std::string config_filename = "param_config.txt";
	int num_highest = -1;
	double abs_mass_tol = 0.01, ppm_mass_tol = 10.0;
	double prob_thresh_for_prune = 0.001;
	std::string score_type = "";
	int apply_postprocessing = -1;

	if (argc < 3 || argc > 14 ){
		printUsage();
		exit(1);
	}

	std::string spectrum_file = argv[1];
	std::string target_id = argv[2];
	std::string candidate_file = argv[3];
	if( argc > 4 ){ 
		try{ num_highest = boost::lexical_cast<int>(argv[4]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid num_highest: " << argv[4] << std::endl;
			exit(1);
		}
	}
	if( argc > 5 ){ 
		try{ ppm_mass_tol = boost::lexical_cast<float>(argv[5]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid ppm_tol: " << argv[5] << std::endl;
			exit(1);
		}
	}
	if( argc > 6 ){
		try{ abs_mass_tol = boost::lexical_cast<float>(argv[6]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid abs_tol: " << argv[6] << std::endl;
			exit(1);
		}
	}
	if( argc > 7 ){ 
		try{ prob_thresh_for_prune = boost::lexical_cast<float>(argv[7]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid prob_thresh_for_prune: " << argv[7] << std::endl;
			exit(1);
		}
	}
	if( argc > 8 ) param_filename = argv[8];
	if( argc > 9 ) config_filename = argv[9];
	if( argc > 10 ) score_type = argv[10];
	if( argc > 11 ){ 
		try{ apply_postprocessing = boost::lexical_cast<int>(argv[11]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid apply_postprocessing flag: " << argv[11] << std::endl;
			exit(1);
		}		
	}
	if( argc > 12 ){
		output_filename = argv[12];
		to_stdout = false;
	}

	//Set up for output to msp file
	std::string output_msp_or_mgf = "";
	if( argc > 13 ) output_msp_or_mgf = argv[13];

	//Initialise model configuration
	config_t cfg;
	initConfig( cfg, config_filename );

	//If apply postprocessing or score_type unspecified, set default according to ionization mode
	if( apply_postprocessing == -1 )
		apply_postprocessing = (int)(cfg.ionization_mode != POSITIVE_EI_IONIZATION_MODE );
	if( score_type.size() == 0 ){
		if( cfg.ionization_mode == POSITIVE_EI_IONIZATION_MODE )
			score_type = "DotProduct";
		else
			score_type = "Jaccard";
	}

	//Read in the input spectrum
	MolData targetData( target_id.c_str(), "Unknown", &cfg );
	if( spectrum_file.substr(spectrum_file.size()-4, 4) == ".msp" ){
		MspReader msp( spectrum_file.c_str(), "" );
		targetData.readInSpectraFromMSP( msp );
	}else targetData.readInSpectraFromFile( spectrum_file );

	//Fetch the list of candidates
	std::vector<Candidate> candidates;
	readInCandidates( candidates, candidate_file );

	//Set up the comparator and identifier
	Comparator *cmp;
	if( score_type == "DotProduct" ){ 
		std::cout << "Using DotProduct score function" << std::endl;
		cmp = new DotProduct( ppm_mass_tol, abs_mass_tol );
	}else if( score_type == "OrigSteinDotProduct" ){ 
		std::cout << "Using OrigSteinDotProduct score function" << std::endl;
		cmp = new OrigSteinDotProduct( ppm_mass_tol, abs_mass_tol );
	}else if( score_type == "Jaccard" ){ 
		std::cout << "Using Jaccard score function" << std::endl;
		cmp = new Jaccard( ppm_mass_tol, abs_mass_tol );
	}else if( score_type == "Combined" ){ 
		std::cout << "Using Combined score function" << std::endl;
		cmp = new Combined( ppm_mass_tol, abs_mass_tol );
	}else{ 
		std::cout << "Unknown comparator type: " << score_type << std::endl;
		exit(1);
	}

	//Read in the parameters and set up the identifier
	Identifier *identifier;
	Param *param; NNParam *nn_param;
	if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION ){
		nn_param = new NNParam( param_filename );
		identifier = new Identifier( nn_param, &cfg, cmp, prob_thresh_for_prune );
	}else{
		param = new Param(param_filename);
		identifier = new Identifier( param, &cfg, cmp, prob_thresh_for_prune);
	}

	std::cout << "TARGET ID: " << target_id << std::endl;

	//Rank the candidates
	identifier->rankCandidatesForSpecMatch( candidates, targetData.getSpectra(), output_msp_or_mgf, apply_postprocessing == 1 );
	delete identifier;

	//Keep only the num_highest
	if( num_highest > 0 && num_highest < (int)candidates.size() ) candidates.resize( num_highest );

	//Set up the output (to file or stdout)
	std::streambuf * buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Report the results
	reportResults( candidates, out );
	
	delete cmp;

	return(0);    
}

void readInCandidates( std::vector<Candidate> &candidates, std::string &candidate_file){

	std::string line, smiles_or_inchi, id;
	std::ifstream ifs ( candidate_file.c_str() , std::ifstream::in );

	if( !ifs.good() ){ 
		std::cout << "Could not open input file " << candidate_file << std::endl;
		return;
	}

	while( ifs.good() ){

		getline( ifs, line );
		if( line.size() < 3 ) continue;

		std::stringstream ss(line);
		ss >> id >> smiles_or_inchi;

		candidates.push_back( Candidate( id, smiles_or_inchi ) );
	}

}

void reportResults( std::vector<Candidate> &candidates, std::ostream &out ){

	std::vector<Candidate>::iterator it = candidates.begin();
	out << std::setprecision(8);
	for( int rank = 1; it != candidates.end(); ++it, rank++ ){
		out << rank << " ";
		out << it->getScore() << " ";
		out << *it->getId() << " ";
		out << *it->getSmilesOrInchi() << std::endl;
	}

}