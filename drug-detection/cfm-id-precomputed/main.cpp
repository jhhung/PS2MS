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
#include <map>

void readInCandidates( std::vector<PrecomputedCandidate> &candidates, std::string &candidate_file, std::map<std::string, MspReader> &msps);
void reportResults( std::vector<PrecomputedCandidate> &candidates, std::ostream &out );

int main(int argc, char *argv[])
{
	std::map<std::string, MspReader> msps;
	bool to_stdout = true;
	std::string output_filename;
	int num_highest = -1;
	double abs_mass_tol = 0.01, ppm_mass_tol = 10.0;
	std::string score_type = "Jaccard";

	if (argc < 4 || argc > 9 )
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-id.exe <spectrum_file> <id> <candidate_file> <num_highest> <ppm_mass_tol> <abs_mass_tol> <score_type> <output_filename>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "spectrum_file:" << std::endl << "The filename where the input spectra can be found. This can be a .msp file in which the desired spectrum is listed under a corresponding id (next arg). Or it could be a single file with a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. ";
		std::cout << "e.g." << std::endl << "energy0" << std::endl << "65.02 40.0" << std::endl << "86.11 60.0" << std::endl << "energy1" << std::endl << "65.02 100.0 ... etc" << std::endl;
		std::cout << std::endl << "id:" << std::endl << "An identifier for the target molecule (printed to output, in case of multiple concatenated results, and used to retrieve input spectrum from msp where msp is used)" << std::endl;
		std::cout << std::endl << "candidate_file:" << std::endl << "The filename where the input list of candidate structures can be found - line separated 'id smiles_or_inchi spectrum_file' triples, where the spectrum file stores the precomputed spectra (no spaces allowed)." << std::endl;
		std::cout << std::endl << "num_highest (opt):" << std::endl << "The number of (ranked) candidates to return or -1 for all (if not given, returns all in ranked order)" << std::endl;
		std::cout << std::endl << "ppm_mass_tol (opt):" << std::endl << "The mass tolerance in ppm to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs (if not given defaults to 10ppm)" << std::endl;
		std::cout << std::endl << "abs_mass_tol (opt):" << std::endl << "The mass tolerance in abs Da to use when matching peaks within the dot product comparison - will use higher resulting tolerance of ppm and abs ( if not given defaults to 0.01Da)" << std::endl;
	    std::cout << std::endl << "score_type (opt):" << std::endl << "The type of scoring function to use when comparing spectra. Options: Jaccard (default), DotProduct" << std::endl;
		std::cout << std::endl << "output_filename (opt):" << std::endl << "The filename of the output file to write to (if not given, prints to stdout)" << std::endl;
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
	if( argc > 7 ) score_type = argv[7];
	if( argc > 8 ){
		output_filename = argv[8];
		to_stdout = false;
	}

	//Dummy param and config (since the spectra are already computed)
	config_t cfg; initDefaultConfig(cfg);
	std::vector<std::string> blank_feature_list;
	Param param(blank_feature_list, 1);


	//Read in the input spectrum
	MolData targetData( target_id.c_str(), "Unknown", &cfg );
	if( spectrum_file.substr(spectrum_file.size()-4, 4) == ".msp" ){
		MspReader msp( spectrum_file.c_str(), "" );
		targetData.readInSpectraFromMSP( msp );
	}else targetData.readInSpectraFromFile( spectrum_file );	


	//Fetch the list of candidates
	std::vector<PrecomputedCandidate> candidates;
	readInCandidates( candidates, candidate_file, msps);

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
	Identifier identifier( &param, &cfg, cmp, 0.0 );

	//Rank the candidates
	identifier.rankPrecomputedCandidatesForSpecMatch( candidates, targetData.getSpectra() );

	//Keep only the num_highest
	if( num_highest > 0 && num_highest < (int)candidates.size() ) candidates.resize( num_highest );

	//Set up the output (to file or stdout)
	std::streambuf *buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Report the results
	out << "ID: " << target_id << std::endl;
	reportResults( candidates, out );
	out << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$";
	
	delete cmp;

	return(0);    
}

void readInCandidates( std::vector<PrecomputedCandidate> &candidates, std::string &candidate_file, std::map<std::string, MspReader> &msps){

	std::string line, spectrum_filename, smiles_or_inchi, id;
	std::ifstream ifs ( candidate_file.c_str() , std::ifstream::in );

	if( !ifs.good() ){ 
		std::cout << "Could not open input file " << candidate_file << std::endl;
		return;
	}

	while( ifs.good() ){

		getline( ifs, line );
		if( line.size() < 3 ) continue;

		std::stringstream ss(line);
		ss >> id >> smiles_or_inchi >> spectrum_filename;

		if( spectrum_filename.substr(spectrum_filename.size()-4, 4) == ".msp" ){
			if( msps.find( spectrum_filename ) == msps.end() ){
				msps[spectrum_filename] = MspReader(spectrum_filename.c_str(), "", true);
			}
			candidates.push_back( PrecomputedCandidate( id, smiles_or_inchi, &msps[spectrum_filename] ) );
		}
		else
			candidates.push_back( PrecomputedCandidate( id, smiles_or_inchi, spectrum_filename ) );
		
	}
	ifs.close();

}

void reportResults( std::vector<PrecomputedCandidate> &candidates, std::ostream &out ){

	std::vector<PrecomputedCandidate>::iterator it = candidates.begin();
	out << std::setprecision(8);
	for( int rank = 1; it != candidates.end(); ++it, rank++ ){
		out << rank << " ";
		out << it->getScore() << " ";
		out << *it->getId() << " ";
		out << *it->getSmilesOrInchi() << std::endl;
	}

}