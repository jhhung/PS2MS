/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.c
#
# Description: 	Search for valid fragment configurations using backtrack.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

int main(int argc, char *argv[]);

#include "MolData.h"

#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
	std::string output_filename;
	std::string smiles_or_inchi; 
	int max_depth;

	bool verbose = false;
	bool to_stdout = true;
	bool fullgraph = true;
    
	if (argc != 4 && argc != 5 && argc != 6 )
	{
		std::cout << std::endl << "Usage:" << std::endl;
		std::cout << "fraggraph-gen.exe <smiles or inchi string> <max_depth> <ionization_mode (+,- or *)> ";
		std::cout << "<opt: fullgraph (default) or fragonly> <opt: output_filename (else stdout)>" << std::endl << std::endl;
		std::cout << "ionization_mode: Positive ESI Ionization (+), Negative ESI Ionization (-), Positive EI Ionization (*)" << std::endl << std::endl;
		exit(1);
	}
    smiles_or_inchi = argv[1];
    max_depth = atoi(argv[2]);

	int ionization_mode = 0;
	int ionization_mode_symbol = *(argv[3]);
	if( ionization_mode_symbol == '+' ) ionization_mode = POSITIVE_ESI_IONIZATION_MODE;
	else if( ionization_mode_symbol == '-' ) ionization_mode = NEGATIVE_ESI_IONIZATION_MODE;
	else if( ionization_mode_symbol == '*' ) ionization_mode = POSITIVE_EI_IONIZATION_MODE;
	else{ 
		std::cout << "Unknown setting for ionization mode (expecting +,- or *): " << ionization_mode << std::endl;
		exit(1);
	}
	
	if( argc > 4 ){
		std::string type = argv[4];
		if( type == "fullgraph" ) fullgraph = true;
		else if( type == "fragonly" ) fullgraph = false;
		else std::cout << "Unknown setting for output type: " << type << " (expecting fullgraph or fragonly)" << std::endl;
	}
	if( argc > 5 ){ 
		output_filename = argv[5];
		to_stdout = false;
	}
	
	//Set up the output (to file or stdout)
	std::streambuf * buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Set up the configuration
	config_t cfg;
	initDefaultConfig(cfg);
	cfg.ionization_mode = ionization_mode;
	cfg.fg_depth = max_depth;
	cfg.include_h_losses = true; 

	time_t before = time( NULL );

	//Run the fragmentation procedure
	std::vector<std::string> feature_names;
	feature_names.push_back(std::string("BreakAtomPair"));
	FeatureCalculator fc( feature_names );
	MolData moldata("FRAGGRAPH-GEN", smiles_or_inchi.c_str(), &cfg);
	moldata.computeFragmentGraphAndReplaceMolsWithFVs(&fc, true);
	const FragmentGraph *graph = moldata.getFragmentGraph();

	time_t after = time( NULL );

	//Write to output
	if( fullgraph ) graph->writeFullGraph( out );
	else graph->writeFragmentsOnly( out );

	//If we're not writing output to stdout, report the statistics of the graph.
	if(!to_stdout){
		std::cout << std::endl << "Output written to " << output_filename << std::endl;
		std::cout << "Depth " << max_depth << ": " << graph->getNumFragments() << " Fragments, " << graph->getNumTransitions() << " Transitions";
		std::cout << ". Time Elaspsed = " << (after - before) << " Seconds" << std::endl;
	}

	return(0);    
}
