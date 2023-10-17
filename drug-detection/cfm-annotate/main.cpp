/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Given a spectrum and a structure, (and optionally, a trained 
#				CFM model), annotate the fragments associated with each peak
#				in the spectrum.
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
#include "IPFP.h"

#include <iostream>
#include <fstream>
#include <string>

void concatBeliefs(beliefs_t *all_beliefs, beliefs_t *beliefs);

int main(int argc, char *argv[])
{
	bool to_stdout = true;
	std::string output_filename;
	std::string param_filename = "none";
	std::string config_filename = "param_config.txt";
	int num_highest = -1;
	double ppm_mass_tol = DEFAULT_PPM_MASS_TOL;
	double abs_mass_tol = DEFAULT_ABS_MASS_TOL;

	if (argc < 3 || argc > 9 )
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-annotate.exe <smiles_or_inchi> <spectrum_file> <id> <ppm_mass_tol> <abs_mass_tol> <param_filename> <config_filename> <output_filename>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "smiles_or_inchi: " << std::endl << "The smiles or Inchi string for the input molecule" << std::endl;  
		std::cout << std::endl << "spectrum_file:" << std::endl << "The filename where the input spectra can be found as a list of peaks 'mass intensity' delimited by lines, with either 'low','med' and 'high' lines beginning spectra of different energy levels, or 'energy0', 'energy1', etc. ";
		std::cout << "e.g." << std::endl << "energy0" << std::endl << "65.02 40.0" << std::endl << "86.11 60.0" << std::endl << "energy1" << std::endl << "65.02 100.0 ... etc" << std::endl;
		std::cout << std::endl << "id:" << std::endl << "An identifier for the target molecule (Used to retrieve input spectrum from msp (if used). Otherwise not used but printed to output, in case of multiple concatenated results)" << std::endl;
		std::cout << std::endl << "ppm_mass_tol (opt):" << std::endl << "The mass tolerance in ppm to use when matching peaks - will use higher resulting tolerance of ppm and abs (if not given defaults to value in the config file, or 10ppm if not specified there)" << std::endl;
		std::cout << std::endl << "abs_mass_tol (opt):" << std::endl << "The mass tolerance in abs Da to use when matching peaks - will use higher resulting tolerance of ppm and abs ( if not given defaults to value in the config file, 0.01Da if not specified there)" << std::endl;
		std::cout << std::endl << "param_filename (opt):" << std::endl << "The filename where the parameters of a trained cfm model can be found (if not given or set to 'none', assumes no parameters set, so all breaks equally likely)" << std::endl;
		std::cout << std::endl << "config_filename (opt):" << std::endl << "The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory)" << std::endl;
		std::cout << std::endl << "output_filename (opt):" << std::endl << "The filename of the output file to write to (if not given, prints to stdout)" << std::endl;
		exit(1);
	}

	std::string smiles_or_inchi = argv[1];
	std::string spectrum_file = argv[2];
	std::string target_id = argv[3];
	if( argc > 4 ){ 
		try{ ppm_mass_tol = boost::lexical_cast<float>(argv[4]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid ppm_tol: " << argv[4] << std::endl;
			exit(1);
		}
	}
	if( argc > 5 ){
		try{ abs_mass_tol = boost::lexical_cast<float>(argv[5]); }
		catch(boost::bad_lexical_cast e){
			std::cout << "Invalid abs_tol: " << argv[5] << std::endl;
			exit(1);
		}
	}
	if( argc > 5 ) param_filename = argv[6];
	if( argc > 6 ) config_filename = argv[7];
	if( argc > 8 ){
		output_filename = argv[8];
		to_stdout = false;
	}

	//Initialise model configuration
	config_t cfg;
	initConfig( cfg, config_filename );

	//Read in the input spectrum
	MolData moldata( target_id.c_str(), smiles_or_inchi.c_str(), &cfg );
	if( spectrum_file.substr(spectrum_file.size()-4, 4) == ".msp" ){
		MspReader msp( spectrum_file.c_str(), "" );
		moldata.readInSpectraFromMSP( msp );
	}else moldata.readInSpectraFromFile( spectrum_file );

	//If the mass tolerances are specified on the command line, overwrite the
	//ones in the config file
	if( argc > 4 ) cfg.ppm_mass_tol = ppm_mass_tol;
	if( argc > 5 ) cfg.abs_mass_tol = abs_mass_tol;

	//Read in the parameters or create a blank set
	Param *param;
	if( param_filename == "none" ){
		std::vector<std::string> fnames;	//empty feature set
		param = new Param(fnames, moldata.getNumSpectra());
	}else{ 
		if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
			param = new NNParam( param_filename );
		else
			param = new Param(param_filename);	
		param = new Param( param_filename );
	}

	if( param->getNumEnergyLevels() != moldata.getNumSpectra() ){
		std::cout << "Mismatch between parameter set and spectra. Parameters have ";
		std::cout << param->getNumEnergyLevels() << " energy levels whereas there are ";
		std::cout << moldata.getNumSpectra() << " spectra" << std::endl;
		exit(1);
	}

	std::cout << "TARGET ID: " << target_id << std::endl;

	//Compute the fragmentation graph with transition probabilities
	FeatureCalculator fc( *param->getFeatureNames() );
	moldata.computeFragmentGraphAndReplaceMolsWithFVs(&fc, true);
	moldata.computeTransitionThetas(*param);
	moldata.computeTransitionProbabilities();

	moldata.removePeaksWithNoFragment( cfg.abs_mass_tol, cfg.ppm_mass_tol );

	//Apply the peak evidence and compute the beliefs
	beliefs_t beliefs;
	if( cfg.use_single_energy_cfm ){
		
		//Concat single energy beliefs into one since subsequent
		//functions just check for any beliefs above threshold
		for( int energy = 0; energy < cfg.dv_spectrum_depths.size(); energy++ ){
			config_t se_cfg;
			initSingleEnergyConfig(se_cfg, cfg, energy);
			Inference infer( &moldata, &se_cfg );
			beliefs_t sbeliefs;
			infer.calculateBeliefs( sbeliefs );
			concatBeliefs(&beliefs, &sbeliefs);
		}
	}
	else{
		IPFP ipfp( &moldata, &cfg); 
		beliefs_t *sbeliefs = ipfp.calculateBeliefs();
		concatBeliefs(&beliefs, sbeliefs);
	}

	//Process the beliefs to extract the fragmentation tree that occurred
	double log_belief_thresh = std::log(0.00001);
	moldata.computeEvidenceFragmentGraph(&beliefs, log_belief_thresh);

	//Match peaks to fragments in the reduced tree
	moldata.readInSpectraFromFile( spectrum_file );
	moldata.annotatePeaks(cfg.abs_mass_tol, cfg.ppm_mass_tol);

	//Set up the output stream
	std::streambuf * buf;
	std::ofstream of;
	if( !to_stdout ) {
		of.open(output_filename.c_str());
		if( !of.is_open() ) 
			std::cout << "Warning: Trouble opening output file" << std::endl;
		buf = of.rdbuf();
	} else buf = std::cout.rdbuf();
	std::ostream out(buf);

	//Output the peak annotations and fragment graph
	moldata.outputSpectra( out, "Experimental", true );
	out << std::endl;
	moldata.getEvidenceFragmentGraph()->writeFullGraph(out);

	delete param;
	return(0);    
}


void concatBeliefs(beliefs_t *all_beliefs, beliefs_t *beliefs){

	if( all_beliefs->ps.size() == 0 ){
		all_beliefs->ps.reserve( beliefs->ps.size() );
		all_beliefs->tn.reserve( beliefs->tn.size() );
		all_beliefs->ps.insert(all_beliefs->ps.end(), beliefs->ps.begin(), beliefs->ps.end() );
		all_beliefs->tn.insert(all_beliefs->tn.end(), beliefs->tn.begin(), beliefs->tn.end() );
	}
	else{
		std::vector<std::vector<double> >::iterator it1, it2;
		it1 = all_beliefs->ps.begin(); it2 = beliefs->ps.begin();
		for( ; it2 != beliefs->ps.end(); ++it1, ++it2 )
			it1->insert( it1->end(), it2->begin(), it2->end() );
		it1 = all_beliefs->tn.begin(); it2 = beliefs->tn.begin();
		for( ; it2 != beliefs->tn.end(); ++it1, ++it2 )
			it1->insert( it1->end(), it2->begin(), it2->end() );
	}


}