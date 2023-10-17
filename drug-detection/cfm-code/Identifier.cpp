/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Identifier.cpp
#
# Description: 	Identifier class for ranking candidate structures according
#				to their match with a set of target spectra.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Identifier.h"

#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <fstream>
#include <sstream>
#include <iostream>

static bool sort_candidates(const Candidate &u, const Candidate &v){
	return u.getScore() > v.getScore();
}

static bool sort_pcandidates(const PrecomputedCandidate &u, const PrecomputedCandidate &v){
	return u.getScore() > v.getScore();
}


//Ranks the list of candidates according to the match between their predicted spectra and the target
void Identifier::rankCandidatesForSpecMatch( std::vector<Candidate> &candidates, const std::vector<Spectrum> *target_spectra, std::string &output_spectra_filename, bool post_process_spectra, bool output_all_scores ){

	int output_mode = NO_OUTPUT_MODE;

	//Set up the output MSP or MGF file (if selected)
	std::ostream *spec_out;
	std::ofstream of_spec; 
	if( output_spectra_filename.size() > 4 ){
		if( output_spectra_filename.substr(output_spectra_filename.size()-4, 4) == ".msp" )
			output_mode = MSP_OUTPUT_MODE;
		else if( output_spectra_filename.substr(output_spectra_filename.size()-4, 4) == ".mgf" )
			output_mode = MGF_OUTPUT_MODE;
		else{
			std::cout << "Unknown output file type (expecting .msp or .mgf):" << output_spectra_filename << std::endl;
			throw std::exception();		
		}
		of_spec.open( output_spectra_filename.c_str());
		if( !of_spec.is_open() ){
			std::cout << "Error: Could not open candidate spectrum output file " << output_spectra_filename << std::endl;
			throw std::exception();
		}
		std::streambuf *buf = of_spec.rdbuf();
		spec_out = new std::ostream(buf);
	}

	double ppm_mass_tol = cmp->getPPMTol();
	double abs_mass_tol = cmp->getAbsTol();

	//All Comparators
	boost::ptr_vector<Comparator> cmps;
	if( output_all_scores ){
		cmps.push_back( new Recall( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new Precision( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new WeightedRecall( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new WeightedPrecision( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new Jaccard( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new WeightedJaccard( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new DotProduct( ppm_mass_tol, abs_mass_tol ));
		cmps.push_back( new OrigSteinDotProduct( ppm_mass_tol, abs_mass_tol ));
	}

	//Compute the scores for each candidate
	std::vector<Candidate>::iterator it = candidates.begin();
	for( ; it != candidates.end(); ++it ){
	
		LikelyFragmentGraphGenerator *fgen;
		double score = 0.0;
		try{

			//Create the MolData structure with the input
			MolData moldata( it->getId()->c_str(), it->getSmilesOrInchi()->c_str(), cfg );
	
			//Calculate the pruned FragmentGraph
			if( cfg->theta_function == NEURAL_NET_THETA_FUNCTION )
				fgen = new LikelyFragmentGraphGenerator(nn_param, cfg, prob_thresh_for_prune ); 
			else
				fgen = new LikelyFragmentGraphGenerator(param, cfg, prob_thresh_for_prune ); 
			moldata.computeLikelyFragmentGraphAndSetThetas(*fgen, prob_thresh_for_prune, cfg->include_isotopes);
			double precursor_mass = moldata.getFragmentGraph()->getFragmentAtIdx(0)->getMass();

			//Predict the spectra (and post-process, use existing thetas)
			moldata.computePredictedSpectra( *param, false, true );
			if( output_all_scores ) std::cout << *it->getId() << ":";
			
			score = 0.0;
			for( int postprocess = 0; postprocess <= 1; postprocess++ ){

				if(postprocess) moldata.postprocessPredictedSpectra();
				if(abs_mass_tol == 0.5 && !postprocess ) moldata.quantisePredictedSpectra(0);	//Integer precision data, combine integer masses.
			
				if(postprocess == (int)post_process_spectra){
				
					//Write the predicted spectra to the output file
					if( output_mode == MSP_OUTPUT_MODE ) moldata.writePredictedSpectraToMspFileStream( *spec_out );
					else if( output_mode == MGF_OUTPUT_MODE ) moldata.writePredictedSpectraToMgfFileStream( *spec_out );

					//Compute the score for the current candidate:
					// - The score is the sum of the comparison scores between all the target and predicted spectra
					for( unsigned int energy = 0; energy < target_spectra->size(); energy++ ){
						score += cmp->computeScore( &((*target_spectra)[energy]), moldata.getPredictedSpectrum(energy));
					}
				}

				//Compute and report all comparator scores: (hack to save computation during testing - enabled by output_all_scores flag (disabled by default))
				boost::ptr_vector<Comparator>::iterator itc = cmps.begin();
				for( ; itc != cmps.end(); ++itc ){
					double tmp_score = 0.0;
					for( unsigned int energy = 0; energy < target_spectra->size(); energy++ ){
						tmp_score += itc->computeScore( &((*target_spectra)[energy]), moldata.getPredictedSpectrum(energy));
					}
					std::cout << " " << tmp_score;
				}
			}
			if( output_all_scores ) std::cout << std::endl;

		}
		catch( RDKit::MolSanitizeException e ){
			std::cout << "Could not sanitize " << *it->getSmilesOrInchi() << std::endl;
		}
		catch( RDKit::SmilesParseException pe ){
			std::cout << "Could not parse " << *it->getSmilesOrInchi() << std::endl;		
		}
		catch( FragmentGraphGenerationException e ){
			std::cout << "Could not compute fragmentation graph for " << *it->getSmilesOrInchi() << std::endl;
		}		
		catch( FragmentGraphTimeoutException te ){
			std::cout << "Timeout computing fragmentation graph for input: " << *it->getSmilesOrInchi() << std::endl;	
		}
		catch( std::exception e ){
			std::cout << "Exception occurred:" << e.what() << std::endl;
		}

		it->setScore(score);
		delete fgen;
	}
	if( output_mode == MSP_OUTPUT_MODE || output_mode == MGF_OUTPUT_MODE ){ 
		of_spec.close();
		delete spec_out;
	}

	//Sort the candidates by decreasing score
	std::sort( candidates.begin(), candidates.end(), sort_candidates );
}



//Ranks the list of candidates according to the match between their predicted spectra and the target
void Identifier::rankPrecomputedCandidatesForSpecMatch( std::vector<PrecomputedCandidate> &candidates, const std::vector<Spectrum> *target_spectra ){


	//Compute the scores for each candidate
	std::vector<PrecomputedCandidate>::iterator it = candidates.begin();
	for( ; it != candidates.end(); ++it ){
	
		double score = 0.0;

		//Load the predicted spectra
		if( !it->hasSpectra() ){ 
			MolData moldata( it->getId()->c_str(), it->getSmilesOrInchi()->c_str(), cfg );
			moldata.readInSpectraFromFile(it->getSpectrumFilename()->c_str(), true);
			for( unsigned int energy = 0; energy < target_spectra->size(); energy++ )
				score += cmp->computeScore( &((*target_spectra)[energy]), moldata.getPredictedSpectrum(energy) );
		}
		else{
			for( unsigned int energy = 0; energy < target_spectra->size(); energy++ )
				score += cmp->computeScore( &((*target_spectra)[energy]), &(*it->getSpectra())[energy]);			
		}

		//Compute the score for the current candidate:
		// - The score is the sum of the comparison scores between all the target and predicted spectra
		it->setScore(score);

	}

	//Sort the candidates by decreasing score
	std::sort( candidates.begin(), candidates.end(), sort_pcandidates );

}
