/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Identifier.h
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

#ifndef __IDENTIFIER_H__
#define __IDENTIFIER_H__

#include "MolData.h"
#include "Param.h"
#include "Config.h"
#include "Comparators.h"
#include "MspReader.h"

class Candidate{
public:
	Candidate(){};
	Candidate( std::string &an_id, std::string &a_smiles_or_inchi ) :
	  id( an_id ), smiles_or_inchi( a_smiles_or_inchi ), score( 0.0 ) {};

	//Access Functions
	std::string *const getId(){ return &id; };
	std::string *const getSmilesOrInchi(){ return &smiles_or_inchi; };
	double getScore() const { return score; };
	void setScore(double val){ score = val; };

private:
	std::string id;
	std::string smiles_or_inchi;
	double score;
};

class PrecomputedCandidate{
public:
	PrecomputedCandidate(){};
	PrecomputedCandidate( std::string &an_id, std::string &a_smiles_or_inchi, std::string &a_filename ) :
	  id( an_id ), smiles_or_inchi( a_smiles_or_inchi ), spectrum_filename( a_filename ), score( 0.0 ), has_spectra(false) {};
	PrecomputedCandidate( std::string &an_id, std::string &a_smiles_or_inchi, MspReader *a_msp ) :
	  id( an_id ), smiles_or_inchi( a_smiles_or_inchi ), score( 0.0 ), has_spectra(true) { spectra = a_msp->fetchSpectrumForId( an_id.c_str() ); };

	//Access Functions
	std::string *const getId(){ return &id; };
	std::string *const getSpectrumFilename(){ return &spectrum_filename; };
	std::string *const getSmilesOrInchi(){ return &smiles_or_inchi; };
	const std::vector<Spectrum> *getSpectra(){ return spectra; };
	double getScore() const { return score; };
	void setScore(double val){ score = val; };
	bool hasSpectra(){ return has_spectra; };

private:
	std::string id;
	std::string spectrum_filename;
	double score;
	std::string smiles_or_inchi;
	const std::vector<Spectrum> *spectra;
	bool has_spectra;

};

class Identifier{
public:
	//Constructor
	Identifier( NNParam *a_nnparam, config_t *a_cfg, Comparator *a_cmp, double a_prob_thresh_for_prune ) : 
	  nn_param( a_nnparam ), cfg( a_cfg ), cmp( a_cmp ), prob_thresh_for_prune(a_prob_thresh_for_prune) {};
	Identifier( Param *a_param, config_t *a_cfg, Comparator *a_cmp, double a_prob_thresh_for_prune ) : 
	  param( a_param ), cfg( a_cfg ), cmp( a_cmp ), prob_thresh_for_prune(a_prob_thresh_for_prune) {};

	//Ranks the list of candidates according to the match between their predicted spectra and the target
	void rankCandidatesForSpecMatch( std::vector<Candidate> &candidates, const std::vector<Spectrum> *target_spectra, std::string &output_spectra_filename, bool post_process_spectra, bool output_all_scores = false );
	void rankPrecomputedCandidatesForSpecMatch( std::vector<PrecomputedCandidate> &candidates, const std::vector<Spectrum> *target_spectra );

private:

	Param *param;
	NNParam *nn_param;
	config_t *cfg;
	Comparator *cmp;
	double prob_thresh_for_prune;

};

#endif // __IDENTIFIER_H__