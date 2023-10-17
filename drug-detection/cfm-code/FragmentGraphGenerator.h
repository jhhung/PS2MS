/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraphGenerator.h
#
# Description: 	FragmentGraphGenerator class for generating a fragment tree.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FRAG_GEN_H__
#define __FRAG_GEN_H__

#include "FragmentGraph.h"
#include "FragmentTreeNode.h"
#include "Features.h"
#include "Param.h"
#include "NNParam.h"
#include "Config.h"
#include <GraphMol/RWMol.h>

static const int MAX_FRAGMENTS_PER_MOLECULE = 100000;
static const int MAX_TRANSITIONS_PER_MOLECULE = 1000000;

class FragmentGraphGenerationException: public std::exception{

	virtual const char* what() const throw(){
		return "Error during fragment graph generation, unable to proceed.";
	}
};

class IonizationException: public std::exception{

	virtual const char* what() const throw(){
		return "Error during molecule ionization, unable to proceed.";
	}
};

class FragmentGraphMaxSizeExceededException: public std::exception{

	virtual const char* what() const throw(){
		return "Molecule exceeded maximum number of fragments or transitions, unable to proceed.";
	}
};

class FragmentGraphTimeoutException: public std::exception{

	virtual const char* what() const throw(){
		return "Molecule fragment graph computation exceeded timeout, unable to proceed.";
	}
};

//Class for generating a FragmentGraph given a starting node in the
//enumerated fragmentation tree and a depth
class FragmentGraphGenerator{
public:
	//Constructor
	FragmentGraphGenerator() : verbose(0), mols_to_fv(false) { fh = new FeatureHelper(); };
	FragmentGraphGenerator(int a_verbose) : verbose(a_verbose), mols_to_fv(false) { fh = new FeatureHelper();};

	//Constructor to use if the transition molecules are to be replaced by a feature vector as
	//soon as they are created (for less memory usage).
	FragmentGraphGenerator(FeatureCalculator *a_fc, bool a_mols_to_fv = true ) :
	fc(a_fc), verbose(0), mols_to_fv(a_mols_to_fv) { fh = new FeatureHelper(fc);};

	//Destructor
	~FragmentGraphGenerator(){delete fh;};

	//Start a graph. Compute can then add to this graph, but it is the caller's 
	//responsibility to delete it
	FragmentGraph *createNewGraph( config_t *cfg );

	//Create the starting node from a smiles or inchi string - responsibility of caller to delete
	FragmentTreeNode *createStartNode( std::string &smiles_or_inchi, int ionization_mode, RDKit::RWMol* rwmol );


	//Compute a FragmentGraph starting at the given node and computing to the depth given.
	//The output will be appended to the current_graph
	void compute( FragmentTreeNode &startnode, int remaining_depth, int parentid, int remaining_ring_breaks );

protected:
	FeatureCalculator *fc;
	FeatureHelper *fh;
	FragmentGraph *current_graph;
	bool mols_to_fv;
	int verbose;

private:

	//Record of previous computations so we know to what depth each fragment has been computed
	std::map<int, int> id_depth_computed_cache;	
	//Helper function - check if the fragment has already been computed to at least this depth
	int alreadyComputed(int id, int remaining_depth);
	
	//Static Helper functions
	static int countExtraElectronPairs( RDKit::RWMol *rwmol, std::vector<int> &output_e_loc  );
	static void applyIonization( RDKit::RWMol *rwmol, int ionization_mode );

};

//Class to be used if a graph is to be pruned as it's created 
//removing anything with probability below a given threshold
class LikelyFragmentGraphGenerator : public FragmentGraphGenerator {
public:
	//Constructor
	LikelyFragmentGraphGenerator(Param *a_param, config_t *a_cfg, double a_prob_thresh ) :
	  cfg( a_cfg ), param(a_param), log_prob_thresh( log(a_prob_thresh) ), is_nn_params(false) 
	  { fc = new FeatureCalculator(*param->getFeatureNames());  mols_to_fv = true; fh = new FeatureHelper(fc);};
	LikelyFragmentGraphGenerator(NNParam *a_param, config_t *a_cfg, double a_prob_thresh ) :
	  cfg( a_cfg ), nnparam(a_param), log_prob_thresh( log(a_prob_thresh) ),  is_nn_params(true) 
	  { fc = new FeatureCalculator(*nnparam->getFeatureNames());  mols_to_fv = true; fh = new FeatureHelper(fc);};

	//Start a graph. Compute can then add to this graph, but it is the caller's 
	//responsibility to delete it
	FragmentGraph *createNewGraph( config_t *cfg );

	//Compute a FragmentGraph starting at the given node and computing to the depth given.
	//The output will be appended to the current_graph
	void compute( FragmentTreeNode &startnode, int remaining_depth, int parentid, double parent_log_prob, int remaining_ring_breaks );

private:
	Param *param;
	NNParam *nnparam;
	bool is_nn_params;
	config_t *cfg;
	double log_prob_thresh;
	time_t start_time;

	//Record of previous computations so we know to what probability each fragment has been computed at
	std::map<int, double> id_prob_computed_cache;	
	int alreadyComputed(int id, double prob_offset);

};

#endif // __FRAG_GEN_H__
