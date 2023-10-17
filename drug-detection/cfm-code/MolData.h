/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MolData.h
#
# Description: 	Class to hold the input data belonging to a molecule:
#					 - An ID and smiles/inchi
#					 - (optional) A cross-validation group identifier
#					 - (optional) A computed fragmentation graph
#					 - (optional) A computed set of features corresponding to that graph
#					 - (optional) A computed set of theta values for that graph
#				     - (optional) A computed set of transition probabilities using those thetas.
#					 - (optional) A set of spectra
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MOLDATA_H__
#define __MOLDATA_H__

#include <DataStructs/ExplicitBitVect.h>

#include "Features.h"
#include "FragmentGraph.h"
#include "FragmentGraphGenerator.h"
#include "MspReader.h"
#include "Param.h"
#include "NNParam.h"
#include "Config.h"
#include "Message.h"
#include "Spectrum.h"

class MolData{
public:
	MolData( std::string &an_id, std::string &an_smiles_or_inchi, int a_group, config_t *a_cfg  ) :
	  id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(a_group), graph_computed(0), ev_graph_computed(0), cfg(a_cfg), mol(nullptr) {}; 
	MolData( const char *an_id, const char *an_smiles_or_inchi, config_t *a_cfg  ) :
	  id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(0), graph_computed(0), ev_graph_computed(0), cfg(a_cfg), mol(nullptr) {}; 
	MolData( const std::string &an_id, const std::string &an_smiles_or_inchi, config_t *a_cfg, RDKit::RWMOL_SPTR a_mol  ) :
	  id(an_id), smiles_or_inchi(an_smiles_or_inchi), group(0), graph_computed(0), ev_graph_computed(0), cfg(a_cfg), mol(new RDKit::RWMol(*a_mol, true)), predicted_spectra(1) {}; 
		
	//Access functions
	const FragmentGraph *getFragmentGraph() const {return fg;};
	int hasComputedGraph() const { return graph_computed; };
	const EvidenceFragmentGraph *getEvidenceFragmentGraph() const {return ev_fg;};
	const Spectrum *getSpectrum(int energy) const {return &(spectra[energy]);};
	const std::vector<Spectrum> *getSpectra() const {return &spectra;};
	const Spectrum *getPredictedSpectrum(int energy) const {return &(predicted_spectra[energy]);};
	unsigned int getNumSpectra() const {return spectra.size();};
	unsigned int getNumPredictedSpectra() const {return predicted_spectra.size();};
	const FeatureVector *getFeatureVectorForIdx(int index)  const { return fvs[index]; };
	double getThetaForIdx(int energy, int index)  const
	{ return thetas[energy][index];};
	double getLogTransitionProbForIdx(int energy, int index)  const
	{ return log_probs[energy][index];};
	double getLogPersistenceProbForIdx(int energy, int index)  const
	{ return log_probs[energy][fg->getNumTransitions() + index];};
	int getGroup() const { return group; };
	void setGroup( int val ){ group = val; };
	std::string getId() const { return id; };
	std::string getSmilesOrInchi() const { return smiles_or_inchi; };
	double getMolecularWeight() const;
	int getIonizationMode() const { return cfg->ionization_mode;};
	RDKit::RWMol get_mol() const { return *mol; };

	void readInSpectraFromFile( const std::string &filename, bool readToPredicted = false );
	void readInSpectraFromMSP( MspReader &msp, bool readToPredicted = false );
	std::vector<Spectrum>* readInSpectraFromMSPFileStream( std::istream &msp, int num_peaks, bool readToPredicted = false );
	void cleanSpectra(double abs_tol, double ppm_tol);
	void freeSpectra(){ std::vector<Spectrum>().swap(spectra); 
			std::vector<Spectrum>().swap(predicted_spectra); };

	//Spectrum Related Functions
	void removePeaksWithNoFragment( double abs_tol, double ppm_tol );
	bool hasEmptySpectrum() const;
	void writePredictedSpectraToFile( std::string &filename );
	void writePredictedSpectraToMspFileStream( std::ostream &out );
	void writePredictedSpectraToMgfFileStream( std::ostream &out );
	void writeFullEnumerationSpectrumToFile( std::string &filename );
	void writeFullEnumerationSpectrumToMspFileStream( std::ostream &out );
	void outputSpectra( std::ostream &out, const char*spec_type, bool do_annotate = false );
	
	//More memory efficient alternative to calling computeFragmentGraph and 
	//then computeFeatureVectors with deleteMols = true
	void computeFragmentGraphAndReplaceMolsWithFVs( FeatureCalculator *fc, bool retain_smiles = false);
	
	//Save/load state functions
	void readInFVFragmentGraph( std::string &fv_filename);
	void readInFVFragmentGraphFromStream( std::istream &ifs);
	void writeFVFragmentGraph( std::string &fv_filename);
	void writeFVFragmentGraphToStream( std::ofstream &out);

	//Replaces computeFragmentGraph, computeFeatureVectors and computeTransitionThetas
	//below (delteMols = true), pruning according to prob_thresh_for_prune value.
	void computeLikelyFragmentGraphAndSetThetas( LikelyFragmentGraphGenerator &fgen, double prob_thresh_for_prune, bool retain_smiles = false );

	//Note that the following should be called in this order
	//since each one assumes all previous have already been called.
	void computeFragmentGraph( );
	void computeFragmentGraph( FeatureCalculator *fc  );	//Use this option if computing features next
	void computeFeatureVectors( FeatureCalculator *fc, bool deleteMols = false );
	void computeTransitionThetas( Param &param );
	void computeTransitionProbabilities();
	void computePredictedSpectra( Param &param, bool postprocess = false, bool use_existing_thetas = false );

	void postprocessPredictedSpectra(double perc_thresh = 80.0, int min_peaks = 5, int max_peaks = 30);
	void quantisePredictedSpectra(int num_dec_places);
	void quantiseMeasuredSpectra(int num_dec_places);

	//Function to compute a much reduced fragment graph containing only those
	//fragmentations as actually occur in the spectra, based on a computed set of beliefs
	//thresholding inclusion in the graph by the provided belief_thresh value (log domain)
	void computeEvidenceFragmentGraph( beliefs_t *beliefs, double log_belief_thresh );
	void annotatePeaks(double abs_tol, double ppm_tol, bool prune_deadends = true);

	~MolData();

	std::string fp;
	std::string predicted_fp;

protected:	//These items are protected rather than private for access during tests.
	int group;
	std::string id;
	std::string smiles_or_inchi;
	FragmentGraph *fg;
	EvidenceFragmentGraph *ev_fg;
	int graph_computed;
	int ev_graph_computed;
	std::vector<Spectrum> spectra;
	std::vector<Spectrum> predicted_spectra;
	std::vector<FeatureVector *> fvs;
	std::vector<std::vector<double> > thetas;
	std::vector<std::vector<double> > log_probs;
	config_t *cfg;
	RDKit::RWMol* mol;

	//General utilty functions
	void computeGraphWithGenerator( FragmentGraphGenerator &fgen );
	void getEnumerationSpectraMasses( std::vector<double> &output_masses );
	void computePredictedSingleEnergySpectra( Param &param, bool postprocess, bool use_existing_thetas );
	void translatePeaksFromMsgToSpectra( Spectrum &out_spec, Message *msg );
	void translatePeaksFromMsgToSpectraWithIsotopes( Spectrum &out_spec, Message *msg );
	void computeFragmentEvidenceValues(std::vector<double> &evidence, int frag_idx, const beliefs_t *beliefs);

};

#endif // __MOLDATA_H__
