/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MolData.cpp
#
# Description: 	Class to hold the input data belonging to a molecule:
#					 - An ID and smiles/inchi
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

#include "MolData.h"
#include "FragmentGraphGenerator.h"
#include "Inference.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <RDGeneral/types.h>
#include <INCHI-API/inchi.h>

#include <string>
#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>

double MolData::getMolecularWeight() const{
	romol_ptr_t mol = createMolPtr( smiles_or_inchi.c_str() );
	return getMonoIsotopicMass( mol );
}

void MolData::readInFVFragmentGraphFromStream( std::istream &ifs ){

	FragmentGraphGenerator fgen;
	fg = fgen.createNewGraph(cfg);
	fg->readFeatureVectorGraph( ifs );
	
	//Copy all the feature vector pointers up into the mol data
	fvs.resize( fg->getNumTransitions() );
	for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
		fvs[i] = fg->getTransitionAtIdx(i)->getTmpFV();	
	
	graph_computed = 1;
}

void MolData::readInFVFragmentGraph( std::string &fv_filename){

	std::ifstream ifs ( fv_filename.c_str(), std::ifstream::in | std::ios::binary );

	if(!ifs) std::cout << "Could not open file " << fv_filename << std::endl;
	readInFVFragmentGraphFromStream(ifs);
}

void MolData::writeFVFragmentGraphToStream( std::ofstream &out){
	fg->writeFeatureVectorGraph(out, cfg->include_isotopes);
}

void MolData::writeFVFragmentGraph( std::string &fv_filename){

	std::ofstream out;
	out.open(fv_filename.c_str(), std::ios::out | std::ios::binary );
	if( !out.is_open() ){
		std::cout << "Warning: Trouble opening output fv fragment graph file: " << fv_filename << std::endl;
	}else{
		writeFVFragmentGraphToStream( out );
	}
}


void MolData::computeGraphWithGenerator( FragmentGraphGenerator &fgen ){
	
	fg = fgen.createNewGraph(cfg);
	FragmentTreeNode *startnode = fgen.createStartNode( smiles_or_inchi, cfg->ionization_mode, mol );
	if( cfg->do_prelim_bfs && !cfg->allow_frag_detours && cfg->fg_depth > 1 ) 
		fgen.compute( *startnode, 1, -1, cfg->max_ring_breaks ); //Initial breadth-first to depth 1 (faster?)
	fgen.compute( *startnode, cfg->fg_depth, -1, cfg->max_ring_breaks );
	if( !cfg->allow_frag_detours) fg->removeDetours();
	delete startnode;
	graph_computed = 1;

}

void MolData::computeFragmentGraph(){
	//Note: Do not use this if you want to compute fvs after - FeatureHelper will not have run correctly. Use other computeFragmentGraph instead.
	FragmentGraphGenerator fgen;
	computeGraphWithGenerator( fgen );
}

void MolData::computeFragmentGraph( FeatureCalculator *fc ){
	
	FragmentGraphGenerator fgen(fc, false);
	computeGraphWithGenerator( fgen );
}

void MolData::computeFragmentGraphAndReplaceMolsWithFVs( FeatureCalculator *fc, bool retain_smiles){
	
	//Compute the fragment graph, replacing the transition molecules with feature vectors
	FragmentGraphGenerator fgen( fc );
	computeGraphWithGenerator( fgen );

	//Copy all the feature vector pointers up into the mol data
	fvs.resize( fg->getNumTransitions() );
	for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
		fvs[i] = fg->getTransitionAtIdx(i)->getTmpFV();

	//Delete all the fragment smiles (we only need these while we're computing the graph)
	if(!retain_smiles) fg->clearAllSmiles();

}

void MolData::computeLikelyFragmentGraphAndSetThetas( LikelyFragmentGraphGenerator &fgen, double prob_thresh_for_prune, bool retain_smiles ){
	
	//Compute the fragment graph, replacing the transition molecules with feature vectors
	fg = fgen.createNewGraph(cfg);
	boost::shared_ptr<FragmentTreeNode> startnode(fgen.createStartNode( smiles_or_inchi, cfg->ionization_mode, mol ));
	fgen.compute( *startnode, cfg->fg_depth, -1, 0.0, cfg->max_ring_breaks );
	if( !cfg->allow_frag_detours) fg->removeDetours();
	// delete startnode;
	graph_computed = 1;

	//Copy all the theta values up into the mol data and delete the tmp thetas
	const unsigned int num_levels = cfg->spectrum_depths.size();
	thetas.resize(num_levels);
	for( unsigned int energy = 0; energy < num_levels; energy++ ){
		thetas[energy].resize( fg->getNumTransitions() );
		for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
			thetas[energy][i] = (*(fg->getTransitionAtIdx(i)->getTmpThetas()))[energy];
	}
	
	//Delete all the fragment smiles (we only need these while we're computing the graph)
	if(!retain_smiles) fg->clearAllSmiles();

}


//Function to compute a much reduced fragment graph containing only those
//fragmentations as actually occur in the spectra, based on a computed set of beliefs
//thresholding inclusion in the graph by the provided belief_thresh value (log domain)
void MolData::computeEvidenceFragmentGraph( beliefs_t *beliefs, double log_belief_thresh ){

	unsigned int num_transitions = fg->getNumTransitions();
	unsigned int num_fragments = fg->getNumFragments();
	ev_fg = new EvidenceFragmentGraph( cfg );
	
	//Add the root fragment
	Transition null_t;
	std::vector<int> id_lookup(num_fragments,-1);
	std::vector<double> main_ev;
	computeFragmentEvidenceValues(main_ev, 0, beliefs);
	id_lookup[0] = ev_fg->addToGraphDirectNoCheck(EvidenceFragment(*(fg->getFragmentAtIdx(0)),-1,main_ev),&null_t,-1);


	//Add the fragments and transitions if the beliefs are high enough
	std::vector<int> t_added_flags(num_transitions, 0);
	for( int depth = 0; depth < beliefs->tn[0].size(); depth++ ){

		for( unsigned int i = 0; i < num_transitions; i++ ){
			if(t_added_flags[i]) continue;

			if( beliefs->tn[i][depth] >= log_belief_thresh ){
				const Transition *t = fg->getTransitionAtIdx(i);
				const Fragment *f = fg->getFragmentAtIdx( t->getToId() );
				int from_id = id_lookup[t->getFromId()];
				if( id_lookup[t->getToId()] == -1 ){
					std::vector<double> evidence;
					computeFragmentEvidenceValues(evidence, t->getToId(), beliefs );
					int to_id = ev_fg->addToGraphDirectNoCheck(EvidenceFragment(*f,-1,evidence),t,from_id);
					id_lookup[t->getToId()] = to_id;
				}
				else if( id_lookup[t->getFromId()] != -1)
					ev_fg->addTransition(from_id, id_lookup[t->getToId()], t->getNLSmiles());
				t_added_flags[i] = 1;
			}
		}
	}
	ev_graph_computed = 1;
}

void MolData::computeFragmentEvidenceValues(std::vector<double> &evidence, int frag_idx, const beliefs_t *beliefs ){

	const tmap_t *tomap = fg->getToIdTMap();

	evidence.resize(cfg->dv_spectrum_depths.size());
	for( int energy = 0; energy < cfg->dv_spectrum_depths.size(); energy++ ){
	
		//Find out the depth of the spectrum of interest in the beliefs
		int depth = cfg->dv_spectrum_depths[energy] - 1;
		if(cfg->use_single_energy_cfm )
			for( int i = 0; i < energy;  i++ ) 
				depth += cfg->dv_spectrum_depths[i];
		
		//Compute the accumulated belief for the fragment of interest at the depth of interest
		evidence[energy] = beliefs->ps[frag_idx][depth];
		std::vector<int>::const_iterator it = (*tomap)[frag_idx].begin();
		for( ; it != (*tomap)[frag_idx].end(); ++it )
			evidence[energy] = logAdd( evidence[energy], beliefs->tn[*it][depth]);
	}

}

void MolData::annotatePeaks( double abs_tol, double ppm_tol, bool prune_deadends){

	if(!ev_graph_computed){
		std::cout  << "Warning: called annotate peaks without first computing evidence graph" << std::endl;
		return;
	}

	unsigned int num_fragments = ev_fg->getNumFragments();
	std::vector<int> frag_flags(num_fragments,0);

	//Assign fragment annotations to peaks
	for( int energy = 0; energy < spectra.size(); energy++ ){
		Spectrum::const_iterator it = spectra[energy].begin();
		for( int pkidx = 0; it != spectra[energy].end();  ++it, pkidx++ ){
		
			double mass_tol = getMassTol( abs_tol, ppm_tol, it->mass );
			int num_fragments = ev_fg->getNumFragments();
			for( int fidx = 0; fidx < num_fragments; fidx++ ){
			
				const EvidenceFragment *f = ev_fg->getFragmentAtIdx(fidx);
				if( ev_fg->hasIsotopesIncluded() ){
					const Spectrum *isotope_spec = f->getIsotopeSpectrum();
					Spectrum::const_iterator itp = isotope_spec->begin();
					for( ; itp != isotope_spec->end(); ++itp ){
						if( fabs( itp->mass - it->mass ) <= mass_tol ){
							spectra[energy].addPeakAnnotation( pkidx, annotation_t(f->getId(),exp(f->getEvidence(energy))) );
							frag_flags[fidx] = 1;
						}
					}

				}
				else if( fabs( f->getMass() - it->mass ) <= mass_tol ){
					spectra[energy].addPeakAnnotation( pkidx, annotation_t(f->getId(),exp(f->getEvidence(energy))) );
					frag_flags[fidx] = 1;
				}
			}
		}
	}

	//Sort the annotations by evidence score
	for( int energy = 0; energy < spectra.size(); energy++ )
		spectra[energy].sortAndNormalizeAnnotations();

	if(!prune_deadends) return;

	std::vector<int> delete_flags(num_fragments);
	for(unsigned int i = 0; i < num_fragments; i++ ){	
		//If there is no peak for this fragment, check there is
		//an annotated descendent, else flag for deletion
		std::vector<int> direct_flags(num_fragments);
		ev_fg->setFlagsForDirectPaths( direct_flags, 0, frag_flags );
		if( ev_fg->fragmentIsRedundant(i, frag_flags, direct_flags) )
			delete_flags[i] = 1;
		else delete_flags[i] = 0;
	}

	//Rebuild the graph again without deleted fragments
	Transition null_t;
	EvidenceFragmentGraph *new_ev_fg = new EvidenceFragmentGraph( cfg );
	std::vector<int> id_lookup(num_fragments, -1);
	id_lookup[0] = new_ev_fg->addToGraphDirectNoCheck( *ev_fg->getFragmentAtIdx(0), &null_t, -1);
	for(unsigned int i = 0; i < ev_fg->getNumTransitions(); i++ ){	
		const Transition *t = ev_fg->getTransitionAtIdx(i);
		int fromid = t->getFromId();
		int toid = t->getToId();
		if( delete_flags[toid] || delete_flags[fromid] ) continue;
		if( id_lookup[toid] == -1 ) 
			id_lookup[toid] = new_ev_fg->addToGraphDirectNoCheck(*ev_fg->getFragmentAtIdx(toid), t, id_lookup[fromid]);
		else
			new_ev_fg->addTransition(id_lookup[fromid], id_lookup[toid], t->getNLSmiles() );
	}
	delete ev_fg;
	ev_fg = new_ev_fg;

	//Update annotations with new ids
	for( int energy = 0; energy < spectra.size(); energy++ ){
		Spectrum::const_iterator it = spectra[energy].begin();
		for( int peak_idx = 0; it != spectra[energy].end();  ++it, peak_idx++ ){
			for( unsigned int i = 0; i < it->annotations.size(); i++)
				spectra[energy].updateAnnotationId( peak_idx, i, id_lookup[it->annotations[i].first] );
		}
	}

}



void MolData::computeFeatureVectors( FeatureCalculator *fc, bool deleteMols ){

	fvs.resize( fg->getNumTransitions() );
	for( unsigned int i = 0; i < fg->getNumTransitions(); i++ ){
		const Transition *t = fg->getTransitionAtIdx(i);
		fvs[i] = fc->computeFV( t->getIon(), t->getNeutralLoss() );
		if( deleteMols ) fg->deleteMolsForTransitionAtIdx(i);
	}
	if( !deleteMols ) return;
	fg->clearAllSmiles();

}

void MolData::computeTransitionThetas( Param &param ){

	const unsigned int num_levels = param.getNumEnergyLevels();
	thetas.resize(num_levels);
	for( unsigned int energy = 0; energy < num_levels; energy++ ){
		
		//Compute the theta value for each feature vector
		thetas[energy].resize( fg->getNumTransitions() );
		for( unsigned int i = 0; i < fg->getNumTransitions(); i++ )
			thetas[energy][i] = param.computeTheta( *fvs[i], energy );		
	}
}

void MolData::computeTransitionProbabilities(){

	log_probs.resize( thetas.size() );
	for( unsigned int energy = 0; energy < thetas.size(); energy++ ){
		
		//Will store an entry for all transitions, followed by all persistences
		log_probs[energy].resize( fg->getNumTransitions() + fg->getNumFragments() );

		//Compute all the denominators
		std::vector<double> denom_cache( fg->getNumFragments() );
		const tmap_t *from_id_map = fg->getFromIdTMap();
		for( unsigned int i = 0; i < fg->getNumFragments(); i++ ){
		
			double denom = 0.0;
			std::vector<int>::const_iterator it = (*from_id_map)[i].begin();
			for( ; it != (*from_id_map)[i].end(); ++it )
				denom = logAdd( denom, thetas[energy][*it] );
			denom_cache[i] = denom;
		}
	
		//Set the transition log probabilities
		for( unsigned int i = 0; i < fg->getNumTransitions(); i++ ){
			const Transition *t = fg->getTransitionAtIdx(i);
			log_probs[energy][i] = thetas[energy][i] - denom_cache[ t->getFromId() ];
		}

		//Set the persistence log probabilities
		int offset = fg->getNumTransitions();
		for( unsigned int i = 0; i < fg->getNumFragments(); i++ )
			log_probs[energy][offset + i] = -denom_cache[i];

	}
	return;
}

void MolData::readInSpectraFromFile( const std::string &peak_filename, bool readToPredicted ){

	//Set the spectra that we are writing to
	std::vector<Spectrum> *spec_dest = &spectra;
	if( readToPredicted ) spec_dest = &predicted_spectra;

	//Clear any existing entries
	spec_dest->clear();
	spec_dest->resize(0);

	Spectrum *curr_spec;
	std::string line;
	std::ifstream ifs ( peak_filename.c_str() , std::ifstream::in );
	if(!ifs) std::cout << "Warning: Could not open file " << peak_filename << std::endl;
	int first = 1;
	while( ifs.good() ){
		getline( ifs, line );
		if( line.size() < 3 ) continue;
		//Check for the energy specifier - start a new spectrum if found 
		//or start one anyway if there is no energy specifier
		if( line.substr(0,3) == "low" 
			|| line.substr(0,3) == "med"
			|| line.substr(0,4) == "high" 
			|| line.substr(0,6) == "energy" ){
				spec_dest->push_back( Spectrum() );
				curr_spec = &(spec_dest->back());
				first = 0;
				continue;
		}
		else if( first ){
			spec_dest->push_back( Spectrum() );
			curr_spec = &(spec_dest->back());	
		}
		first = 0;

		//Otherwise allocate peak to the current spectrum
		std::stringstream ss(line);
		double mass, intensity;
		ss >> mass >> intensity;
		curr_spec->push_back( Peak(mass, intensity) );
	}
	std::vector<Spectrum>::iterator it = spec_dest->begin();
	for( ; it != spec_dest->end(); ++it ) it->normalizeAndSort();
	ifs.close();
}

void MolData::readInSpectraFromMSP( MspReader &msp, bool readToPredicted ){
	
	//Set the spectra that we are writing to
	std::vector<Spectrum> *spec_dest = &spectra;
	if( readToPredicted ) spec_dest = &predicted_spectra;

	//Clear any existing entries
	spec_dest->clear();
	spec_dest->resize(0);

	//Copy the spectra in
	const std::vector<Spectrum> *spec = msp.fetchSpectrumForId( id.c_str() );
	for( int energy = 0; energy < spec->size(); energy++ ){
		spec_dest->push_back( (*spec)[energy] );
		(*spec_dest)[energy].normalizeAndSort();
	}
}

std::vector<Spectrum>* MolData::readInSpectraFromMSPFileStream( std::istream &msp, int num_peaks, bool readToPredicted ){
	
	//Set the spectra that we are writing to
	std::vector<Spectrum> *spec_dest = &spectra;
	if( readToPredicted ) spec_dest = &predicted_spectra;

	//Clear any existing entries
	spec_dest->clear();
	spec_dest->resize(0);

	//Copy the spectra in
	Spectrum spec;
	std::string line;
	std::vector<std::string> split;
	for (int i = 0; i < num_peaks; ++i)
	{
		getline(msp, line);
		boost::split(split, line, boost::is_any_of(" "));
		spec.push_back(Peak(std::stod(split[0]), std::stod(split[1])));
	}
	spec_dest->push_back(spec);
	spec_dest->back().normalizeAndSort();
	return spec_dest;
}

void MolData::cleanSpectra(double abs_tol, double ppm_tol){

	std::vector<Spectrum>::iterator it = spectra.begin();
	for( ; it != spectra.end(); ++it )
		it->clean(abs_tol, ppm_tol);

}

void MolData::removePeaksWithNoFragment( double abs_tol, double ppm_tol ){

	std::vector<double> all_masses;
	getEnumerationSpectraMasses( all_masses );

	std::vector<Spectrum>::iterator it = spectra.begin();
	for( ; it != spectra.end(); ++it )
		it->removePeaksWithNoFragment( all_masses, abs_tol, ppm_tol );
}
	
bool MolData::hasEmptySpectrum() const{

	bool result = false;
	std::vector<Spectrum>::const_iterator it = spectra.begin();
	for( ; it != spectra.end(); ++it ){
		if( it->size() == 0 ) result = true;
	}
	return result;
}

void MolData::computePredictedSpectra( Param &param, bool postprocess, bool use_existing_thetas ){

	//Divert to other function if doing single energy CFM
	if( cfg->use_single_energy_cfm ){ 
		computePredictedSingleEnergySpectra( param, postprocess, use_existing_thetas );
		return;
	}

	//Compute the transition probabilities using this parameter set
	if(!use_existing_thetas) computeTransitionThetas( param );
	computeTransitionProbabilities();

	//Run forward inference
	std::vector<Message> msgs;
	Inference infer(this, cfg);
	infer.runInferenceDownwardPass( msgs, cfg->model_depth );

	//Generate and collect the peak results
	std::string outfilename;
	predicted_spectra.resize( cfg->spectrum_depths.size() );
	for( unsigned int energy = 0; energy < cfg->spectrum_depths.size(); energy++ ){
		predicted_spectra[energy].clear();

		//Extract the peaks from the relevant message
		int msg_depth = cfg->spectrum_depths[energy]-1;		
		Message *msg = &(msgs[msg_depth]);
		if( cfg->include_isotopes )
			translatePeaksFromMsgToSpectraWithIsotopes( predicted_spectra[energy], msg );
		else 
			translatePeaksFromMsgToSpectra( predicted_spectra[energy], msg );

	}

	if( postprocess ) postprocessPredictedSpectra();
	else{
		for( unsigned int energy = 0; energy < cfg->spectrum_depths.size(); energy++ ){
			predicted_spectra[energy].normalizeAndSort();
			predicted_spectra[energy].quantisePeaksByMass(10);
		}
	}
}

void MolData::computePredictedSingleEnergySpectra( Param &param, bool postprocess, bool use_existing_thetas ){

	//Compute the transition probabilities using this parameter set
	if( !use_existing_thetas ) computeTransitionThetas( param );
	computeTransitionProbabilities();

	//Generate and collect the peak results
	std::string outfilename;
	predicted_spectra.resize( cfg->spectrum_depths.size() );
	for( unsigned int energy = 0; energy < cfg->spectrum_depths.size(); energy++ ){
		predicted_spectra[energy].clear();

		config_t se_cfg;
		initSingleEnergyConfig( se_cfg, *cfg, energy );

		//Run forward inference
		std::vector<Message> msgs;
		Inference infer(this, &se_cfg);
		infer.runInferenceDownwardPass( msgs, se_cfg.model_depth );

		//Extract the peaks from the relevant message
		int msg_depth = se_cfg.spectrum_depths[0]-1;		
		Message *msg = &(msgs[msg_depth]);
		if( cfg->include_isotopes )
			translatePeaksFromMsgToSpectraWithIsotopes( predicted_spectra[energy], msg );
		else 
			translatePeaksFromMsgToSpectra( predicted_spectra[energy], msg );
	}

	if( postprocess ) postprocessPredictedSpectra();
	else{
		for( unsigned int energy = 0; energy < cfg->spectrum_depths.size(); energy++ ){
			predicted_spectra[energy].normalizeAndSort();
			predicted_spectra[energy].quantisePeaksByMass(10);
		}
	}
}

void MolData::translatePeaksFromMsgToSpectra( Spectrum &out_spec, Message *msg ){

	//Create the peaks
	std::map<double, Peak> peak_probs;
	Message::const_iterator itt = msg->begin();
	for( ; itt != msg->end() ; ++itt ){
		double mass = fg->getFragmentAtIdx(itt.index())->getMass();
		double intensity_contrib = exp(*itt);
		if( peak_probs.find(mass) != peak_probs.end() )
			peak_probs[mass].intensity += intensity_contrib;
		else
			peak_probs[mass].intensity = intensity_contrib;
		peak_probs[mass].mass = mass;
		peak_probs[mass].annotations.push_back( annotation_t(itt.index(),intensity_contrib) );
	}		
		
	//Add the peaks to the spectra
	std::map<double, Peak>::iterator itm = peak_probs.begin();
	for( ; itm != peak_probs.end(); ++itm )
		out_spec.push_back( itm->second );

}

void MolData::translatePeaksFromMsgToSpectraWithIsotopes( Spectrum &out_spec, Message *msg ){

	//Create the peaks
	std::map<double, Peak> peak_probs;
	Message::const_iterator itt = msg->begin();
	for( ; itt != msg->end() ; ++itt ){
		const Spectrum *isotope_spec = fg->getFragmentAtIdx(itt.index())->getIsotopeSpectrum();
		Spectrum::const_iterator itp = isotope_spec->begin();
		for( ; itp != isotope_spec->end(); ++itp ){
			double intensity_contrib = 0.01*itp->intensity * exp(*itt );
			double mass = itp->mass;	
			if( peak_probs.find(mass) != peak_probs.end() )
				peak_probs[mass].intensity += intensity_contrib;
			else
				peak_probs[mass].intensity = intensity_contrib;
			peak_probs[mass].mass = mass;
			peak_probs[mass].annotations.push_back( annotation_t(itt.index(),intensity_contrib) );
		}
	}		
		
	//Add the peaks to the spectra
	std::map<double, Peak>::iterator itm = peak_probs.begin();
	for( ; itm != peak_probs.end(); ++itm )
		out_spec.push_back( itm->second );

}

void MolData::writePredictedSpectraToFile( std::string &filename ){

	std::ofstream of;
	of.open(filename.c_str());
	if( !of.is_open() ){
		std::cout << "Warning: Trouble opening predicted spectrum file" << std::endl;
		return;
	}
	std::streambuf *buf = of.rdbuf();	
	std::ostream out(buf);
	outputSpectra( out, "Predicted" );
	of.close();
}

void MolData::writePredictedSpectraToMspFileStream( std::ostream &out ){
	for( unsigned int energy = 0; energy < getNumPredictedSpectra(); energy++ ){
		const Spectrum *spec = getPredictedSpectrum( energy );
		spec->outputToMspStream( out, id, cfg->ionization_mode, energy);
	}
}

void MolData::writePredictedSpectraToMgfFileStream( std::ostream &out ){
	double mw = getMolecularWeight();
	for( unsigned int energy = 0; energy < getNumPredictedSpectra(); energy++ ){
		const Spectrum *spec = getPredictedSpectrum( energy );
		spec->outputToMgfStream( out, id, cfg->ionization_mode, energy, mw);
	}
}

void MolData::writeFullEnumerationSpectrumToFile( std::string &filename ){
	
	std::ofstream of;
	of.open(filename.c_str());
	if( !of.is_open() ){
		std::cout << "Warning: Trouble opening enumerated spectrum file" << std::endl;
		return;
	}
	std::streambuf *buf = of.rdbuf();	
	std::ostream out(buf);
	std::vector<double> all_masses;
	getEnumerationSpectraMasses( all_masses );

	//All peaks of uniform intensity at each possible fragment mass
	for( unsigned int energy = 0; energy < spectra.size(); energy++ ){
		
		out << "energy" << energy << std::endl;

		out << std::setprecision(10);
		double peak_height = 100.0/(double)all_masses.size();
		for( int i = 0; i < all_masses.size(); i++ )
			out << all_masses[i] << " " << peak_height << std::endl;
	}

	of.close();
}

void MolData::writeFullEnumerationSpectrumToMspFileStream( std::ostream &out ){

	std::vector<double> all_masses;
	getEnumerationSpectraMasses( all_masses );

	out << "ID: " << id << std::endl;
	out << "Num peaks: " << all_masses.size() << std::endl; 
	out << std::setprecision(10);
	double peak_height = 100.0/(double)all_masses.size();
	for( int i = 0; i < all_masses.size(); i++ )
		out << all_masses[i] << " " << peak_height << std::endl;	
	out << std::endl;

}

void MolData::getEnumerationSpectraMasses( std::vector<double> &output_masses ){

	//Get and sort the fragment masses
	unsigned int numf = fg->getNumFragments();
	std::vector<double> all_masses;
	if( fg->hasIsotopesIncluded() ){
		for( unsigned int i = 0; i < numf; i++ ){
			const Fragment *f = fg->getFragmentAtIdx(i);
			const Spectrum *isotope_spec = f->getIsotopeSpectrum();
			Spectrum::const_iterator itp = isotope_spec->begin();
			for( ; itp != isotope_spec->end(); ++itp )
				all_masses.push_back( itp->mass );
		}		
	}
	else{
		all_masses.resize(numf);
		for( unsigned int i = 0; i < numf; i++ ){
			const Fragment *f = fg->getFragmentAtIdx(i);
			all_masses[i] = f->getMass();
		}
	}
	std::sort( all_masses.begin(), all_masses.end() );

	//Remove repeats
	output_masses.resize(all_masses.size());
	int i = 0;
	double prev_mass = -1.0, tol = 1e-10;
	std::vector<double>::iterator it = all_masses.begin();
	for( ; it != all_masses.end(); ++it ){ 
		if( fabs(*it - prev_mass) > tol )
			output_masses[i++] = *it;
		prev_mass = *it;
	}
	int num_unique = i;
	output_masses.resize( num_unique );

}

void MolData::outputSpectra( std::ostream &out, const char*spec_type, bool do_annotate ){
	
	std::vector<Spectrum> *spectra_to_output;
	if(std::string(spec_type) == "Predicted") spectra_to_output = &predicted_spectra;
	else if(std::string(spec_type) == "Experimental") spectra_to_output = &spectra;
	else std::cout << "Unknown spectrum type to output: " << spec_type << std::endl;

	std::vector<Spectrum>::iterator it = spectra_to_output->begin();
	for( int energy = 0; it != spectra_to_output->end(); ++it, energy++ ){
		out << "energy" << energy << std::endl;
		it->outputToStream( out, do_annotate );
	}
}

void MolData::postprocessPredictedSpectra(double perc_thresh, int min_peaks, int max_peaks){
	
	std::vector<Spectrum>::iterator it = predicted_spectra.begin();
	for( ; it != predicted_spectra.end(); ++it){
		it->quantisePeaksByMass(10);	
		it->postProcess( perc_thresh, min_peaks, max_peaks );
		it->normalizeAndSort();
		it->sortAndNormalizeAnnotations();
	}
}

void MolData::quantisePredictedSpectra(int num_dec_places){
	
	std::vector<Spectrum>::iterator it = predicted_spectra.begin();
	for( ; it != predicted_spectra.end(); ++it)
		it->quantisePeaksByMass(num_dec_places);
}

void MolData::quantiseMeasuredSpectra(int num_dec_places){
	
	std::vector<Spectrum>::iterator it = spectra.begin();
	for( ; it != spectra.end(); ++it)
		it->quantisePeaksByMass(num_dec_places);
}

MolData::~MolData(){

	if( graph_computed ) delete fg;
	if( ev_graph_computed ) delete ev_fg;

	//Delete any computed feature vectors
	std::vector<FeatureVector *>::iterator it = fvs.begin();
	for( ; it != fvs.end(); ++it ) delete *it;

}
