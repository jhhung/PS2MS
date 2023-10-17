/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraphGenerator.cpp
#
# Description: 	FragmentGraphGenerator class for generating a fragment tree.
#				Also contains Break and FragmentTreeNode classes, for use in
#				this process.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentGraphGenerator.h"
#include "MILP.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <INCHI-API/inchi.h>




//Start a graph. Compute can then add to this graph, but it is the caller's 
//responsibility to delete it
FragmentGraph *FragmentGraphGenerator::createNewGraph(config_t *cfg){
	current_graph = new FragmentGraph(cfg);
	id_depth_computed_cache.clear();	//The graph is empty, so clear all computation records
	return current_graph;
}

//Start a graph. Compute can then add to this graph, but it is the caller's 
//responsibility to delete it
FragmentGraph *LikelyFragmentGraphGenerator::createNewGraph( config_t *cfg ){
	current_graph = new FragmentGraph(cfg);
	id_prob_computed_cache.clear();	//The graph is empty, so clear all computation records
	return current_graph;
}

//Create the starting node from a smiles or inchi string - responsibility of caller to delete
FragmentTreeNode *FragmentGraphGenerator::createStartNode( std::string &smiles_or_inchi, int ionization_mode, RDKit::RWMol* rwmol){
	
	//Create the RDKit mol - this will be the ion
	if (!rwmol)
	{
		if( smiles_or_inchi.substr(0,6) == "InChI=" ){
			RDKit::ExtraInchiReturnValues rv;
			rwmol = RDKit::InchiToMol( smiles_or_inchi, rv);
		}else
			rwmol = RDKit::SmilesToMol( smiles_or_inchi ); 
	}

	//This is dirty, but for some reason RDKit doesn't throw the exception...
	if(!rwmol) throw RDKit::SmilesParseException("Error occurred - assuming Smiles Parse  Exception");
	
	//Remove stereochemistry
	RDKit::MolOps::removeStereochemistry( *rwmol );
	
	//Compute and label anything required by features that won't be present once the molecule breaks
	fh->addLabels( rwmol );

	//Kekulize and count the excess electron pairs
	RDKit::MolOps::Kekulize( *rwmol );
	std::vector<int> e_loc(rwmol->getNumAtoms());
	int num_ep = countExtraElectronPairs( rwmol, e_loc );

	//Check for multiple fragments - flag ionic atoms
	std::vector<int> mapping;
	int num_frags = RDKit::MolOps::getMolFrags( *rwmol, mapping );

	//Initialize some properties of the molecule
	RDKit::MolOps::findSSSR( *rwmol );	
	RDKit::RingInfo *rinfo = rwmol->getRingInfo();
	RDKit::ROMol::AtomIterator ai;

	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){ 
		(*ai)->setProp( "FragIdx", 0 );
		(*ai)->setProp( "NumUnbrokenRings", rinfo->numAtomRings((*ai)->getIdx()) );
	}
	int num_ionic = addIonicChargeLabels(rwmol);
	if( num_frags - num_ionic != 1 ){
		std::cout << "Unsupported input molecule: Too many starting fragments in " << smiles_or_inchi << std::endl;
		throw FragmentGraphGenerationException();
	}

	//Ionize the molecule
	applyIonization( rwmol, ionization_mode );

	return( new FragmentTreeNode(romol_ptr_t(rwmol), num_ep, 0, fh, e_loc ) );
}

int FragmentGraphGenerator::countExtraElectronPairs( RDKit::RWMol *rwmol, std::vector<int> &output_e_loc ){

	//Compute the total number of bond electrons
	double total_bond_es = 0;
	RDKit::ROMol::AtomIterator ai;
	for( ai = rwmol->beginAtoms(); ai != rwmol->endAtoms(); ++ai ){
	
		int Hs = (*ai)->getTotalNumHs();
		int implied_val = (*ai)->getExplicitValence() + (*ai)->getImplicitValence();
		total_bond_es += (implied_val - Hs);
		output_e_loc[(*ai)->getIdx()] = implied_val - Hs - (*ai)->getDegree();	//Record where the extra electrons were
	}

	//Adjust to figure out how many are extras (non-single bond electrons)
	//and divide by 2 to count pairs
	double extra_bond_eps = (total_bond_es - rwmol->getNumBonds()*2)/2;
	return (int)extra_bond_eps;

}

//Helper function - check if the fragment has already been computed to at least this depth
int FragmentGraphGenerator::alreadyComputed(int id, int remaining_depth){
	if( id_depth_computed_cache.find(id) == id_depth_computed_cache.end()	//Not found
		|| id_depth_computed_cache[id] < remaining_depth){					//Or computed previously to less depth
		id_depth_computed_cache[id] = remaining_depth;
		return 0;
	}
	return 1;
}

//Compute a FragmentGraph starting at the given node and computing to the depth given.
//The output will be appended to the current_graph
void FragmentGraphGenerator::compute( FragmentTreeNode &node, int remaining_depth, int parentid, int remaining_ring_breaks ){
	
	if( current_graph->getNumFragments() > MAX_FRAGMENTS_PER_MOLECULE 
		|| current_graph->getNumTransitions() > MAX_TRANSITIONS_PER_MOLECULE ){
		std::cout << "Maximum number of fragments or transitions exceeded." << std::endl;
		throw FragmentGraphMaxSizeExceededException();
	}

	if(verbose) std::cout << current_graph->getNumFragments() << ":" << RDKit::MolToSmiles( *node.ion.get() ) << std::endl;

	//Add the node to the graph, and return a fragment id
	int id = -1;
	if( mols_to_fv ) 
		id = current_graph->addToGraphAndReplaceMolWithFV( node, parentid, fc );
	else 
		id = current_graph->addToGraph( node, parentid );

	//Only compute to the desired depth
	if( remaining_depth <= 0 ) return;

	//If the node was already in the graph at sufficient depth, skip any further computation
	if( alreadyComputed(id, remaining_depth) ){ 
		if(verbose)	std::cout << "Node already computed: Skipping" << std::endl;
		return;
	}

	//Generate Breaks
	std::vector<Break> breaks;
	bool h_loss_allowed = false;
	if( parentid < 0 )	//Break from Precursor
		h_loss_allowed = current_graph->includesHLossesPrecursorOnly() || current_graph->includesHLosses();
	else				//Break from Non-Precursor
		h_loss_allowed = !(current_graph->includesHLossesPrecursorOnly()) && current_graph->includesHLosses();
	node.generateBreaks( breaks, h_loss_allowed );
		
	//Iterate over the possible breaks
	std::vector<Break>::iterator it = breaks.begin();
	for( ; it != breaks.end(); ++it ){

		if( it->isRingBreak() && remaining_ring_breaks == 0 ) continue;

		for( int ifrag_idx = 0; ifrag_idx < it->getNumIonicFragAllocations(); ifrag_idx++ ){
			
			node.applyBreak( *it, ifrag_idx );
			node.generateChildrenOfBreak( *it );

			int child_remaining_ring_breaks = remaining_ring_breaks;
			if( it->isRingBreak() && remaining_ring_breaks > 0 ) child_remaining_ring_breaks--;

			//Recur over children
			std::vector<FragmentTreeNode>::iterator itt = node.children.begin();
			for( ; itt != node.children.end(); ++itt ){
				compute( *itt, remaining_depth-1, id, child_remaining_ring_breaks );
			}

			//Undo and remove children
			node.undoBreak( *it, ifrag_idx );
		}
		node.children = std::vector<FragmentTreeNode>();
	}

}


//Compute a FragmentGraph starting at the given node and computing to the depth given.
//The output will be appended to the current_graph
void LikelyFragmentGraphGenerator::compute( FragmentTreeNode &node,  int remaining_depth, int parentid, double parent_log_prob, int remaining_ring_breaks ){
	
	//Check Timeout
	if( parentid < 0 ) start_time = time( NULL );
	time_t current_time = time(NULL);
	if( cfg->fragraph_compute_timeout_in_secs > 0 ){
		if( (current_time - start_time ) > cfg->fragraph_compute_timeout_in_secs )
			throw FragmentGraphTimeoutException();
	}



	//Add the node to the graph, and return a fragment id: note, no mols or fv will be set,
	//but the precomputed theta value will be used instead
	int id = -1;
	id = current_graph->addToGraphWithThetas( node, node.getAllTmpThetas(), parentid );

	//Reached max depth? 
	if( remaining_depth <= 0 ) return;

	//If we've already run the fragmentation on this fragment with an equal or higher
	//probability offset, we don't need to run again, unless it is persisting
	if( alreadyComputed( id, parent_log_prob )) return;

	//Generate Children 
	std::vector<Break> breaks;
	bool h_loss_allowed = false;
	if( parentid < 0 )	//Break from Precursor
		h_loss_allowed = current_graph->includesHLossesPrecursorOnly() || current_graph->includesHLosses();
	else				//Break from Non-Precursor
		h_loss_allowed = !(current_graph->includesHLossesPrecursorOnly()) && current_graph->includesHLosses();
	node.generateBreaks( breaks, h_loss_allowed );
	std::vector<Break>::iterator it = breaks.begin();
	std::vector<int> children_isring;
	for( ; it != breaks.end(); ++it ){
		
		//Record the index where the children for this break start (if there are any)
		if( remaining_ring_breaks == 0 && it->isRingBreak() ) continue;

		int start_idx = node.children.size();

		//Generate the children
		for( int iidx = 0; iidx < it->getNumIonicFragAllocations(); iidx++ ){
			node.applyBreak( *it, iidx );
			node.generateChildrenOfBreak( *it );
			node.undoBreak( *it, iidx);
		}

		//Record whether each child was generated from a ring break
		for( int i = start_idx; i < node.children.size(); i++ ) 
			children_isring.push_back(it->isRingBreak());
	} 

	//Compute child thetas
	std::vector<FragmentTreeNode>::iterator itt = node.children.begin();
	for( ; itt != node.children.end(); ++itt ){ 
		Transition tmp_t( -1, -1, itt->nl, itt->ion ); 
		FeatureVector *fv = fc->computeFV( tmp_t.getIon(), tmp_t.getNeutralLoss() );
		for( int engy = cfg->spectrum_depths.size()-1; engy >= 0; engy-- ){
			if( is_nn_params ) itt->setTmpTheta( nnparam->computeTheta(*fv, engy), engy );
			else itt->setTmpTheta( param->computeTheta(*fv, engy), engy );
		}
		delete fv;
	}

	//Compute child probabilities (including persistence) - for all energy levels
	std::vector<double> denom( cfg->spectrum_depths.size() );
	for(int i = 0; i < denom.size(); i++ ) denom[i] = 0.0;
	for( itt = node.children.begin(); itt != node.children.end(); ++itt ){
		for(int energy = 0; energy < denom.size(); energy++ )
			denom[energy] = logAdd( denom[energy], itt->getTmpTheta(energy) );
	}

	//Add and recur over likely children if above threshold for any energy level
	double max_child_prob; int child_idx = 0;
	for( itt = node.children.begin(); itt != node.children.end(); ++itt, child_idx++ ){
		max_child_prob = log_prob_thresh - 10.0;
		for(int energy = 0; energy < denom.size(); energy++ ){		
			double child_log_prob = itt->getTmpTheta(energy) - denom[energy] + parent_log_prob;
			if( child_log_prob > max_child_prob ) max_child_prob = child_log_prob;
		}
		if( max_child_prob >= log_prob_thresh ){
			int child_remaining_ring_breaks = remaining_ring_breaks;
			if( children_isring[child_idx] && child_remaining_ring_breaks > 0 ) 
				child_remaining_ring_breaks--;
			compute( *itt, remaining_depth-1, id, max_child_prob, child_remaining_ring_breaks );
		}
	}

	//Clear the children
	node.children = std::vector<FragmentTreeNode>();
}

void FragmentGraphGenerator::applyIonization( RDKit::RWMol *rwmol, int ionization_mode ){
	
	int rad_side = -1; if( ionization_mode == POSITIVE_EI_IONIZATION_MODE ) rad_side = 0;
	bool is_neg = (ionization_mode == NEGATIVE_ESI_IONIZATION_MODE );

	boost::tuple<int,int,int> pindx_nidx_ridx(-1, -1, -1);
	boost::tuple<bool,bool,bool> alreadyq_oktogo = FragmentTreeNode::findAlreadyChargedOrSplitCharge(pindx_nidx_ridx, *rwmol, 0, rad_side, is_neg );
	if( boost::get<0>(alreadyq_oktogo) && !boost::get<1>(alreadyq_oktogo)){
		std::cout << "Could not ionize - already charged molecule and didn't know what to do here" << std::endl;
		throw IonizationException();	
	}
	else if( !boost::get<0>(alreadyq_oktogo)){
		std::pair<int,int> qidx_ridx = FragmentTreeNode::findChargeLocation( *rwmol, 0, rad_side, is_neg );
		if( qidx_ridx.first < 0 ){
			std::cout << "Could not ionize - no location found for charge" << std::endl;
			throw IonizationException();
		}
		else if( qidx_ridx.second < 0 && rad_side >= 0 ){
			std::cout << "Could not ionize - no location found for radical" << std::endl;
			throw IonizationException();		
		}
		try{
			FragmentTreeNode::assignChargeAndRadical( *rwmol, qidx_ridx.first, qidx_ridx.second, is_neg  );
			RDKit::MolOps::sanitizeMol(*rwmol);	//Re-sanitize...sometimes RDKit only throws the exception the second time...
		}
		catch( RDKit::MolSanitizeException e ){
			std::cout << "Could not ionize - sanitization failure" << std::endl;
			throw IonizationException();		
		}
	}

}


//Helper function - check if the fragment has already been computed with at least this probability offset
int LikelyFragmentGraphGenerator::alreadyComputed(int id, double prob_offset){
	if( id_prob_computed_cache.find(id) == id_prob_computed_cache.end()
		|| id_prob_computed_cache[id] < prob_offset){					//Or computed previously with lower prob
		id_prob_computed_cache[id] = prob_offset;
		return 0;
	}
	return 1;
}
