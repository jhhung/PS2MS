/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGraph.cpp
#
# Description: 	FragmentGraph class for holding the results of a generated
#				fragment graph.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentGraph.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/MolOps.h>
 
Transition::Transition( int a_from_id, int a_to_id, const romol_ptr_t &a_nl, const romol_ptr_t &an_ion ){

	RDKit::Atom *root = NULL, *other_root = NULL;
	RDKit::Atom *first_atom = an_ion.get()->getAtomWithIdx(0);
	if( first_atom->hasProp( "Root" ) && first_atom->hasProp( "OtherRoot" )  ){
		root = getLabeledAtom( an_ion, "Root" );
		other_root = getLabeledAtom( an_ion, "OtherRoot" );
	}
	else std::cout << "Warning: Ion Root atoms not defined" << std::endl;
	ion = RootedROMolPtr( an_ion, root, other_root );
	
	root = NULL; other_root = NULL;
	first_atom = a_nl.get()->getAtomWithIdx(0);
	if( first_atom->hasProp( "Root" ) && first_atom->hasProp( "OtherRoot" )  ){
		root = getLabeledAtom( a_nl, "Root" );
		other_root = getLabeledAtom( a_nl, "OtherRoot" );
	}
	else std::cout << "Warning: NL Root atoms not defined" << std::endl;
	nl = RootedROMolPtr( a_nl, root, other_root );
	
	to_id = a_to_id;
	from_id = a_from_id;
	nl_smiles = RDKit::MolToSmiles( *(nl.mol.get()) );
}

Transition::Transition( int a_from_id, int a_to_id, const RootedROMolPtr &a_nl, const RootedROMolPtr &an_ion  ) :
	to_id( a_to_id ), from_id( a_from_id ), nl( a_nl ), ion( an_ion ){
	nl_smiles = RDKit::MolToSmiles( *(nl.mol.get()) );
}

//Function to remove detour transitions from the graph (used if !cfg.allow_frag_detours)
void FragmentGraph::removeDetours(){

	std::vector<int> id_map(transitions.size());	//Map old id to new id (or -1 if deleting)

	//Remove detour transitions, and record a mapping of old->new transition ids
	int next_id = 0;
	for( int i = 0; i < transitions.size(); i++ ){
		Transition *t = &transitions[i];
		Fragment *f_parent = &fragments[t->getFromId()];
		Fragment *f_child = &fragments[t->getToId()];
	
		if( f_parent->getDepth() >= f_child->getDepth() )
			id_map[i] = -1;
		else{
			transitions[next_id] = transitions[i];
			id_map[i] = next_id++;
		}
	}
	transitions.resize(next_id);

	//Update the id maps
	for( int i = 0; i < fragments.size(); i++ ){
		
		//Update the To Id Tmap
		int count = 0;
		for( int j = 0; j < to_id_tmap[i].size(); j++ ){
			if( id_map[ to_id_tmap[i][j] ] >= 0 )
				to_id_tmap[i][count++] = id_map[ to_id_tmap[i][j] ];
		}
		to_id_tmap[i].resize(count);

		//Update the From Id Tmap
		count = 0;
		for( int j = 0; j < from_id_tmap[i].size(); j++ ){
			if( id_map[ from_id_tmap[i][j] ] >= 0 )
				from_id_tmap[i][count++] = id_map[ from_id_tmap[i][j] ];
		}
		from_id_tmap[i].resize(count);
	}
}


//Add a fragment node to the graph (should be the only way to modify the graph)
//	-- Add a fragment, or return an id, if it already exists
//  -- Add a transition to this fragment, based on the provided parent fragment id
//  -- Update the relevant tmaps
// Note: If parentid < 0 (i.e. starting node) doesn't add transition.
int FragmentGraph::addToGraph( const FragmentTreeNode &node, int parentid ){

	//If the fragment doesn't exist, add it
	double mass = getMonoIsotopicMass( node.ion );
	int id = addFragmentOrFetchExistingId( node.ion, mass );

	if( parentid < 0 || fragments[id].getDepth() == -1 ) fragments[id].setDepth( node.depth );	//Set start fragment depth

	//Add a transition, IF:
	// - This transition does not exist (i.e. this parent_id -> id)
	// - There is no (already added) shorter path to this fragment (we track the minimum depth at which each fragment is seen)
	//  (NOTE: This will not remove previously added longer paths, we will remove them all at once at the end)
	if( parentid >= 0 && findMatchingTransition( parentid, id ) < 0 && (allow_frag_detours || node.depth <= fragments[id].getDepth()) ){
		
		if( node.depth < fragments[id].getDepth() )
			fragments[id].setDepth( node.depth );	//Update the depth of the fragment

		int idx = transitions.size();
		transitions.push_back( Transition( parentid, id, node.nl, node.ion ) );
	
		//Update the tmaps
		from_id_tmap[parentid].push_back( idx );
		to_id_tmap[id].push_back( idx );
	}

	return id;

}

//As for previous function, but delete the mols in the transition and compute and store a feature vector instead
int FragmentGraph::addToGraphAndReplaceMolWithFV( const FragmentTreeNode &node, int parentid, FeatureCalculator *fc ){

	//If the fragment doesn't exist, add it
	double mass = getMonoIsotopicMass( node.ion );
	int id = addFragmentOrFetchExistingId( node.ion, mass );

	if( parentid < 0 || fragments[id].getDepth() == -1 ) fragments[id].setDepth( node.depth );	//Set start fragment depth

	//Add a transition, if one does not exist
	if( parentid >= 0 && findMatchingTransition( parentid, id ) < 0 && (allow_frag_detours || node.depth <= fragments[id].getDepth()) ){
		
		if( node.depth < fragments[id].getDepth() )
			fragments[id].setDepth( node.depth );	//Update the depth of the fragment	

		int idx = transitions.size();
		transitions.push_back( Transition( parentid, id, node.nl, node.ion ) );
	
		//Update the tmaps
		from_id_tmap[parentid].push_back( idx );
		to_id_tmap[id].push_back( idx );

		//Compute a feature vector
		Transition *t = &(transitions.back());
		FeatureVector *fv;
		try{
			fv = fc->computeFV( t->getIon(), t->getNeutralLoss() );
			t->setTmpFV(fv);
		}
		catch( FeatureCalculationException fe ){
			//If we couldn't compute the feature vector, set a dummy feature vector with bias only.
			fv = new FeatureVector();
			fv->addFeatureAtIdx(1.0,idx);
			fv->addFeatureAtIdx(0.0, fc->getNumFeatures()-1);
		}

		//Delete the mols
		t->deleteIon();
		t->deleteNeutralLoss();

	}

	return id;
}

int FragmentGraph::addToGraphWithThetas(const FragmentTreeNode &node, const std::vector<double> *thetas, int parentid ){

	//If the fragment doesn't exist, add it
	double mass = getMonoIsotopicMass( node.ion );
	int id = addFragmentOrFetchExistingId( node.ion, mass );

	if( parentid < 0 || fragments[id].getDepth() == -1 ) fragments[id].setDepth( node.depth );	//Set start fragment depth

	//Add a transition, if one does not exist
	if( parentid >= 0 && findMatchingTransition( parentid, id ) < 0 && (allow_frag_detours || node.depth <= fragments[id].getDepth()) ){
		
		int idx = transitions.size();
		transitions.push_back( Transition( parentid, id, node.nl, node.ion ) );
	
		//Update the tmaps
		from_id_tmap[parentid].push_back( idx );
		to_id_tmap[id].push_back( idx );

		//Set the theta values and delete the mols
		Transition *t = &(transitions.back());
		t->setTmpThetas( thetas );
		t->deleteIon();
		t->deleteNeutralLoss();
	}

	return id;

}

//Direct constructor that bipasses the mols altogether and directly sets the nl_smiles
int EvidenceFragmentGraph::addToGraphDirectNoCheck( const EvidenceFragment &fragment, const Transition *transition, int parentid ){

	int id = fragments.size();
	fragments.push_back( EvidenceFragment( fragment, id ) );
	from_id_tmap.resize(id+1);
	to_id_tmap.resize(id+1);
	if(parentid >= 0 ) addTransition( parentid, id, transition->getNLSmiles() );
	return id;
}

void EvidenceFragmentGraph::addTransition(int from_id, int to_id, const std::string *nl_smiles ){
	unsigned int idx = transitions.size();
	transitions.push_back( Transition( from_id, to_id, nl_smiles ) );
	from_id_tmap[from_id].push_back( idx );
	to_id_tmap[to_id].push_back( idx );
}

int FragmentGraph::addFragmentOrFetchExistingId( romol_ptr_t ion, double mass ){

	std::string reduced_smiles;
	
	//Round the mass to 5 decimal places and use that as an initial filter
	double rounded_mass = floor(mass*10000.0 + 0.5)/10000.0;	
	if( frag_mass_lookup.find( rounded_mass ) != frag_mass_lookup.end() ){
		
		//Create a copy of the ion and then reduce it, making all bonds single and filling in hydrogens
		RDKit::RWMol f1_copy = *ion.get();
		reduceMol( f1_copy );
		RDKit::MolOps::sanitizeMol(f1_copy);

		//Found an entry with this mass, check the linked fragments for a match
		std::vector<int>::iterator it;
		it = frag_mass_lookup[rounded_mass].begin();
		for( ; it != frag_mass_lookup[rounded_mass].end(); ++it ){
			try{
				RDKit::RWMol *f2_reduced = RDKit::SmilesToMol( *fragments[*it].getReducedSmiles() );
				if( areMatching( &f1_copy, f2_reduced ) ){ 
					delete f2_reduced;
					return *it;
				}
				delete f2_reduced;
			}
		    catch( RDKit::MolSanitizeException e ){
				std::cout << "Could not sanitize " << *fragments[*it].getReducedSmiles() << std::endl;
				throw e;
			}
		}
		reduced_smiles = RDKit::MolToSmiles(f1_copy);
	}
	else frag_mass_lookup[rounded_mass];
	
	//No match found, create the fragment
	if( reduced_smiles.size() == 0 ){
		RDKit::RWMol f1_copy = *ion.get();
		reduceMol( f1_copy );
		reduced_smiles = RDKit::MolToSmiles(f1_copy);
	}
	std::string smiles = RDKit::MolToSmiles(*ion.get());
	int newid = fragments.size();
	if( include_isotopes ){
		Spectrum isotope_spectrum;
		long charge = RDKit::MolOps::getFormalCharge( *ion.get() );
		isotope->computeIsotopeSpectrum( isotope_spectrum, ion, charge );
		fragments.push_back( Fragment( smiles, reduced_smiles, newid, mass, isotope_spectrum) );
	}
	else fragments.push_back( Fragment( smiles, reduced_smiles, newid, mass) );
	frag_mass_lookup[rounded_mass].push_back( newid );
	from_id_tmap.resize(newid+1);
	to_id_tmap.resize(newid+1);
	return newid;
}

bool FragmentGraph::areMatching( RDKit::ROMol *f1_reduced_ion, RDKit::ROMol *f2_reduced_ion ){

	//Quick preliminary check to throw away non-matches
	if( f1_reduced_ion->getNumAtoms() != f2_reduced_ion->getNumAtoms() )
		return false;

	//Note: this will fail if one mol is a substruct
	//of another but not an exact match, however given we've just
	//checked that the number of atoms match, this can't happen
	RDKit::MatchVectType match;
	return RDKit::SubstructMatch( *f1_reduced_ion, *f2_reduced_ion, match, false, false );
}

void FragmentGraph::reduceMol( RDKit::RWMol &rwmol ){

	//Ensure that the OrigValence tags are set, otherwise add them
	int israd = moleculeHasSingleRadical(&rwmol);
	RDKit::ROMol::AtomIterator ai; 
	if( !rwmol.getAtomWithIdx(0)->hasProp("OrigValence") ){
		for( ai = rwmol.beginAtoms(); ai != rwmol.endAtoms(); ++ai ){
			int origval = (*ai)->getExplicitValence() + (*ai)->getImplicitValence() + (*ai)->getNumRadicalElectrons() - (*ai)->getFormalCharge();
			//For radicals generated by removing a multibond, charge and radical are separate, so charge needs to be added, not subtracted
			if( israd && (*ai)->getFormalCharge() && !(*ai)->getNumRadicalElectrons() && (*ai)->getSymbol() == "C") origval += 2*(*ai)->getFormalCharge();
			(*ai)->setProp("OrigValence", origval );
		}
	}

	//Set all bonds to SINGLE
	RDKit::ROMol::BondIterator bi; 
	for( bi = rwmol.beginBonds(); bi != rwmol.endBonds(); ++bi ){
		(*bi)->setBondType( RDKit::Bond::SINGLE );
		(*bi)->setIsAromatic( false );
	}
	
	//Set all Hydrogens to give full valence and no charge or radicals
	for( ai = rwmol.beginAtoms(); ai != rwmol.endAtoms(); ++ai ){
		(*ai)->setFormalCharge(0);
		(*ai)->setNumRadicalElectrons(0);
		int valence;
		(*ai)->getProp("OrigValence", valence);
		int numH = valence - (*ai)->getDegree();
		if( numH < 0 ){ 
			(*ai)->setFormalCharge(1);
			numH = 0;
		}
		(*ai)->setNumExplicitHs(numH);
		(*ai)->setIsAromatic(false);
	}

}

//Find the id for an existing transition that matches the input ids
//or -1 in the case where no such transition is found
int FragmentGraph::findMatchingTransition( int from_id, int to_id ){
	std::vector<int>::iterator it;
	it = from_id_tmap[ from_id ].begin();
	for( ; it != from_id_tmap[ from_id ].end(); ++it )
		if( transitions[*it].getToId() == to_id ) return *it;
	return -1;
}


//Write the Fragments only to file (formerly the backtrack output - without extra details)
void FragmentGraph::writeFragmentsOnly( std::ostream &out ) const{

	std::vector<Fragment>::const_iterator it = fragments.begin();
	for( ; it != fragments.end(); ++it ){
		out << it->getId() << " ";
		out << std::setprecision(10) << it->getMass() << " ";
		out << *it->getIonSmiles() << std::endl;
	}
}

//Write the FragmentGraph to file (formerly the transition output - without feature details)
void FragmentGraph::writeFullGraph( std::ostream &out ) const{

	//Fragments
	out << fragments.size() << std::endl;
	writeFragmentsOnly( out );
	out << std::endl;

	//Transitions
	std::vector<Transition>::const_iterator itt = transitions.begin();
	for( ; itt != transitions.end(); ++itt ){
		out << itt->getFromId() << " ";
		out << itt->getToId() << " ";
		out << *itt->getNLSmiles() << std::endl;
	}
}

void FragmentGraph::writeFeatureVectorGraph( std::ostream &out, bool include_isotopes ) const{

	//Fragments
	unsigned int inlcis = include_isotopes;
	out.write(reinterpret_cast<const char *>(&inlcis), sizeof(inlcis));
	unsigned int numf = fragments.size();
	out.write(reinterpret_cast<const char *>(&numf), sizeof(numf));
	std::vector<Fragment>::const_iterator it = fragments.begin();
	for( ; it != fragments.end(); ++it ){
		int id = it->getId();
		out.write(reinterpret_cast<const char *>(&id), sizeof(id));
		double mass = it->getMass();
		out.write(reinterpret_cast<const char *>(&mass), sizeof(mass));
		if( include_isotopes ){
			const Spectrum *isospec = it->getIsotopeSpectrum();
			unsigned int iso_size = isospec->size();
			out.write(reinterpret_cast<const char *>(&iso_size), sizeof(iso_size));
			Spectrum::const_iterator itp = isospec->begin();
			for( ; itp != isospec->end(); ++itp ){
				out.write(reinterpret_cast<const char *>(&itp->mass), sizeof(itp->mass));
				out.write(reinterpret_cast<const char *>(&itp->intensity), sizeof(itp->intensity));
			}
		}
	}
	//Transitions
	unsigned int numt = transitions.size();
	out.write(reinterpret_cast<const char *>(&numt), sizeof(numt));
	std::vector<Transition>::const_iterator itt = transitions.begin();
	for( ; itt != transitions.end(); ++itt ){
		int fromid = itt->getFromId(); int toid = itt->getToId();
		out.write(reinterpret_cast<const char *>(&fromid), sizeof(fromid));
		out.write(reinterpret_cast<const char *>(&toid), sizeof(toid));
		FeatureVector *fv = itt->getTmpFV();
		unsigned int num_set =  fv->getNumSetFeatures();
		unsigned int fv_len = fv->getTotalLength();
		out.write(reinterpret_cast<const char *>(&num_set), sizeof(num_set));
		out.write(reinterpret_cast<const char *>(&fv_len), sizeof(fv_len));
		std::vector<feature_t>::const_iterator fvit = fv->getFeatureBegin();
		for( ; fvit != fv->getFeatureEnd(); ++fvit ){
			unsigned int fidx = *fvit;
			out.write(reinterpret_cast<const char *>(&fidx), sizeof(fidx));
		}
	}
}

void FragmentGraph::readFeatureVectorGraph( std::istream &ifs ){

	std::string null = "";
	unsigned int include_isotopes; ifs.read(reinterpret_cast<char *>(&include_isotopes), sizeof(include_isotopes));
	unsigned int numf; ifs.read(reinterpret_cast<char *>(&numf), sizeof(numf));
	for( int i = 0; i < numf; i++ ){
		int id; ifs.read(reinterpret_cast<char *>(&id), sizeof(id));
		double mass; ifs.read(reinterpret_cast<char *>(&mass), sizeof(mass));
		if( include_isotopes ){
			Spectrum isospec;
			unsigned int iso_size; ifs.read(reinterpret_cast<char *>(&iso_size), sizeof(iso_size));
			for( int j = 0; j < iso_size; j++ ){
				double pmass; ifs.read(reinterpret_cast<char *>(&pmass), sizeof(pmass));
				double pintensity; ifs.read(reinterpret_cast<char *>(&pintensity), sizeof(pintensity));
				isospec.push_back( Peak(pmass, pintensity) );
			}
			fragments.push_back( Fragment( null, null, id, mass, isospec ) );
		}else fragments.push_back( Fragment( null, null, id, mass ) );
	}
	//Transitions
	unsigned int numt; ifs.read(reinterpret_cast<char *>(&numt), sizeof(numt));
	for( int i = 0; i < numt; i++ ){
		int fromid; ifs.read(reinterpret_cast<char *>(&fromid), sizeof(fromid));
		int toid; ifs.read(reinterpret_cast<char *>(&toid), sizeof(toid));
		FeatureVector *fv = new FeatureVector();
		unsigned int num_set; ifs.read(reinterpret_cast<char *>(&num_set), sizeof(num_set));
		unsigned int fv_len; ifs.read(reinterpret_cast<char *>(&fv_len), sizeof(fv_len));
		for( unsigned int j = 0; j < num_set; j++ ){
			unsigned int fidx; ifs.read(reinterpret_cast<char *>(&fidx), sizeof(fidx));
			fv->addFeatureAtIdx( 1.0, fidx );
		}
		if( fv->getTotalLength() != fv_len )
			fv->addFeatureAtIdx( 0.0, fv_len-1 );
		transitions.push_back( Transition(fromid, toid, &null) );
		transitions[i].setTmpFV(fv);
	}

	//Create from_id and to_id maps
	to_id_tmap.resize(fragments.size());
	from_id_tmap.resize(fragments.size());
	std::vector<Transition>::const_iterator it = transitions.begin();
	for( int idx = 0; it != transitions.end(); ++it, idx++ ){
		to_id_tmap[it->getToId()].push_back(idx);
		from_id_tmap[it->getFromId()].push_back(idx);
	}

}

void EvidenceFragmentGraph::writeFragmentsOnly( std::ostream &out ) const{

	std::vector<EvidenceFragment>::const_iterator it = fragments.begin();
	for( ; it != fragments.end(); ++it ){
		out << it->getId() << " ";
		out << std::setprecision(10) << it->getMass() << " ";
		out << *it->getIonSmiles() << std::endl;
	}
}

//Write the FragmentGraph to file (formerly the transition output - without feature details)
void EvidenceFragmentGraph::writeFullGraph( std::ostream &out ) const{

	//Fragments
	out << fragments.size() << std::endl;
	writeFragmentsOnly( out );
	out << std::endl;

	//Transitions
	std::vector<Transition>::const_iterator itt = transitions.begin();
	for( ; itt != transitions.end(); ++itt ){
		out << itt->getFromId() << " ";
		out << itt->getToId() << " ";
		out << *itt->getNLSmiles() << std::endl;
	}
}

void FragmentGraph::deleteMolsForTransitionAtIdx( int index ){
	transitions[index].deleteNeutralLoss();
	transitions[index].deleteIon();
}

void FragmentGraph::clearAllSmiles(){
	std::vector<Fragment>::iterator it = fragments.begin();
	for( ; it != fragments.end(); ++it ) it->clearSmiles();
};

bool EvidenceFragmentGraph::fragmentIsRedundant( unsigned int fidx, std::vector<int> &annotated_flags, std::vector<int> &direct_flags ) const{

	if( annotated_flags[fidx] ) return false;
	std::vector<int>::const_iterator it = from_id_tmap[fidx].begin();
	for( ; it != from_id_tmap[fidx].end(); ++it ){
		int cidx = transitions[*it].getToId();
		if( !fragmentIsRedundant( cidx, annotated_flags, direct_flags ) && !direct_flags[cidx] ) return false;
	}
	return true;
}

void EvidenceFragmentGraph::setFlagsForDirectPaths(std::vector<int> &direct_flags, unsigned int fidx, std::vector<int> &annotated_flags ) const{
	
	if( !annotated_flags[fidx] ) return;
	direct_flags[fidx] = 1;
	std::vector<int>::const_iterator it = from_id_tmap[fidx].begin();
	for( ; it != from_id_tmap[fidx].end(); ++it ){
		int cidx = transitions[*it].getToId();
		setFlagsForDirectPaths( direct_flags, cidx, annotated_flags );  
	}
}




