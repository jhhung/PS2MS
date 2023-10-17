/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentTreeNode.cpp
#
# Description: 	Contains Break and FragmentTreeNode classes, for use in
#				the fragment generation process.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "FragmentTreeNode.h"
#include "MILP.h"
#include "Features.h"
#include "Config.h"

#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SanitException.h>
#include <GraphMol/PeriodicTable.h>

#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>


void FragmentTreeNode::generateChildrenOfBreak( Break &brk ){
	
	bool verbose = false;

	int f0_max_e = 0, f1_max_e = 0, f0_max_out = 0, f1_max_out = 0;
	std::vector<int> f0_output_bmax, f1_output_bmax;
	
	int isringbrk = brk.isRingBreak();
	int brk_ringidx = isringbrk*brk.getRingIdx() - (1-isringbrk);
	
	//C how many electrons were allocated to each side originally
	std::pair<int,int> orig_epairs = computeOrigFreeElectronsPerFrag();

	//Determine the upper limit for how many electron pairs can be allocated
	int total_free_epairs = ion_free_epairs + !brk.isIonicBreak() + isringbrk;
	int f0_max_limit = std::min( total_free_epairs, orig_epairs.first + MAX_E_MOVE );
	int f1_max_limit = std::min( total_free_epairs, orig_epairs.second + MAX_E_MOVE );

	//Compute the max electron assignment for F0
	MILP f0_solver(ion.get(), 0, brk_ringidx, verbose );
	f0_max_e = f0_solver.runSolver(f0_output_bmax, true, f0_max_limit);
		
	//Compute the max electron assignment for F1
	if( brk.getBondIdx() != -1 ){
		MILP f1_solver(ion.get(), 1, brk_ringidx, verbose );
		f1_max_e = f1_solver.runSolver(f1_output_bmax, true, f1_max_limit);

		//Combine the electron allocations of the two fragments
		unsigned int N = f0_output_bmax.size()-2;
		for( unsigned int i = 0; i < N; i++ )
			f0_output_bmax[i] += f1_output_bmax[i];
		f0_output_bmax[N+1] = f1_output_bmax[N+1];
	}

	//Remove single bonds
	int num_bonds = f0_output_bmax.size()/2 - 1;
	for( unsigned int i = 0; i < num_bonds; i++ ){
		if( f0_output_bmax[i] >= 1 ) f0_output_bmax[i] -= 1;
	}

	//Iterate through the possible solutions, adding them to the node as children
	int min_f0 = total_free_epairs - f1_max_e;
	int max_f0 = f0_max_e;
	for( int charge_frag = 0; charge_frag <= 1; charge_frag++ ){
		for( int e_f0 = min_f0; e_f0 <= max_f0; e_f0++ ){
			addChild( e_f0, total_free_epairs, f0_output_bmax, brk, charge_frag);
		}
	}
}

void FragmentTreeNode::addChild( int e_f0, int e_to_allocate, std::vector<int> &output_bmax, Break &brk, int charge_frag ){

	RDKit::RWMol rwmol = *ion;	//Copy the ion
	int is_radical, is_negative;
	ion->getProp( "isRadical", is_radical );
	ion->getProp( "isNegative", is_negative );
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();

	//Array to record the total connected bond orders for each atom
	std::vector<int> ctd_bond_orders(rwmol.getNumAtoms(),0);

	//Set the correct bond orders on the fragments
	int broken, fragidx;
	int numbonds = rwmol.getNumBonds();
	int remaining_e[2] = {e_f0, e_to_allocate - e_f0};
	int allocated_e[2] = {e_f0, e_to_allocate - e_f0};
	std::vector<std::pair<int,int> > broken_specs;
	for( int i = 0; i < numbonds; i++ ){
		RDKit::Bond *bond = rwmol.getBondWithIdx(i);
		bond->getProp("Broken", broken);
		if( broken ){ 
			//Save the bond to remove at the end, removing it now causes problems with the iteration/indexing
			broken_specs.push_back( std::pair<int,int>(bond->getBeginAtomIdx(), bond->getEndAtomIdx()) );
			continue;
		}
		bond->setIsAromatic(false);
		bond->getBeginAtom()->getProp("FragIdx", fragidx);
		int num_to_add = 0;
		if( remaining_e[fragidx] > 0){
			num_to_add = std::min( output_bmax[i], remaining_e[fragidx] );
			bond->setBondType( RDKit::Bond::BondType(1 + num_to_add ) );
			remaining_e[fragidx] -= num_to_add;
		}
		else bond->setBondType( RDKit::Bond::SINGLE );
		
		//Accumulate the connected bond orders (for use during hydrogen allocation)
		ctd_bond_orders[ bond->getBeginAtomIdx() ] += (1 + num_to_add);
		ctd_bond_orders[ bond->getEndAtomIdx() ] += (1 + num_to_add);	
	}

	//If there's still remaining_e, allocate as lp bond (only use lp bond if necessary)
	if( remaining_e[0] > 0 || remaining_e[1] > 0 ){
		for( int i = 0; i < numbonds; i++ ){
			RDKit::Bond *bond = rwmol.getBondWithIdx(i);
			bond->getBeginAtom()->getProp("FragIdx", fragidx);
			if( output_bmax[i + numbonds] > 0 && remaining_e[fragidx] > 0 ){
				bond->setBondType( RDKit::Bond::BondType(1 + (int)(bond->getBondTypeAsDouble()) ) );
				ctd_bond_orders[ bond->getBeginAtomIdx() ] += 1;
				ctd_bond_orders[ bond->getEndAtomIdx() ] += 1;
				remaining_e[fragidx] -= 1;
				if( remaining_e[0] == 0 && remaining_e[1] == 0 ) break;
			}
		}
	}

	//If this requires H loss charge on the neutral loss abort
	if( remaining_e[1-charge_frag] > 0 ) return;

	//Remove the broken bonds
	std::vector<std::pair<int,int> >::iterator itr = broken_specs.begin();
	for( ; itr != broken_specs.end(); ++itr )
		rwmol.removeBond(itr->first, itr->second);

	//Fill in the Hydrogens
	int orig_val;
	int split_charge_fragidx = -1, nl_ionic_q = 0;
	for( unsigned int i = 0; i < rwmol.getNumAtoms(); i++ ){
		RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
		atom->setNoImplicit(true);
		atom->setNumExplicitHs( 0 );
		atom->setNumRadicalElectrons( 0 );
		atom->setIsAromatic(false);

		atom->getProp( "OrigValence", orig_val );
		int fragidx; atom->getProp("FragIdx",fragidx);
		
		//Ionic fragments must retain their charge and not attract hydrogens (for simplicity)
		int ionic_frag_q; atom->getProp("IonicFragmentCharge", ionic_frag_q);
		atom->setFormalCharge( ionic_frag_q );
		
		//Track any pre-existing charge in the neutral loss
		if( ionic_frag_q != 0 && fragidx != charge_frag ) 
			nl_ionic_q += ionic_frag_q;

		int val = ctd_bond_orders[i];
		int numH = orig_val - val - abs(ionic_frag_q);
		//Check for charge via lone pair bond
		if( numH < 0 ){ 
			numH = 0;	//In case of charge via lone pair bond
			atom->setFormalCharge( 1 );	
			if( fragidx != charge_frag ) 
				split_charge_fragidx = fragidx;
		}
		//Check for charge via H loss
		if( remaining_e[fragidx] > 0 && output_bmax[2*numbonds + fragidx] == i ){
			numH -= 1;
			atom->setFormalCharge( 1 );	
			remaining_e[fragidx] -= 1;
		}
		atom->setNumExplicitHs( numH );
		val = atom->calcExplicitValence();
	}

	//If there should be a split charge on the neutral loss 
	//(that isn't already covered by a negative ionic fragment), 
	//find somewhere for the negative charge
	if( (split_charge_fragidx == (1-charge_frag) && nl_ionic_q == 0) || nl_ionic_q > 0 ){
		bool qloc_found = false;
		for( unsigned int i = 0; i < rwmol.getNumAtoms(); i++ ){
			RDKit::Atom *atom = rwmol.getAtomWithIdx(i);		
			int fragidx; atom->getProp("FragIdx",fragidx);
			if( fragidx == (1-charge_frag) && atom->getFormalCharge() == 0 && atom->getTotalNumHs() > 0 ){
				atom->setFormalCharge(-1);
				alterNumHs( atom, -1 );
				atom->calcExplicitValence();
				qloc_found = true;
				break;
			}
		}
		//If we couldn't find anywhere, give up
		if( !qloc_found) return;
	}

	//For Hydrogen loss only, add the hydogen neutral loss
	if( brk.isHydrogenOnlyBreak() ){
		RDKit::RWMol *h2mol = RDKit::SmilesToMol("[HH]");
		fh->addLabels( h2mol );
		RDKit::Atom *atom = h2mol->getAtomWithIdx(0);
		atom->setProp("Root", 1);
		atom->setProp("FragIdx", 1);
		atom->setProp("NumUnbrokenRings",(unsigned int)0);
		atom->setProp("IonicFragmentCharge", 0 );
		rwmol.addAtom(atom, true );
		delete h2mol;
	}

	RDKit::MolOps::sanitizeMol(rwmol);
	//std::cout << "DEBUG:" << RDKit::MolToSmiles( rwmol ) << std::endl;

	int min_radical_frag = -1; int max_radical_frag = -1;
	if( is_radical ){ min_radical_frag = 0; max_radical_frag = 1; }
	for( int radical_frag = min_radical_frag; radical_frag <= max_radical_frag; radical_frag++ ){

		//Return to a kekulized molecule here before we start messing with it again
		RDKit::MolOps::Kekulize( rwmol );

		//Check for an already charge molecule (e.g. due to N lone pair bonding), or a molecule with split charge (in which case, assign charge there)
		boost::tuple<int,int,int> pidx_nidx_ridx( -1, -1, -1 );
		std::pair<int,int> qidx_ridx(-1, -1);
		boost::tuple<bool,bool,bool> qfound_oktogo_switchq = findAlreadyChargedOrSplitCharge( pidx_nidx_ridx, rwmol, charge_frag, radical_frag, is_negative );
		if( boost::get<0>(qfound_oktogo_switchq) && !boost::get<1>(qfound_oktogo_switchq) ) continue;
		if( !boost::get<0>(qfound_oktogo_switchq) || boost::get<2>(qfound_oktogo_switchq) ){
			
			//Assign the charge (and radical) location
			if( boost::get<2>(qfound_oktogo_switchq) )
				qidx_ridx = findChargeLocation( rwmol, 1-charge_frag, radical_frag, is_negative );
			else
				qidx_ridx = findChargeLocation( rwmol, charge_frag, radical_frag, is_negative );
			if( qidx_ridx.first < 0 ) continue;				   //If we can't find a charge location, don't include this one.
			if( is_radical && qidx_ridx.second < 0 ) continue; //If we can't find a radical location (if needed), don't include this one.
			assignChargeAndRadical( rwmol, qidx_ridx.first, qidx_ridx.second, is_negative );
		}

		//Separate the ion and nl mols and add the child node
		recordOrigAtomIdxs( rwmol );
		std::vector<romol_ptr_t> mols = RDKit::MolOps::getMolFrags( rwmol, false );
		if( mols.size() == 2 ){
			int mol0_q = RDKit::MolOps::getFormalCharge( *mols[0].get() );
			int ion_idx = (mol0_q == 0);

			//Record some of the properties of the break in the neutral loss
			//e.g. ring break? aromaticity etc.
			romol_ptr_t nl = mols[1-ion_idx];
			labelBreakPropertiesInNL( nl, ion, brk );

			//Add the child node
			std::vector<int> child_e_loc; createChildIonElectronLocRecord( child_e_loc, mols[ion_idx] );
			children.push_back( FragmentTreeNode( mols[ion_idx], nl, allocated_e[charge_frag], depth+1, fh, child_e_loc ) );
		}
		else if( mols.size() > 2 ){

			//Collect the NL and ion parts (which may be multiple fragments each)
			romol_ptr_t acc_mols[2];  bool set_flags[2] = {false, false};
			for( int mol_idx = 0; mol_idx < mols.size(); mol_idx++ ){
				int q =  RDKit::MolOps::getFormalCharge( *mols[mol_idx].get() );
				int fragidx; mols[mol_idx]->getAtomWithIdx(0)->getProp("FragIdx",fragidx);
				if( set_flags[fragidx] ) acc_mols[fragidx].reset(RDKit::combineMols( *acc_mols[fragidx], *mols[mol_idx] ));
				else{ set_flags[fragidx] = true; acc_mols[fragidx] = mols[mol_idx]; }
			}
			romol_ptr_t nl = acc_mols[1-charge_frag];
			labelBreakPropertiesInNL( nl, ion, brk );
			std::vector<int> child_e_loc; createChildIonElectronLocRecord( child_e_loc, acc_mols[charge_frag] );
			children.push_back( FragmentTreeNode( acc_mols[charge_frag], nl, allocated_e[charge_frag], depth+1, fh, child_e_loc ) );
		}

		//Undo the charge (and radical) assignment
		if( !boost::get<0>(qfound_oktogo_switchq) || boost::get<2>(qfound_oktogo_switchq) )
			undoChargeAndRadical( rwmol, qidx_ridx.first, qidx_ridx.second, is_negative );
		else //(Split charge will still require undo)
			undoAlreadyChargedOrSplitCharge( rwmol, pidx_nidx_ridx );
		
	}

	//For Hydrogen loss only, remove the added hydrogens (should be just the last atom, but doesn't assume...)
	if( brk.getBondIdx() == -1 ){
		for( int i = rwmol.getNumAtoms()-1 ; i >= 0; i-- ){
			RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
			if( atom->getSymbol() == "H" ) rwmol.removeAtom( atom );
		}
	}
}
	
	
std::pair<int,int> FragmentTreeNode::findChargeLocation( RDKit::RWMol &rwmol, int charge_side, int radical_side, bool is_negative ){ 
	
	std::pair<int,int> qidx_ridx(-1, -1);

	//Find a non-ring Nitrogen, Sulfur, Oxygen, Carbon, or a ring N,S,O,C (preferenced in that order)
	qidx_ridx.first = findAtomChargeLocationNSOC( rwmol, charge_side, is_negative );
	if( qidx_ridx.first < 0 ) return qidx_ridx;

	//Radical locations
	if( radical_side < 0) return qidx_ridx;
	if( radical_side == charge_side ){
		
		//For non-C charge locations, we will remove one electron from that atom to create the charged radical
		if( rwmol.getAtomWithIdx( qidx_ridx.first )->getSymbol() != "C" ){
			qidx_ridx.second = qidx_ridx.first;
			return qidx_ridx;
		}

		//Otherwise we need to find a multiple bond to remove an electron from 
		//to leave a charge on one side and a radical on the other
		RDKit::ROMol::BondIterator bi; 
		for( bi = rwmol.beginBonds(); bi != rwmol.endBonds(); ++bi ){
			//(Note: We assume here that any broken bonds have been removed beforehand so we don't check for broken-ness).
			int fragidx; (*bi)->getBeginAtom()->getProp("FragIdx", fragidx);
			if( fragidx != charge_side ) continue;
			if( (*bi)->getBondTypeAsDouble() > 1.5 
				&& (*bi)->getBeginAtom()->getSymbol() == "C"	//Must be between C atoms
				&& (*bi)->getEndAtom()->getSymbol() == "C" ){
				qidx_ridx.first = (*bi)->getBeginAtomIdx();
				qidx_ridx.second = (*bi)->getEndAtomIdx();
				return qidx_ridx;
			}
		}

		//Otherwise pull the electron from a sigma bond with hydrogen
		qidx_ridx.second = qidx_ridx.first;
		return qidx_ridx;

	}
	else{
		//Allocate the radical somewhere on the neutral loss (anywhere with a H to remove)
		for( int i = 0; i < rwmol.getNumAtoms(); i++ ){
			RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
			int fragidx; atom->getProp("FragIdx", fragidx);
			if( fragidx != radical_side ) continue;
			if( atom->getIsAromatic() ) continue;	//Don't allow on aromatic atom
			if( atom->getTotalNumHs() > 0 ){
				qidx_ridx.second = i;
				return qidx_ridx;
			} 
		}
	}

	return qidx_ridx;

}

boost::tuple<bool,bool,bool> FragmentTreeNode::findAlreadyChargedOrSplitCharge( boost::tuple<int,int,int> &out_pidx_nidx_ridx, RDKit::RWMol &rwmol, int charge_side, int radical_side, bool is_negative ){ 
	
	boost::tuple<bool,bool,bool> qfound_oktogo_switchq( false, false, false );
	boost::get<0>(out_pidx_nidx_ridx) = -1;
	boost::get<1>(out_pidx_nidx_ridx) = -1;
	boost::get<2>(out_pidx_nidx_ridx) = -1;

	//Disallow molecules with 2+ total charge -> can otherwise occur if lone pair and ionic occur together.
	if( RDKit::MolOps::getFormalCharge( rwmol ) > 1 ){
		boost::get<0>(qfound_oktogo_switchq) = true;
		boost::get<1>(qfound_oktogo_switchq) = false;
		return qfound_oktogo_switchq;
	}

	//Locate any charges and radicals on the molecule
	int frag_total_qs[2] = {0,0};
	std::vector<int> neg_idxs;
	int rad_idx = -1;
	RDKit::ROMol::AtomIterator ai; 
	for( ai = rwmol.beginAtoms(); ai != rwmol.endAtoms(); ++ai ){
		int fragidx; (*ai)->getProp("FragIdx", fragidx);
		int q = (*ai)->getFormalCharge();
		frag_total_qs[ fragidx ] += q;
		if( q < 0 ) neg_idxs.push_back((*ai)->getIdx());
		if( (*ai)->getNumRadicalElectrons() == 1 ) rad_idx = (*ai)->getIdx();	//Note: only an assigned rad idx is added to out_pidx_nidx_ridx, this one is existing.
	}

	//No charges found (or existing charges balance)
	if( frag_total_qs[charge_side] == 0 && frag_total_qs[1-charge_side] == 0 ) return qfound_oktogo_switchq;
	boost::get<0>(qfound_oktogo_switchq) = true;

	//Positive ESI or EI: Charges balance correctly, just need to check radicals
	if( !is_negative && frag_total_qs[charge_side] == 1 && frag_total_qs[1-charge_side] == 0 ){
		if( radical_side < 0 && rad_idx == -1 ) boost::get<1>(qfound_oktogo_switchq) = true;  //No radical required or found
		else if( radical_side >= 0 && rad_idx >= 0 ){
			int rfragidx; rwmol.getAtomWithIdx(rad_idx)->getProp("FragIdx", rfragidx);
			if( rfragidx == radical_side ) boost::get<1>(qfound_oktogo_switchq) = true;
		}
		else if( radical_side >= 0 ){	//Still need to add a radical somewhere.
			for( int i = 0; i < rwmol.getNumAtoms(); i++ ){
				RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
				int fragidx; atom->getProp("FragIdx", fragidx);
				if( fragidx != radical_side ) continue;
				if( atom->getFormalCharge() != 0 ) continue;
				if( atom->getIsAromatic() ) continue;
				if( atom->getTotalNumHs() > 0 ){
					alterNumHs( atom, -1 );
					atom->setNumRadicalElectrons(1);
					atom->calcExplicitValence();
					boost::get<1>(qfound_oktogo_switchq) = true;
					boost::get<2>(out_pidx_nidx_ridx) = i;	//Assigned rad_idx
					break;
				} 
			}
		}
		return qfound_oktogo_switchq;
	}

	//Negative ESI and already charged
	if( is_negative && frag_total_qs[charge_side] == -1 && frag_total_qs[1-charge_side] == 0){
		boost::get<1>(qfound_oktogo_switchq) = true;
		return qfound_oktogo_switchq;
	}

	//Split Charge (overall neutral molecule)
	if( !is_negative && frag_total_qs[charge_side] == 1 && frag_total_qs[1-charge_side] == -1 ){
		
		//Find a non-ionic negative idx on the neutral loss to add a charge/radical to
		std::vector<int>::iterator it = neg_idxs.begin();
		for( ; it != neg_idxs.end(); ++it ){
			RDKit::Atom *neg_atom = rwmol.getAtomWithIdx(*it);
			int nfragidx; neg_atom->getProp("FragIdx", nfragidx);
			int nionic_q; neg_atom->getProp("IonicFragmentCharge", nionic_q);
			if( nfragidx == (1-charge_side) && nionic_q == 0 ){
				neg_atom->setFormalCharge(0);
				if( radical_side < 0 )
					alterNumHs( neg_atom, +1 );
				else
					neg_atom->setNumRadicalElectrons(1);
				neg_atom->calcExplicitValence();
				RDKit::MolOps::sanitizeMol(rwmol);
				RDKit::MolOps::Kekulize(rwmol);
				boost::get<1>(out_pidx_nidx_ridx) = *it;
				boost::get<1>(qfound_oktogo_switchq) = true;
			}
		}

		//If that didn't work - prepare to try adding the charge on the neutral loss side instead of the ion
		if(!boost::get<1>(qfound_oktogo_switchq)){
			boost::get<1>(qfound_oktogo_switchq) = true;
			boost::get<2>(qfound_oktogo_switchq) = true;
		}
	}

	return qfound_oktogo_switchq;
}

void FragmentTreeNode::undoAlreadyChargedOrSplitCharge( RDKit::RWMol &rwmol, boost::tuple<int,int, int> &pidx_nidx_ridx ){
	
	//Positive Charge only (check if radical was added - and if so, remove)
	if( boost::get<2>(pidx_nidx_ridx) >= 0 ){
		RDKit::Atom *rad_atom = rwmol.getAtomWithIdx( boost::get<2>(pidx_nidx_ridx) );
		rad_atom->setNumRadicalElectrons(0);
		alterNumHs( rad_atom, +1 );
	}

	//Undo Split Charge 
	if( boost::get<1>(pidx_nidx_ridx) >= 0 ){
		RDKit::Atom *neg_atom = rwmol.getAtomWithIdx( boost::get<1>(pidx_nidx_ridx) );
		neg_atom->setFormalCharge(-1);
		if( neg_atom->getNumRadicalElectrons() > 0 )
			neg_atom->setNumRadicalElectrons(0);
		else
			alterNumHs( neg_atom, -1 );
		neg_atom->calcExplicitValence();
		RDKit::MolOps::sanitizeMol(rwmol);
	}
}

int FragmentTreeNode::findAtomChargeLocationNSOC(RDKit::RWMol &rwmol, int charge_side, bool is_negative){

	int output = -1;

	//Search for possible atom locations
	int Nidx = -1, Pidx = -1, Sidx = -1, Oidx = -1, Hidx = -1, Cidx = -1;
	int NidxR = -1, PidxR = -1, SidxR = -1, OidxR = -1, CidxR = -1;
	for( unsigned int i = 0; i < rwmol.getNumAtoms(); i++ ){
		RDKit::Atom *atom = rwmol.getAtomWithIdx(i);
		
		int fragidx; atom->getProp("FragIdx", fragidx);
		if( fragidx != charge_side ) continue;
		if( is_negative && atom->getTotalNumHs() == 0) continue;
		if( atom->getFormalCharge() != 0 ) continue;	//If it's already charged (due to split charge), don't use it again.

		unsigned int numrings; atom->getProp("NumUnbrokenRings", numrings);		
		if( numrings > 0 ){
			if( atom->getIsAromatic() ) continue;
			if( atom->getSymbol() == "N" && NidxR < 0 ) NidxR = i;
			else if( atom->getSymbol() == "P" && PidxR < 0 ) PidxR = i;
			else if( atom->getSymbol() == "S" && SidxR < 0 ) SidxR = i;
			else if( atom->getSymbol() == "O" && OidxR < 0 ) OidxR = i;
			else if( atom->getSymbol() == "C" && CidxR < 0 ) CidxR = i;
		}
		else{
			if( atom->getSymbol() == "N" && Nidx < 0 ) Nidx = i;
			else if( atom->getSymbol() == "P" && Pidx < 0 ) Pidx = i;
			else if( atom->getSymbol() == "S" && Sidx < 0 ) Sidx = i;
			else if( atom->getSymbol() == "O" && Oidx < 0 ) Oidx = i;
			else if( Hidx < 0 && (atom->getSymbol() == "Cl" || atom->getSymbol() == "I" || 
				     atom->getSymbol() == "Br" || atom->getSymbol() == "F") ) Hidx = i;
			else if( atom->getSymbol() == "C" && Cidx < 0 ) Cidx = i;
		}
	}

	//Decide on the location
	if( Nidx >= 0 ) output = Nidx;		 //Non-Ring Nitrogen
	else if( Pidx >= 0 ) output = Pidx;  //Non-Ring Phosphorus
	else if( Sidx >= 0 ) output = Sidx;  //Non-Ring Sulfur
	else if( Oidx >= 0 ) output = Oidx;  //Non-Ring Oxygen
	else if( Hidx >= 0 ) output = Hidx;  //Non-Ring Halogen
	else if( Cidx >= 0 ) output = Cidx;  //Non-Ring Carbon
	else if( NidxR >= 0 ) output = NidxR; //Non-aromatic ring Nitrogen
	else if( PidxR >= 0 ) output = PidxR; //Non-aromatic ring Phosphorus
	else if( SidxR >= 0 ) output = SidxR; //Non-aromatic ring Sulfur
	else if( OidxR >= 0 ) output = OidxR; //Non-aromatic ring Oxygen
	else if( CidxR >= 0 ) output = CidxR; //Non-aromatic ring Carbon

	return output;

}

void FragmentTreeNode::assignChargeAndRadical( RDKit::RWMol &rwmol, int charge_idx, int radical_idx, bool is_negative ){

	RDKit::Atom *atom = rwmol.getAtomWithIdx(charge_idx);
	if( is_negative ){
		atom->setFormalCharge(-1);
		alterNumHs( atom, -1 ); //Even-[H+] -> Even-[H+] + NL
	}else{
		atom->setFormalCharge(1);
		if( radical_idx < 0 ) alterNumHs( atom, +1 ); //Even+[H+] -> Even+[H+] + NL
		else{ 
			RDKit::Atom *rad_atom = rwmol.getAtomWithIdx( radical_idx );
			rad_atom->setNumRadicalElectrons(1); 
			if( radical_idx != charge_idx ){	//Else Odd+. -> Odd+. (lone pair or C.H) + NL
				int rad_fragidx; rad_atom->getProp("FragIdx", rad_fragidx);
				int charge_fragidx; atom->getProp("FragIdx", charge_fragidx);
				if( rad_fragidx != charge_fragidx ){ //Odd+. -> Even+[H+] + NL-[H.]
					alterNumHs( rad_atom, -1);
					alterNumHs( atom, +1);
				}
				else{	//Odd+. -> Odd+. (multiple bond) + NL
					RDKit::Bond *mult_bond = rwmol.getBondBetweenAtoms( radical_idx, charge_idx );
					mult_bond->setBondType( RDKit::Bond::BondType((int)(mult_bond->getBondTypeAsDouble() - 1.0) ) );
					alterNumHs( atom, 0);	//Leep num Hs the same but make then explicit
				}
			}
			rad_atom->calcExplicitValence();
		}
	}
	atom->calcExplicitValence();
	RDKit::MolOps::sanitizeMol(rwmol);
	RDKit::MolOps::Kekulize(rwmol);

}

void FragmentTreeNode::undoChargeAndRadical( RDKit::RWMol &rwmol, int charge_idx, int radical_idx, bool is_negative  ){

	RDKit::Atom *atom = rwmol.getAtomWithIdx(charge_idx);
	atom->setFormalCharge(0);
	if( is_negative ) alterNumHs( atom, +1); //Even-[H+] -> Even-[H+] + NL
	else{
		if( radical_idx < 0 ) alterNumHs( atom, -1 ); //Even+[H+] -> Even+[H+] + NL
		else{ 
			RDKit::Atom *rad_atom = rwmol.getAtomWithIdx( radical_idx );
			rad_atom->setNumRadicalElectrons(0); 
			if( radical_idx != charge_idx ){	//Else Odd+. -> Odd+. (lone pair) + NL
				int rad_fragidx; rad_atom->getProp("FragIdx", rad_fragidx);
				int charge_fragidx; atom->getProp("FragIdx", charge_fragidx);
				if( rad_fragidx != charge_fragidx ){ //Odd+. -> Even+[H+] + NL-[H.]
					alterNumHs( rad_atom, +1 );
					alterNumHs( atom, -1 );
				}
				else{	//Odd+. -> Odd+. (multiple bond) + NL
					RDKit::Bond *mult_bond = rwmol.getBondBetweenAtoms( radical_idx, charge_idx );
					mult_bond->setBondType( RDKit::Bond::BondType( (int)(mult_bond->getBondTypeAsDouble() + 1.0) ) );
				}
			}
			rad_atom->calcExplicitValence();
		}
	}
	atom->calcExplicitValence();
}

void FragmentTreeNode::labelBreakPropertiesInNL( romol_ptr_t &nl, romol_ptr_t &parent_ion, Break &brk ){

	//Ring Properties
	nl.get()->setProp("IsRingBreak",brk.isRingBreak());
	int is_arom = 0, is_dbl_arom = 0;
	if( brk.getBondIdx() != -1 ){
		RDKit::Bond *bond = parent_ion->getBondWithIdx( brk.getBondIdx() );
		bond->getProp("InAromaticRing", is_arom);
		bond->getProp("InDblAromaticRing", is_dbl_arom);
	}
	nl.get()->setProp("IsAromaticRingBreak", is_arom);
	nl.get()->setProp("IsAromaticDblRingBreak", is_dbl_arom);

	//Broken Bond Type Properties
	if( fh->getExecFlag(3) ){
		if( brk.isIonicBreak() )  nl.get()->setProp("BrokenOrigBondType", 6);
		else if( brk.isHydrogenOnlyBreak() ) nl.get()->setProp("BrokenOrigBondType", 7);
		else{
			RDKit::Bond *brokenbond = parent_ion.get()->getBondWithIdx( brk.getBondIdx() );
			int bondtype; brokenbond->getProp( "OrigBondType", bondtype );
			nl.get()->setProp("BrokenOrigBondType", bondtype);
		}
	}

}

void FragmentTreeNode::generateBreaks(std::vector<Break> &breaks, bool include_H_only_loss){

	int israd = moleculeHasSingleRadical( ion.get() );
	int num_ionic = countNumIonicFragments( ion.get() );
	RDKit::PeriodicTable *pt = RDKit::PeriodicTable::getTable();

	//Populate the Ring Info for the ion and set the NumUnbrokenRings and OrigValence properties
	RDKit::MolOps::findSSSR( *ion.get() );	
	RDKit::RingInfo *rinfo = ion.get()->getRingInfo();
	RDKit::ROMol::AtomIterator ai; 
	for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai ){
		(*ai)->setProp("NumUnbrokenRings", rinfo->numAtomRings((*ai)->getIdx()) );
		
		//Fetch or compute the valence of the atom in the input molecule (we disallow valence changes for now)
		int origval = -1;
		unsigned int num_val = pt->getValenceList( (*ai)->getSymbol() ).size();
		int def_val = pt->getDefaultValence( (*ai)->getSymbol() );
		if( num_val == 1 && def_val != -1 ) origval = def_val; //Hack to cover many cases - which can otherwise get complicated
		else{
			//This seems to work in most cases....
			origval = (*ai)->getExplicitValence() + (*ai)->getImplicitValence() + (*ai)->getNumRadicalElectrons();
			if( 4 - pt->getNouterElecs((*ai)->getAtomicNum()) > 0) origval += (*ai)->getFormalCharge();
			else origval -= (*ai)->getFormalCharge();
		}

		(*ai)->setProp("OrigValence", origval);
		(*ai)->setProp("Root", 0);
		(*ai)->setProp("OtherRoot", 0);

		//Create breaks for any implied ionic bonds (-1 bond_idx, and ring_idx overloaded with the atom idx)
		int ionic_frag_q; (*ai)->getProp( "IonicFragmentCharge", ionic_frag_q );
		if( ionic_frag_q != 0 && num_ionic != ion.get()->getNumAtoms() ) //(must have at least one non-ionic atom to break further)
			breaks.push_back( Break((*ai)->getIdx(), true, computeNumIonicAlloc(num_ionic-1)) );
	}

	//Generate Non-Ring Breaks
	for( unsigned int bidx = 0; bidx < ion.get()->getNumBonds(); bidx++ ){
		RDKit::Bond *bond = ion.get()->getBondWithIdx(bidx);
		bond->setProp("Broken",0);
		bond->setProp("NumUnbrokenRings", rinfo->numBondRings(bidx));
		if( rinfo->numBondRings(bidx) == 0 )
			breaks.push_back( Break(bidx, false, computeNumIonicAlloc(num_ionic)) );
	}
	
	//Ring Breaks
	RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
	RDKit::RingInfo::VECT_INT_VECT::iterator bit = brings.begin();
	for( int ringidx = 0; bit != brings.end(); ++bit, ringidx++ ){

		//Only include rings with size less than MAX_BREAKABLE_RING_SIZE 
		//(larger rings can exist, but will not break since it blows out the computation)
		if( bit->size() > MAX_BREAKABLE_RING_SIZE ) continue;

		//All pairs of bonds within the ring, that don't
		//belong to any other ring
		RDKit::RingInfo::INT_VECT::iterator it1, it2;
		for( it1 = bit->begin(); it1 != bit->end(); ++it1 ){
			if( rinfo->numBondRings(*it1) != 1 ) continue;
			for( it2 = it1 + 1; it2 != bit->end(); ++it2 ){
				if( rinfo->numBondRings(*it2) != 1 ) continue;

				breaks.push_back( Break( std::pair<int, int>(*it1, *it2), ringidx, computeNumIonicAlloc(num_ionic)) );
			}
		}
	}

	//Hydrogen only breaks (-1 bond_idx, and -1 ring_idx)
	if( include_H_only_loss ) breaks.push_back( Break()  );
}

void FragmentTreeNode::applyBreak(Break &brk, int ionic_allocation_idx){

	//Initialise fragment idxs
	RDKit::ROMol::AtomIterator ai; 
	for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai )
		(*ai)->setProp( "FragIdx", 0 );

	//Label the broken bond, root atoms and fragment indexes
	if( brk.isHydrogenOnlyBreak() ){
		//Hydrogen Only Break - No explicit bond to break, set arbitrary root atom (pick a C if there is one)
		int root_idx = 0;
		for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai )
			if( (*ai)->getSymbol() == "C" ){ root_idx = (*ai)->getIdx(); break; }
		ion.get()->getAtomWithIdx(root_idx)->setProp("Root",1);	
	}
	else if( brk.isIonicBreak() ){
		//Implied Ionic Break - No explicit bond to break, set root atoms
		RDKit::Atom *ionic_atom = ion.get()->getAtomWithIdx( brk.getIonicIdx() );
		ionic_atom->setProp("Root",1);
		ionic_atom->setProp("FragIdx", 1);
		int ionic_q; ionic_atom->getProp("IonicFragmentCharge", ionic_q);
		int root_idx = -1;	//Pick a root atom on the other fragment
		for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai ){
			int a_ionic_q; (*ai)->getProp("IonicFragmentCharge", a_ionic_q);
			if( a_ionic_q != 0.0 ) continue;	//Don't pick another ionic fragment
			if( (*ai)->getFormalCharge() == -ionic_q ){ root_idx = (*ai)->getIdx(); break; }
			if( root_idx < 0 && (*ai)->getIdx() != brk.getIonicIdx() ) root_idx = (*ai)->getIdx();
		}
		ion.get()->getAtomWithIdx(root_idx)->setProp("Root",1);
	}
	else{ //Standard Bond Break or Ring Break
		RDKit::Bond *broken_bond = ion.get()->getBondWithIdx( brk.getBondIdx() );
		broken_bond->setProp( "Broken", 1 );
		broken_bond->getBeginAtom()->setProp("Root",1);
		broken_bond->getEndAtom()->setProp("Root",1);
		if( brk.isRingBreak() ){
			RDKit::Bond *broken_bond2 = ion.get()->getBondWithIdx( brk.getSecondBondIdx() );
			broken_bond2->setProp( "Broken", 1 );
			broken_bond2->getBeginAtom()->setProp("OtherRoot",1);
			broken_bond2->getEndAtom()->setProp("OtherRoot",1);
		}
		allocatedCtdToFragment( ion.get(), broken_bond->getBeginAtom() );
	}

	//Assign fragment indexes for any (other) ionic fragments (which are otherwise always FragIdx=0)
	int ionic_idx = ionic_allocation_idx;
	for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms() && ionic_idx > 0; ++ai ){
		int ionic_q; (*ai)->getProp( "IonicFragmentCharge", ionic_q );
		if( ionic_q != 0 && (!brk.isIonicBreak() || (*ai)->getIdx() != brk.getIonicIdx()) ){
			if( ionic_idx & 0x1 ) (*ai)->setProp("FragIdx", 1);
			ionic_idx = ionic_idx >> 1;
		}
	}

	//Decrement the NumUnbrokenRings setting for atoms and bonds in a broken ring
	if( brk.isRingBreak() ){
		int ringidx = brk.getRingIdx();
		RDKit::RingInfo *rinfo = ion.get()->getRingInfo();
		RDKit::RingInfo::INT_VECT::iterator it;		

		RDKit::RingInfo::VECT_INT_VECT arings = rinfo->atomRings();
		for( it = arings[ringidx].begin(); it != arings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Atom *atom = ion.get()->getAtomWithIdx(*it);
			atom->getProp("NumUnbrokenRings", numrings);
			atom->setProp("NumUnbrokenRings", numrings-1);
		}

		RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
		for( it = brings[ringidx].begin(); it != brings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Bond *bond = ion.get()->getBondWithIdx(*it);
			bond->getProp("NumUnbrokenRings", numrings);
			bond->setProp("NumUnbrokenRings", numrings-1);
		}
	}
}

void FragmentTreeNode::undoBreak(Break &brk, int ionic_allocation_idx){

	//Label the broken bond and root atoms
	if( brk.isHydrogenOnlyBreak() ){
		int root_idx = 0;
		RDKit::ROMol::AtomIterator ai;
		for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai )
			if( (*ai)->getSymbol() == "C" ){ root_idx = (*ai)->getIdx(); break; }
		ion.get()->getAtomWithIdx(root_idx)->setProp("Root",0);
	}
	else if( brk.isIonicBreak() ){
		
		//Label the ionic atom
		RDKit::Atom *ionic_atom = ion.get()->getAtomWithIdx( brk.getIonicIdx() );
		ionic_atom->setProp("Root", 0);
		ionic_atom->setProp("FragIdx", 0);
		int ionic_q = ionic_atom->getFormalCharge();
		
		//Pick a root atom on the other fragment (ideally one of opposite charge to the ionic atom, else any)
		int root_idx = -1;	RDKit::ROMol::AtomIterator ai;
		for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms(); ++ai ){
			if( (*ai)->getFormalCharge() == -ionic_q ){ root_idx = (*ai)->getIdx(); break; }
			if( root_idx < 0 && (*ai)->getIdx() != brk.getIonicIdx() ) root_idx = (*ai)->getIdx();
		}
		ion.get()->getAtomWithIdx(root_idx)->setProp("Root",0);
	}
	else{
		RDKit::Bond *broken_bond = ion.get()->getBondWithIdx( brk.getBondIdx() );
		broken_bond->setProp( "Broken", 0 );
		broken_bond->getBeginAtom()->setProp("Root",0);
		broken_bond->getEndAtom()->setProp("Root",0);	
	}

	if( brk.isRingBreak() ){
		RDKit::Bond *broken_bond2 = ion.get()->getBondWithIdx( brk.getSecondBondIdx() );
		broken_bond2->setProp( "Broken", 0 );
		broken_bond2->getBeginAtom()->setProp("OtherRoot",0);
		broken_bond2->getEndAtom()->setProp("OtherRoot",0);
	}

	//Un-assign fragment indexes for any ionic fragments (which are otherwise always FragIdx=0)
	int ionic_idx = ionic_allocation_idx; RDKit::ROMol::AtomIterator ai; 
	for( ai = ion.get()->beginAtoms(); ai != ion.get()->endAtoms() && ionic_idx > 0; ++ai ){
		int ionic_q; (*ai)->getProp( "IonicFragmentCharge", ionic_q );
		if( ionic_q != 0 && (!brk.isIonicBreak() || (*ai)->getIdx() != brk.getIonicIdx()) ){
			if( ionic_idx & 0x1 ) (*ai)->setProp("FragIdx", 0);
			ionic_idx = ionic_idx >> 1;
		}
	}

	//Increment the NumUnbrokenRings setting for atoms and bonds in a broken ring
	if( brk.isRingBreak() ){
		int ringidx = brk.getRingIdx();
		RDKit::RingInfo *rinfo = ion.get()->getRingInfo();
		RDKit::RingInfo::INT_VECT::iterator it;		

		RDKit::RingInfo::VECT_INT_VECT arings = rinfo->atomRings();
		for( it = arings[ringidx].begin(); it != arings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Atom *atom = ion.get()->getAtomWithIdx(*it);
			atom->getProp("NumUnbrokenRings", numrings);
			atom->setProp("NumUnbrokenRings", numrings+1);
		}

		RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
		for( it = brings[ringidx].begin(); it != brings[ringidx].end(); ++it ){
			unsigned int numrings;
			RDKit::Bond *bond = ion.get()->getBondWithIdx(*it);
			bond->getProp("NumUnbrokenRings", numrings);
			bond->setProp("NumUnbrokenRings", numrings+1);
		}
	}
}

void FragmentTreeNode::allocatedCtdToFragment( RDKit::ROMol *romol, RDKit::Atom *atom ){
	
	int broken, fragidx;
	atom->setProp( "FragIdx", 1 );
	RDKit::ROMol::ADJ_ITER_PAIR itp = romol->getAtomNeighbors( atom );
	for( ; itp.first != itp.second; ++itp.first ){
		RDKit::Atom *nbr_atom = romol->getAtomWithIdx(*itp.first);
		nbr_atom->getProp("FragIdx",fragidx);
		RDKit::Bond *bond = romol->getBondBetweenAtoms( atom->getIdx(), nbr_atom->getIdx() );
		bond->getProp("Broken", broken);
		if(!fragidx && !broken) allocatedCtdToFragment( romol, nbr_atom );
	}
}

//Utility function to label properties of the ion: is it negative? is it a radical?
void FragmentTreeNode::labelIonProperties(){

	int is_radical = moleculeHasSingleRadical( ion.get() );
	int is_negative = ( RDKit::MolOps::getFormalCharge( *ion ) < 0 );
	ion->setProp( "isRadical", is_radical );
	ion->setProp( "isNegative", is_negative );
}

//Utility function to count the number of ionic fragments in a molecule
int FragmentTreeNode::countNumIonicFragments( const RDKit::ROMol *romol ){
	int num_ionic = 0;
	for(RDKit::ROMol::ConstAtomIterator ait =romol->beginAtoms(); ait!=romol->endAtoms(); ++ait){
        int ionic_frag_q; (*ait)->getProp( "IonicFragmentCharge", ionic_frag_q );
		if( ionic_frag_q != 0 )num_ionic++;
    }
	return num_ionic;
}

int FragmentTreeNode::computeNumIonicAlloc( int num_ionic_fragments ){
	return 1 << num_ionic_fragments;
}

std::pair<int,int> FragmentTreeNode::computeOrigFreeElectronsPerFrag(){
	
	std::pair<int,int> output(0,0);
	RDKit::ROMol::AtomIterator ait =ion.get()->beginAtoms();
	for( ; ait!=ion.get()->endAtoms(); ++ait){
		int fragidx; (*ait)->getProp("FragIdx", fragidx);
		if( fragidx == 0 ) output.first += e_loc[(*ait)->getIdx()]; 
		else output.second += e_loc[(*ait)->getIdx()]; ;
	}
	output.first /= 2;
	output.second /= 2;
	return output;
}

void FragmentTreeNode::createChildIonElectronLocRecord( std::vector<int> &child_e_loc, romol_ptr_t childmol ){

	child_e_loc.resize(childmol.get()->getNumAtoms());
	RDKit::ROMol::AtomIterator ait = childmol.get()->beginAtoms();
	for( ; ait!=childmol.get()->endAtoms(); ++ait){	
		unsigned int idx_in_parent; (*ait)->getProp("OrigIdx", idx_in_parent);
		child_e_loc[(*ait)->getIdx()] = e_loc[idx_in_parent];
	}
}

void FragmentTreeNode::recordOrigAtomIdxs( RDKit::RWMol &rwmol ){
	RDKit::ROMol::AtomIterator ait = rwmol.beginAtoms();
	for( ; ait!= rwmol.endAtoms(); ++ait){
		unsigned int idx = (*ait)->getIdx();
		(*ait)->setProp("OrigIdx", idx);
	}
}