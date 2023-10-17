/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# milp.c
#
# Description: 	Functions for running the milp solver looking for
#				valid assignment of bonds and hydrogens.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "MILP.h"
#include "lp_lib.h"
#include <GraphMol/RingInfo.h>
#include <GraphMol/ROMol.h>

//#ifndef __DEBUG_CONSTRAINTS__
//#define __DEBUG_CONSTRAINTS__

int MILP::runSolver( std::vector<int> &output_bmax, bool allow_lp_q, int max_free_pairs){
   //Uses lp_solve: based on demonstration code provided.
  unsigned int numunbroken;
  int broken, fragidx, origval;
  int Ncol, *colno = NULL, i, ret = 0, output_max_e=0;
  REAL *row = NULL;
  lprec *lp;

  // Variables are: num electron pairs
  // added per bond-> i.e. 1 = single, 2 = double, 3 = triple,
  // lone pair added per bond due to beginning atom,
  // lone pair added per bond due to end atom (at most one total),
  // place for charge due to H loss per atom (at most one total)
  int num_bonds = mol->getNumBonds();
  int num_atoms = mol->getNumAtoms();
  Ncol = num_bonds*3 + num_atoms;
  lp = make_lp(0, Ncol);
  if(lp == NULL)
    ret = 1; /* couldn't construct a new model... */

  if(ret == 0) {
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL))
      ret = 2;
  }

  //Bond Constraints
  int min_single_bonds = 0;
  if(ret == 0) {
    set_add_rowmode(lp, TRUE);

	//All bonds must be at most TRIPLE bonds ( <= 3 ) 
	//except broken bonds ( <= 0 ) and ring bonds ( <= 2 )
	for( i = 0; i < num_bonds && ret == 0; i++ ){
		set_int(lp, i+1, TRUE); //sets variable to integer

		RDKit::Bond *bond = mol->getBondWithIdx(i);
		int limit = 0;	  //bonds that are broken or in the other fragment are limited to 0
		int end_lp_limit = 0, begin_lp_limit = 0;	//bonds for which there is no lone pair to donate (or broken, or in other fragment) are limited to 0
		bond->getProp("Broken", broken);
		RDKit::Atom *begin_atom = bond->getBeginAtom();
		begin_atom->getProp("FragIdx", fragidx);
		int min_limit = 0;
		if( !broken && fragidx == fragmentidx ){
			limit = 3; min_limit = 1; min_single_bonds++;
			bond->getProp("NumUnbrokenRings", numunbroken);
			if( numunbroken > 0 ) limit = 2;
			else{ 
				begin_lp_limit = allow_lp_q && getAtomLPLimit(begin_atom);
				end_lp_limit = allow_lp_q && getAtomLPLimit(bond->getEndAtom());
			}
		}
		//Valence limit constraint
		colno[0] = i + 1; //variable idx i.e. bond
		row[0] = 1;		  //multiplier
		colno[1] = num_bonds + i + 1; //variable idx i.e. begin lone pair bond
		row[1] = 1;		  //multiplier
		colno[2] = 2*num_bonds + i + 1; //variable idx i.e. end lone pair bond
		row[2] = 1;		  //multiplier

		if(!add_constraintex(lp, 3, row, colno, LE, limit)){
		  ret = 3;
		  break;
		}
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( 3, colno, false, limit );
		#endif

		//Minimum constraint (LP bond + standard bond is at least a single bond)
		colno[0] = i + 1;
		row[0] = 1;		  //multiplier	
		if(!add_constraintex(lp, 1, row, colno, GE, min_limit)){
		  ret = 3;
		  break;
		}
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( 1, colno, true, min_limit );
		#endif
		//Lone pair constraint due to begin atom
		colno[0] = num_bonds + i + 1;
		row[0] = 1;		  //multiplier
		if(!add_constraintex(lp, 1, row, colno, LE, begin_lp_limit)){
		  ret = 3;
		  break;
		}
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( 1, colno, false, begin_lp_limit );
		#endif

		//Lone pair constraint due to end atom
		colno[0] = 2*num_bonds + i + 1;
		row[0] = 1;		  //multiplier
		if(!add_constraintex(lp, 1, row, colno, LE, end_lp_limit)){
		  ret = 3;
		  break;
		}
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( 1, colno, false, end_lp_limit );
		#endif
	}

	//Add constraints for neighbouring ring bonds (can't have two double in a row)
	RDKit::RingInfo *rinfo = mol->getRingInfo();
	RDKit::RingInfo::VECT_INT_VECT brings = rinfo->bondRings();
	RDKit::RingInfo::VECT_INT_VECT::iterator bit = brings.begin();
	for( int ringidx = 0; bit != brings.end() && ret == 0; ++bit, ringidx++ ){

		if( ringidx == broken_ringidx ) continue;
		row[0] = 1; row[1] = 1; row[2] = 1; row[3] = 1;
		
		//Create a vector of flags indicating bonds included in the ring
		std::vector<int> ring_bond_flags(num_bonds);
		for( int i = 0; i < num_bonds;i++ ) ring_bond_flags[i] = 0;
		RDKit::RingInfo::INT_VECT::iterator it;
		for( it = bit->begin(); it != bit->end(); ++it ) ring_bond_flags[*it] = 1;
		
		//Traverse around the ring, creating the constraints
		RDKit::Bond *bond = mol->getBondWithIdx(*(bit->begin())); //Starting Bond	
		RDKit::Bond *start_bond = bond, *prev_bond = bond;
		RDKit::Atom *atom = bond->getBeginAtom();
		int first_flag = 1;
		while( ret == 0 && (first_flag || prev_bond != start_bond) ){
			bond = getNextBondInRing(bond, atom, ring_bond_flags);
			colno[0] = prev_bond->getIdx() + 1;
			colno[1] = bond->getIdx() + 1;
			colno[2] = prev_bond->getIdx() + 1 + num_bonds;
			colno[3] = bond->getIdx() + 1 + num_bonds;
			if(!add_constraintex(lp, 4, row, colno, LE, 3)) ret = 3;
			atom = bond->getOtherAtom( atom );
			prev_bond = bond;
			first_flag = 0;
		}
	}

	//Add atom valence constraints to neighbouring bonds
	int numlp = 0;	std::vector<int> lp_indexes(num_atoms, -1);
	for( i = 0; i < num_atoms && ret == 0; i++ ){
		RDKit::Atom *atom = mol->getAtomWithIdx(i);
		atom->getProp("FragIdx", fragidx);
		atom->getProp("OrigValence", origval);
		int ionic_q; atom->getProp("IonicFragmentCharge", ionic_q);
		int has_lp; atom->getProp("HasLP", has_lp);
		unsigned int num_ur; atom->getProp("NumUnbrokenRings", num_ur);		

		//Charge due to H loss constraints
		colno[0] = 3*num_bonds + i + 1;
		row[0] = 1;
		int hloss_allowed = (!has_lp) && (fragidx == fragmentidx) && (ionic_q == 0) && (num_ur == 0);
		if(!add_constraintex(lp, 1, row, colno, LE, hloss_allowed))
		  ret = 3;
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( 1, colno, false, hloss_allowed);
		#endif
		if( ret != 0 ) break;

		if( fragidx != fragmentidx ) continue;

		//Base valence constraints
                RDKit::ROMol::OEDGE_ITER cbond_beg,cbond_end;
		boost::tie(cbond_beg, cbond_end) = mol->getAtomBonds( atom );
		//RDKit::ROMol::OEDGE_ITER it = ip.first;
		int j = 0; if( hloss_allowed ) j = 1;
                for ( ; cbond_beg != cbond_end; ++cbond_beg ){
		//for( ; it != ip.second; ++it ){
			RDKit::Bond *cbond = (*mol)[*cbond_beg];
			int cbond_broken; cbond->getProp("Broken", cbond_broken);
			if(cbond_broken) continue;

			colno[j] = cbond->getIdx() + 1; //variable idx i.e. bond
			row[j++] = 1;		  //multiplier
			//For atoms that don't contribute the lone pairs, 
			//ensure lone pair bonds are counted towards valence total
			if( cbond->getBeginAtomIdx() != i ){
				colno[j] = cbond->getIdx() + 1 + num_bonds;
				row[j++] = 1;			
			}
			if( cbond->getEndAtomIdx() != i ){
				colno[j] = cbond->getIdx() + 1 + 2*num_bonds;
				row[j++] = 1;			
			}			
		}

		if( j > 0 ){
			int val_limit = origval;
			if( atom->getDegree() > origval )  val_limit = atom->getDegree();
			if(!add_constraintex(lp, j, row, colno, LE, val_limit))
			  ret = 3;
			#ifdef __DEBUG_CONSTRAINTS__
				printConstraint( j, colno, false, val_limit );
			#endif
		}

	}

	//Total Lone Pair bond and H loss constraints - at most one of either used
	if( ret == 0 ){
		int j = 0;
		for( i = num_bonds; i < 3*num_bonds + num_atoms; i++ ){
			colno[j] = i + 1; //variable idx
			row[j++] = 1;     //multiplier
		}
		if(!add_constraintex(lp, j, row, colno, LE, 1))
			ret = 3;
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( j, colno, false, 1 );
		#endif

	}

	//Add the maximum constraint (no point finding a solution with more electrons than we have)
	if( ret == 0 ){
		for( int j = 0; j < Ncol; j++ ){
			colno[j] = j+1; 
			row[j] = 1;	
		}
		if(!add_constraintex(lp, Ncol, row, colno, LE, max_free_pairs+min_single_bonds))
			ret = 3;
		#ifdef __DEBUG_CONSTRAINTS__
			printConstraint( Ncol, colno, false, max_free_pairs+min_single_bonds );
		#endif
	}

  }

  if(ret == 0) {
	// set the objective function = sum bond electrons + sum lone pairs used for bonding
	set_add_rowmode(lp, FALSE);
	for( int j = 0; j < Ncol; j++ ){
		colno[j] = j+1; 
		row[j] = 1;	
	}
	if(!set_obj_fnex(lp, Ncol, row, colno))
		ret = 4;
  }

  //Run optimization
  if(ret == 0) {
    set_maxim(lp);
    set_verbose(lp, IMPORTANT);
    ret = solve(lp);
	if(ret == OPTIMAL) ret = 0;
    else ret = 5;
  }

  //Extract Results
  output_bmax.resize( Ncol );
  if(ret == 0) {
	get_variables(lp, row);
	for( int j = 0; j < Ncol; j++) output_bmax[j] = (int)row[j];
	int obj = get_objective(lp);
	output_max_e = obj - min_single_bonds;
  }
  else for( int j = 0; j < Ncol; j++) output_bmax[j] = 0;

  //Combine the lone pair results
  for( int j=0; j< num_bonds; j++ ) output_bmax[j + num_bonds] += output_bmax[j + 2*num_bonds];
  //Condense the atom H loss charge position results to a single index (or -1 if none)
  int hloss_idx[2] = {-1,-1};
  for( int j=0; j < num_atoms; j++ ){ 
	  if( output_bmax[j + 3*num_bonds] ) hloss_idx[fragmentidx] = j;
  }
  output_bmax.resize(2*num_bonds + 2);
  output_bmax[2*num_bonds] = hloss_idx[0];
  output_bmax[2*num_bonds + 1] =  hloss_idx[1];

  // Free allocated memory
  if(row != NULL) free(row);
  if(colno != NULL) free(colno);
  if(lp != NULL) delete_lp(lp);
  
  return output_max_e;
}

//Helper function - allows traversal of a ring one bond at a time
RDKit::Bond *MILP::getNextBondInRing( RDKit::Bond *bond, RDKit::Atom *atom, std::vector<int> &ring_bond_flags){

	RDKit::ROMol::OEDGE_ITER beg_abond, end_abond;
	boost::tie(beg_abond, end_abond)= atom->getOwningMol().getAtomBonds( atom );
	for( ; beg_abond != end_abond; ++beg_abond ){
		int idx = (*mol)[*beg_abond]->getIdx();
		if( ring_bond_flags[idx] && idx != bond->getIdx() )
			return atom->getOwningMol().getBondWithIdx(idx);
	}
	return NULL;	//Error@!
}

void MILP::printConstraint( int num_terms, int *colno, bool ge, int val  ){
	if( !verbose) return;
	for(int k=0; k <num_terms; k++ ){ std::cout << colno[k]; if(k !=num_terms-1) std::cout << "+"; }
	if(ge) std::cout << " >= " << val << std::endl;
	else  std::cout << " <= " << val << std::endl;

} 

int MILP::getAtomLPLimit(RDKit::Atom *atom){
	int has_lp; atom->getProp("HasLP", has_lp);
	int val; atom->getProp("OrigValence", val);
	unsigned int num_unbroken; atom->getProp("NumUnbrokenRings", num_unbroken);
	return( has_lp && (atom->getDegree() <= val) && (num_unbroken == 0) );
}
