/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FragmentGenerator.h
#
# Description: 	FragmentTree class for holding the results of a generated
#				fragment tree.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MILP_H__
#define __MILP_H__

#include <GraphMol/ROMol.h>
#include <vector>

class MILP{

public:
	MILP( RDKit::ROMol *a_mol, int a_fragmentidx, int a_broken_ringidx, bool a_verbose )
		: mol( a_mol ), fragmentidx(a_fragmentidx), broken_ringidx(a_broken_ringidx), verbose(a_verbose) {};

	MILP( RDKit::ROMol *a_mol, int a_fragmentidx, bool a_verbose )
		: mol( a_mol ), fragmentidx(a_fragmentidx), broken_ringidx(-1), verbose(a_verbose) {};

	int runSolver( std::vector<int> &output_bmax, bool allow_lp_q, int max_free_pairs );
	int status;

private:
	RDKit::ROMol *mol;
	int fragmentidx;
	int broken_ringidx;		//Store the idx of any broken rings (or -1 if there are none).
	bool verbose;

	//Helper functions:
	//Allows traversal of a ring one bond at a time
	RDKit::Bond *getNextBondInRing( RDKit::Bond *bond, RDKit::Atom *atom, std::vector<int> &ring_bond_flags);
	//Checks whether an atom should be allowed a lone pair bond (not including those already using theirs to create an extra single bond)
	static int getAtomLPLimit(RDKit::Atom *atom);

	void printConstraint( int num_terms, int *colno, bool ge, int val  );

};

#endif // __MILP_H__