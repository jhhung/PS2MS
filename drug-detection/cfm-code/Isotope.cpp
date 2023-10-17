/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Isotope.cpp
#
# Description:  Class for generation of isotope peaks for a molecular formula.
#				Provides wrapper to emass from:
#		
# A. Rockwood and P. Haimi, "Efficient calculation of accurate masses of isotopic peaks.", 
# Journal of the American Society for Mass Spectrometry, 17:3 p415-9 2006.				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Isotope.h"
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/MolOps.h>
#include <fstream>

typedef unsigned long ulong;

static const double DUMMY_MASS = -10000000;

void IsotopeCalculator::computeIsotopeSpectrum( Spectrum &output, const romol_ptr_t mol, long charge ){
	
	FormMap fm;
	setFormulaMap( fm, mol );
	
	//Adjust H count for charge
	int mol_q = RDKit::MolOps::getFormalCharge( *mol.get() );
	if( mol_q != charge ){
		if( mol_q == 0 && charge > 0 )
			fm[em["H"]] += 1;
		if( mol_q == 0 && charge < 0 )
			fm[em["H"]] -= 1;
	}

	Pattern result(1);
	Pattern tmp;

    // initialize the result
    result.resize(1);
    result.front().mass = 0.0;
    result.front().rel_area = 1.0;

	calculate(tmp, result, fm, 0, charge);

    if(verbose) print_pattern(result, 10);
	print_to_output(output, result);
	output.normalizeAndSort();

}

//This function sets the formula map (atom_idx -> count) structure for the input molecule
void IsotopeCalculator::setFormulaMap( FormMap &output, const romol_ptr_t mol ){

	if( em.find("H") == em.end() ){
		std::cout << "Could not find symbol H in ElemMap" << std::endl;
		throw IsotopeCalculationException();	
	}
	std::string hstr = "H";
	size_t hidx = em[hstr];
	output[hidx] = 0;

	RDKit::ROMol::AtomIterator ait = mol.get()->beginAtoms();
	for( ; ait!=mol.get()->endAtoms(); ++ait){
		std::string symbol = (*ait)->getSymbol();
		if( em.find(symbol) == em.end() ){
			std::cout << "Could not find symbol " << symbol << " in ElemMap" << std::endl;
			throw IsotopeCalculationException();
		}
		size_t atomidx = em[symbol];
		if( output.find(atomidx) == output.end() ) output[atomidx] = 1;
		else output[atomidx] += 1;

		output[hidx] += (*ait)->getTotalNumHs();
	}
}

//The remainder of the functions are copied directly from emass (with minor mods)
//-------------------------------------------------------------------------------
/*This collective work is Copyright (C)2005 by Perttu Haimi 
Individual portions may be copyright by individual 
contributors, and are included in this collective work with 
permission of the copyright owners.

All rights reserved.

Redistribution and use in source and binary forms, 
with or without modification, are permitted provided 
that the following conditions are met:

    * Redistributions of source code must retain the 
      above copyright notice, this list of conditions 
      and the following disclaimer.
    * Redistributions in binary form must reproduce 
      the above copyright notice, this list of conditions 
      and the following disclaimer in the documentation 
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors 
      may be used to endorse or promote products derived 
      from this software without specific prior written 
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//-------------------------------------------------------------------------------

void IsotopeCalculator::init_data( const char *filename ){
  
  std::ifstream f(filename);
  
  if(f.fail()){
      std::cout << "Error opening ISOTOPE.DAT file. Please ensure it is in the executable directory." << std::endl;
	  throw IsotopeCalculationException();
  }

  sad.clear();
  em.clear();

  std::string line;
  ulong elemindex = 0;
  int state = 0;
  while(getline(f, line)) {
    std::istringstream ist(line);
    std::string element;
    switch(state) {
    case 0: // new element
      ist >> element;
      em[element] = elemindex;
      sad.push_back(SuperAtomList(1));
      sad.back().reserve(8); // reserve room for 8 superatoms
      elemindex++;
      state = 1;
      break;
    case 1: // isotope
      ipeak p;
      Pattern & idist = sad.back()[0];
      if(ist >> p.mass >> p.rel_area) {
	// fill the gaps in the patterns with zero abundancy peaks
	if(idist.size() > 0) {
	  double prevmass = idist.back().mass;
	  for(int i = 0; i < int(p.mass - prevmass - 0.5); i++) {
	    ipeak filler;
	    filler.mass = DUMMY_MASS;
	    filler.rel_area = 0;
	    idist.push_back(filler);
	  }
	}
	// insert the peak
	idist.push_back(p);                                                        
      } else  
	state = 0; // no more isotope data
      break;
    }
  }  
  f.close();
}


// Merge two patterns to one.
void IsotopeCalculator::convolute_basic(Pattern & h, const Pattern & g, const Pattern & f)
{
  h.clear();
  size_t g_n = g.size();
  size_t f_n = f.size();
  if(g_n == 0 || f_n == 0)
     return;
  for(size_t k = 0; k < g_n + f_n - 1; k++) {
    double sumweight = 0, summass = 0;
    size_t start = k < (f_n - 1) ? 0 : k - f_n + 1; // max(0, k-f_n+1)
    size_t end = k < (g_n - 1) ? k : g_n - 1;       // min(g_n - 1, k)
    for(size_t i = start; i <= end; i++) {
      double weight = g[i].rel_area * f[k - i].rel_area;
      double mass = g[i].mass + f[k - i].mass;
      sumweight += weight;
      summass += weight * mass;
    }
    ipeak p;
    if(sumweight == 0)
      p.mass = DUMMY_MASS;
    else
      p.mass = summass / sumweight;
    p.rel_area = sumweight;
    h.push_back(p);
  }
}

// Prune the small peaks from both sides but
// leave them within the pattern.
void IsotopeCalculator::prune(Pattern & f, double limit)
{
  // prune the front
  Pattern::iterator i = f.begin();
  while(i != f.end()) {
    if((*i).rel_area > limit)
      break;
    i++;
  }
  f.erase(f.begin(), i);

  // prune the end
  while(1) {
    if(f.size() == 0)
      break;
    if(f.back().rel_area > limit)
      break;
    f.pop_back();
  } 
}

void IsotopeCalculator::calculate(Pattern & tmp, Pattern & result, FormMap & fm, double limit, long charge)
{
  for(FormMap::iterator i = fm.begin(); i != fm.end(); i++) {
    size_t atom_index = (*i).first;
    SuperAtomList sal = sad[atom_index];
    ulong n = (*i).second;
    ulong j = 0;
    while(n > 0) {
      size_t sz = sal.size();
      if(j == sz) { 
	sal.resize(sz + 1); // Make new superatom from previous
                            // largest superatom. We are trying to
                            // avoid copying on assignment here.
	convolute_basic(sal[j], sal[j - 1], sal[j - 1]);
	prune(sal[j], limit);
      }
      if(n & 1) { // digit is 1, convolute result
	convolute_basic(tmp , result, sal[j]);
	prune(tmp, limit);
	swap(tmp, result);  // Hopefully the swap implementation
                            // will not copy all elements.
      }
      n >>= 1; 
      j++;
    }
  }
  
  // take charge into account
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(charge > 0)
      (*i).mass = (*i).mass / abs(charge) - MASS_ELECTRON;
    else if (charge < 0)
      (*i).mass = (*i).mass / abs(charge) + MASS_ELECTRON;
  }
}


void IsotopeCalculator::print_pattern(Pattern & result, int digits)
{
  // find the maximum
  double max_area = 0;
  double sum_area = 0;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(max_area < (*i).rel_area)
      max_area = (*i).rel_area;
    sum_area += (*i).rel_area;
  }
  if(max_area == 0)
    return; // empty pattern

  std::wcout.setf(std::ios::fixed);
  std::wcout.precision(digits);
  double print_limit = pow(10.0, -digits) / 2;
  //wcout.precision(30);
  //double print_limit = 0.000001 / 2;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    double mass = (*i).mass;
    double rel_area = (*i).rel_area;
    double val_perc = rel_area / max_area * 100;
    //double val_norm = rel_area / sum_area;
    if(mass != DUMMY_MASS && val_perc >= print_limit)
      std::wcout << mass << L" " << val_perc << std::endl;
    //wcout << mass << L" " << val_perc << L" " << val_norm << endl;
  }
}

void IsotopeCalculator::print_to_output(Spectrum & output, Pattern & result)
{
  // find the maximum
  double max_area = 0;
  double sum_area = 0;
  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    if(max_area < (*i).rel_area)
      max_area = (*i).rel_area;
    sum_area += (*i).rel_area;
  }
  if(max_area == 0)
    return; // empty pattern

  for(Pattern::iterator i = result.begin(); i != result.end(); ++i) {
    double mass = (*i).mass;
    double rel_area = (*i).rel_area;
    double val_perc = rel_area / max_area * 100;
	if(mass != DUMMY_MASS && val_perc >= intensity_thresh)
		output.push_back( Peak(mass, rel_area) );
  }

  if( output.size() == 0 ) throw EmptyIsotopeSpectrumException();
  
}