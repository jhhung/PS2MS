/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Isotope.h
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

#ifndef __ISOTOPE_H__
#define __ISOTOPE_H__

#include "Spectrum.h"
#include "Util.h"
#include <map>
#include <string>

class EmptyIsotopeSpectrumException: public std::exception{

	virtual const char* what() const throw(){
		return "Error in isotope spectrum generation: - no peaks generated above threshold.";
	}
};


struct ipeak
{
  double mass;
  double rel_area;
};

// map from element abbreviation to index in the elements table
typedef std::map<std::string, size_t> ElemMap;

// map from element index to the count of occurences in the formula
typedef std::map<size_t, long> FormMap;

typedef std::vector<ipeak> Pattern;              // index: peak_number
typedef std::vector<Pattern> SuperAtomList;        // index: bit_number 
typedef std::vector<SuperAtomList> SuperAtomData;  // index: element_number


//Exception to throw when something goes wrong during the isotope calculation
class IsotopeCalculationException: public std::exception{

	virtual const char* what() const throw(){
		return "Exception occurred computing isotope peaks";
	}
};

class IsotopeCalculator{

public:
	IsotopeCalculator( double a_intensity_thresh ) : verbose(0), intensity_thresh( a_intensity_thresh) { init_data("ISOTOPE.DAT"); };
	void computeIsotopeSpectrum( Spectrum &output, const romol_ptr_t mol, long charge );
	void setVerbose(){ verbose = true; };
	double getIntensityThresh() const { return intensity_thresh; };
private:

	double intensity_thresh;
	SuperAtomData sad;	//atom_idx -> isotope information
	ElemMap em;			//Atom symbol -> atom_idx

	//This function sets the formula map (atom_idx -> count) structure for the input molecule
	void setFormulaMap( FormMap &output, const romol_ptr_t mol );

	//The remainder of the functions are copied directly from emass (with minor mods)
	void init_data( const char *filename );
	void convolute_basic(Pattern & h, const Pattern & g, const Pattern & f);
	void prune(Pattern & f, double limit);
	void calculate(Pattern & tmp, Pattern & result, FormMap & fm, double limit, long charge);
	void print_pattern(Pattern & result, int digits);	

	//This function writes the output to our spectrum output
	void print_to_output(Spectrum & output, Pattern & result);

	bool verbose;
};


#endif // __ISOTOPE_H__
