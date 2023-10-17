/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Comparators.h
#
# Description: 	Functions for comparing spectra and returning a score.
#				So far, just includes the dot product as described in Stein 1994.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __COMPARATORS_H__
#define __COMPARATORS_H__

#include "MolData.h"

#include <iostream>
#include <sstream>
#include <fstream>

class ComparatorException: public std::exception{

	virtual const char* what() const throw(){
		return "Unsorted or unnormalized spectra detected in comparator - cannot proceed.";
	}
};

typedef std::pair<Peak, Peak> peak_pair_t;

//Base class to compute a spectrum comparison score
class Comparator{
public:
	Comparator( double a_ppm_tol, double a_abs_tol ) : ppm_tol(a_ppm_tol), abs_tol( a_abs_tol ) {};
	virtual double computeScore( const Spectrum *measured, const Spectrum *predicted ) const = 0;
	double getPPMTol() const {return ppm_tol; };
	double getAbsTol() const {return abs_tol; };
protected:
	double ppm_tol;
	double abs_tol;

	void getMatchingPeakPairs( std::vector<peak_pair_t> &peak_pairs,  const Spectrum *p, const Spectrum *q ) const;
};

class DotProduct : public Comparator {
public:
	DotProduct( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted) const;
private:
	double getTotalPeakSum( const Spectrum *spectrum ) const;
	virtual double getAdjustedIntensity( double intensity, double mass ) const;
};

class OrigSteinDotProduct : public DotProduct {
public:
	OrigSteinDotProduct( double a_ppm_tol, double a_abs_tol ) : DotProduct( a_ppm_tol, a_abs_tol ) {};
private:
	double getAdjustedIntensity( double intensity, double mass ) const;
};

class Precision : public Comparator {
public:
	Precision( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class WeightedRecall : public Comparator {
public:
	WeightedRecall( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class Recall : public Comparator {
public:
	Recall( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};


class WeightedPrecision: public Comparator {
public:
	WeightedPrecision( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class Jaccard : public Comparator {
public:
	Jaccard( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class WeightedJaccard : public Comparator {
public:
	WeightedJaccard( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class Combined : public Comparator {
public:
	Combined( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const;
private:
};

class ScatterOutput : public Comparator {
public:
	ScatterOutput( double a_ppm_tol, double a_abs_tol ) : Comparator( a_ppm_tol, a_abs_tol ) {
		out.open("matching_peak_pairs.log");
		if( !out.is_open() )
			std::cout << "Warning: Trouble opening output peak pair file" << std::endl;
	};
	double computeScore( const Spectrum *measured, const Spectrum *predicted ) const {return -1.0;};
	void outputData( const Spectrum *measured, const Spectrum *predicted );
	~ScatterOutput(){out.close();};
	std::ofstream out;
private:

};

#endif // __COMPARATORS_H__