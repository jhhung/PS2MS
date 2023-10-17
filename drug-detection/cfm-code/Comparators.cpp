/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Comparators.cpp
#
# Description: 	Functions for comparing spectra and returning a score.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include <cmath> 
#include "Comparators.h"
#include <boost/math/distributions.hpp>

void Comparator::getMatchingPeakPairs( std::vector<peak_pair_t> &peak_pairs,  const Spectrum *p, const Spectrum *q ) const{
	
	//Check that the input spectra are sorted and normalized correctly
	if( !p->isNormalizedAndSorted() || !q->isNormalizedAndSorted() ) throw ComparatorException();

	//Use Dynamic Programming to find the best match between peaks 
	//such that each peak matches at most one other
	std::vector< std::vector<double> > dp_vals(p->size()+1);
	for( unsigned int i = 0; i <= p->size(); i++ ) dp_vals[i].resize(q->size()+1);
	for( unsigned int i = 0; i <= p->size(); i++ ) dp_vals[i][0] = 0.0;
	for( unsigned int j = 0; j <= q->size(); j++ ) dp_vals[0][j] = 0.0;

	Spectrum::const_iterator itp = p->begin();
	for( unsigned int i = 1; itp != p->end(); ++itp, i++ ){
	
		Spectrum::const_iterator itq = q->begin();
		for( unsigned int j = 1; itq != q->end(); ++itq, j++ ){

			//Compute the mass tolerance 
			double mass = itp->mass;
			if( itq->mass > mass ) mass = itq->mass;
			double mass_tol = (mass/1000000.0) * ppm_tol;
			if( mass_tol < abs_tol ) mass_tol = abs_tol;

			//Check for a match
			double best_val = 0.0;
			if( fabs( itp->mass - itq->mass ) < mass_tol ){
				best_val = dp_vals[i-1][j-1] + itp->intensity * itq->intensity + 1;	//Weight the match by the average intensity of the two peaks
			}

			//Check for better values without matching
			if( dp_vals[i-1][j] > best_val ) best_val = dp_vals[i-1][j];
			if( dp_vals[i][j-1] > best_val ) best_val = dp_vals[i][j-1];

			//Store the best result up to this point
			dp_vals[i][j] = best_val;
		}
	}
	// std::cerr << std::endl;
	// for (int i = 0; i <= p->size(); ++i)
	// {
	// 	for (int j = 0; j <= q->size(); ++j)
	// 		std::cerr << dp_vals[i][j] << " ";
	// 	std::cerr << std::endl;
	// }
			

	//Backtrack to collect the matching peaks
	unsigned int i = p->size();
	unsigned int j = q->size();
	while( i > 0 && j > 0 ){
		
		if( fabs(dp_vals[i][j] - dp_vals[i][j-1]) < 1e-12 ) j--;
		else if( fabs(dp_vals[i][j] - dp_vals[i-1][j]) < 1e-12 ) i--;
		else if( dp_vals[i][j] > dp_vals[i-1][j-1] + 1e-12 ){
			peak_pairs.push_back(peak_pair_t());
			peak_pairs.back().first = *p->getPeak(i-1);
			peak_pairs.back().second = *q->getPeak(j-1);			
			i--; j--;
		}
		else {
			i--; j--;
		}
	}

}

double DotProduct::computeScore( const Spectrum *measured, const Spectrum *predicted ) const{

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted );

	if( peak_pairs.size() == 0 ) return 0.0;

	std::vector<peak_pair_t>::iterator it = peak_pairs.begin();
	double num = 0.0, denomp = 0.0, denomq = 0.0;
	for( ; it != peak_pairs.end(); ++it ){
		double p_int = getAdjustedIntensity(it->first.intensity, it->first.mass );
		double q_int = getAdjustedIntensity(it->second.intensity, it->second.mass );
		num += p_int*q_int;
	}
	// std::cerr << "NUM: " << num << std::endl;
	denomp = getTotalPeakSum( measured );
	denomq = getTotalPeakSum( predicted );
	// std::cerr << "denomp: " << denomp << std::endl;
	// std::cerr << "denomq: " << denomq << std::endl;
	// std::cerr << "SCORE: " << num*num/(denomp*denomq) << std::endl;
	return num*num/(denomp*denomq);
}

double DotProduct::getAdjustedIntensity( double intensity, double mass ) const{
	return std::pow(intensity, 0.5)*std::pow(mass, 0.5);
}

double OrigSteinDotProduct::getAdjustedIntensity( double intensity, double mass ) const{
	return std::pow(intensity, 0.6)*std::pow(mass, 3);
}

double DotProduct::getTotalPeakSum( const Spectrum *spectrum ) const{

	double result = 0.0;
	Spectrum::const_iterator it = spectrum->begin();
	for( ; it != spectrum->end(); ++it ){
		double pint = getAdjustedIntensity(it->intensity, it->mass);
		result += pint*pint;
	}
	return result;
}

double Precision::computeScore( const Spectrum *measured, const Spectrum *predicted ) const {

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted );

	if( predicted->size() == 0 ) return 0.0;
	return (double)peak_pairs.size() * 100.0 / (double)predicted->size();

}

double Recall::computeScore( const Spectrum *measured, const Spectrum *predicted ) const {

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted );

	if( measured->size() == 0 ) return 0.0;
	return (double)peak_pairs.size() * 100.0 / (double)measured->size();

}

double WeightedRecall::computeScore( const Spectrum *measured, const Spectrum *predicted ) const {

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted);

	if( measured->size() == 0 ) return 0.0;
	
	//Numerator
	double num = 0.0;
	std::vector<peak_pair_t>::iterator it = peak_pairs.begin();
	for( ; it != peak_pairs.end(); ++it )
		num += it->first.intensity;

	//Denominator
	Spectrum::const_iterator itc = measured->begin();
	double den = 0.0;
	for( ; itc != measured->end(); ++itc )
		den += itc->intensity;

	return num*100.0/den;
}

double WeightedPrecision::computeScore( const Spectrum *measured, const Spectrum *predicted ) const{

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted);

	if( predicted->size() == 0 ) return 0.0;
	
	//Numerator
	double num = 0.0;
	std::vector<peak_pair_t>::iterator it = peak_pairs.begin();
	for( ; it != peak_pairs.end(); ++it )
		num += it->second.intensity;

	//Denominator
	Spectrum::const_iterator itc = predicted->begin();
	double den = 0.0;
	for( ; itc != predicted->end(); ++itc )
		den += itc->intensity;

	return num*100.0/den;
}

double WeightedJaccard::computeScore( const Spectrum *measured, const Spectrum *predicted ) const{

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted);
	double num = 0.0;
	std::vector<peak_pair_t>::iterator it = peak_pairs.begin();
	for( ; it != peak_pairs.end(); ++it ){
		num += it->first.intensity;
		num += it->second.intensity;
	}
	return num/200.0;
}


double Jaccard::computeScore( const Spectrum *measured, const Spectrum *predicted ) const{

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted);
	return 2*(double)peak_pairs.size()/(measured->size() + predicted->size());

}

double Combined::computeScore( const Spectrum *measured, const Spectrum *predicted ) const{

	DotProduct dp( ppm_tol, abs_tol );
	WeightedRecall wr( ppm_tol, abs_tol );
	WeightedPrecision wp( ppm_tol, abs_tol );
	Jaccard ja( ppm_tol, abs_tol );
	double dp_score = 100.0 * dp.computeScore(measured, predicted);
	double ja_score = 100.0 * ja.computeScore(measured, predicted);
	double wp_score = wp.computeScore(measured, predicted);
	double wr_score = wr.computeScore(measured, predicted);
	std::cout << dp_score << " " << ja_score << " " << wp_score << " " << wr_score;
	return dp_score + ja_score + wp_score + wr_score;

}

void ScatterOutput::outputData( const Spectrum *measured, const Spectrum *predicted ){

	std::vector<peak_pair_t> peak_pairs;
	getMatchingPeakPairs( peak_pairs, measured, predicted);
	
	//Write the pairs to file
	std::vector<peak_pair_t>::iterator it = peak_pairs.begin();
	for( ; it != peak_pairs.end(); ++it ){
		out << it->first.mass << " " << it->first.intensity << " ";
		out << it->second.mass << " " << it->second.intensity << std::endl;
	}
}
