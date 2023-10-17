/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Spectrum.h
#
# Description: 	Class for spectrum.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

#include <vector>
#include <algorithm>
#include <ostream>
#include "Config.h"

typedef std::pair<int,double> annotation_t; //<Fragment Id, Score>

static bool sort_annotations_by_score(const annotation_t &u, const annotation_t &v){
   return u.second > v.second;
}

class Peak{
public:
	Peak(){};
	Peak(double a_mass, double an_intensity) : 
	  mass(a_mass), intensity(an_intensity) {};
	double mass;
	double intensity;
	std::vector<annotation_t> annotations; 
};

static bool sort_peaks_by_intensity(const Peak &u, const Peak &v){
   return u.intensity > v.intensity;
}

static bool sort_peaks_by_mass(const Peak &u, const Peak &v){
   return u.mass < v.mass;
}

class Spectrum{
public:
	Spectrum() : is_normalized(false), is_sorted(false) {};

	//Iterating over the spectrum = iterating over the peaks
	typedef std::vector<Peak>::const_iterator const_iterator;
	const_iterator begin() const {return peaks.begin();};
	const_iterator end() const {return peaks.end();};
	void push_back( const Peak &pk ){ peaks.push_back(pk); is_normalized = false; is_sorted = false; };
	void clear() { peaks.clear(); peaks.resize(0); is_normalized = false; is_sorted = false; };
	unsigned int size() const { return peaks.size(); };
	const Peak *getPeak( int peak_idx ) const { return &peaks[peak_idx]; };
	void addPeakAnnotation( int peak_idx, annotation_t annot )
		{ peaks[peak_idx].annotations.push_back( annot ); };
	void updateAnnotationId( int peak_idx, int annot_idx, int new_id )
		{ peaks[peak_idx].annotations[annot_idx].first = new_id; }

	const std::vector<Peak> *getPeaks() const { return &peaks; };
	bool isNormalizedAndSorted() const { return is_normalized && is_sorted; };
	double getLastMass() const { return peaks.back().mass; };
	void postProcess( double perc_thresh, int min_peaks, int max_peaks );
	void normalizeAndSort();
	void roundPeaksToInteger();
	void setIntensityThreshold(double threshold = 0.1, int min_peaks = 5);
	void getHighPeaks(int num = 5);
	void clean(double abs_mass_tol, double ppm_mass_tol);	
	void sortAndNormalizeAnnotations();
	void outputToStream( std::ostream &out, bool do_annotate ) const;
	void outputToMspStream( std::ostream &out, std::string id, int ionization_mode, int energy ) const;
	void outputToMgfStream( std::ostream &out, std::string id, int ionization_mode, int energy, double mw ) const;
	void removePeaksWithNoFragment( std::vector<double> &frag_masses, double abs_tol, double ppm_tol );
	void quantisePeaksByMass( int num_dec_places );
    // brian modify
    void sortedByIntensity() {std::sort(peaks.begin(), peaks.end(), sort_peaks_by_intensity);};
    void sortedByMass() {std::sort(peaks.begin(), peaks.end(), sort_peaks_by_mass);};

private:
	std::vector<Peak> peaks;
	bool is_normalized;
	bool is_sorted;
};



#endif // __SPECTRUM_H__
