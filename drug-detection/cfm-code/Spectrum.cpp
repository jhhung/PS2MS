/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Spectrum.cpp
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

#include "Util.h"
#include "Spectrum.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>

void Spectrum::outputToStream( std::ostream &out, bool do_annotate ) const{
	
	std::vector<Peak>::const_iterator itp = peaks.begin();
	for( ; itp != peaks.end(); ++itp ){
		out << std::setprecision(10) << itp->mass << " " << itp->intensity;
		if(do_annotate){
			std::stringstream ss_values;
			ss_values << std::setprecision(5) << "(";
			std::vector<annotation_t>::const_iterator ita = itp->annotations.begin();
			for( ; ita != itp->annotations.end(); ++ita ){
				out << " " << ita->first;
				if( ita != itp->annotations.begin() ) ss_values << " ";
				ss_values << ita->second*100.0;
			}
			ss_values << ")";
			if( itp->annotations.size() > 0 )
				out << " " << ss_values.str();
		}
		out << std::endl;
	}
}

void Spectrum::outputToMspStream( std::ostream &out, std::string id, int ionization_mode, int energy ) const{

	if(  ionization_mode == POSITIVE_EI_IONIZATION_MODE )
		out << "Name: +ve in-silico MS by CFM-ID for " << id << std::endl;
	else if(  ionization_mode == POSITIVE_ESI_IONIZATION_MODE ) 
		out << "Name: +ve in-silico MS/MS by CFM-ID for " << id << std::endl;
	else
		out << "Name: -ve in-silico MS/MS by CFM-ID for " << id << std::endl;
	out << "ID: " << id << std::endl;
	out << "Comment: Energy" << energy << std::endl;
	out << "Num peaks: " << peaks.size() << std::endl; 
	outputToStream( out, false );
	out << std::endl;
}

void Spectrum::outputToMgfStream( std::ostream &out, std::string id, int ionization_mode, int energy, double mw ) const{

	out << "BEGIN IONS" << std::endl;
	out << "PEPMASS=" << std::setprecision(10) << mw << std::endl;
	if( ionization_mode == POSITIVE_ESI_IONIZATION_MODE || ionization_mode == POSITIVE_EI_IONIZATION_MODE )
		out << "CHARGE=1+" <<  std::endl;
	else if( ionization_mode == NEGATIVE_ESI_IONIZATION_MODE )
		out << "CHARGE=1-" <<  std::endl;
	out << "TITLE=" << id << ";Energy" << energy << ";";
	if( ionization_mode == POSITIVE_ESI_IONIZATION_MODE ) out << "[M+H]+;In-silico MS/MS by CFM-ID;" << std::endl;
	else if( ionization_mode == POSITIVE_ESI_IONIZATION_MODE ) out << "[M]+;In-silico MS by CFM-ID;" << std::endl;
	else if( ionization_mode == NEGATIVE_ESI_IONIZATION_MODE ) out << "[M-H]+;In-silico MS/MS by CFM-ID;" << std::endl;
	outputToStream( out, false );
	out << "END IONS" << std::endl;
}

void Spectrum::quantisePeaksByMass( int num_dec_places ){

	//Combine peaks that have the same mass when reduced to num_dec_places
	//Note: this is mostly used for extreme cases like the NIST data, where masses are given only to integer precision.
	normalizeAndSort();
	long long prev_mass = 0;
	std::vector<Peak>::iterator it = peaks.begin();
	for( ; it != peaks.end(); ++it ){
		long long tmp_mass = (long long)(it->mass * std::pow(10.0, num_dec_places) + 0.5);
		it->mass = tmp_mass*std::pow(10.0, -num_dec_places);
		if( tmp_mass == prev_mass ){ 
			it->intensity += (it-1)->intensity;
			it->annotations.insert( it->annotations.end(), (it-1)->annotations.begin(), (it-1)->annotations.end() );
			it = peaks.erase(it-1);
		}
		prev_mass = tmp_mass;
	}
	normalizeAndSort();

}


void Spectrum::postProcess( double perc_thresh, int min_peaks, int max_peaks ){

	std::sort( peaks.begin(), peaks.end(), sort_peaks_by_intensity );
	double total = 0.0;
	std::vector<Peak>::iterator it = peaks.begin();
	int count = 0;
	for( ; it != peaks.end(); ++it ){
		total += it->intensity;
		count++;
		//e.g. Take the top 80% of energy (assuming at least 5 peaks), 
		//or the highest 30 peaks (whichever comes first)
		if( (total > perc_thresh && count > min_peaks ) || count > max_peaks ) break;	
	}
	peaks.resize( count );
	std::sort( peaks.begin(), peaks.end(), sort_peaks_by_mass );
}

void Spectrum::normalizeAndSort(){
		
	if( !is_normalized ){
		//Compute the normalizer
		double sum = 0.0;
		std::vector<Peak>::iterator itp = peaks.begin();
		for( ; itp != peaks.end(); ++itp )
			sum += itp->intensity;
		double norm = 1.0;
		if( sum > 0.0 ) norm = 100.0/sum;

		//Adjust the values
		for( itp = peaks.begin(); itp != peaks.end(); ++itp )
			itp->intensity *= norm;
	}

	//Ensure the peaks are sorted by mass
	if( !is_sorted )
		std::sort( peaks.begin(), peaks.end(), sort_peaks_by_mass );

	is_sorted = true;
	is_normalized = true;
}

void Spectrum::roundPeaksToInteger()
{
	if( !is_sorted )
		std::sort( peaks.begin(), peaks.end(), sort_peaks_by_mass );
	std::vector<decltype(peaks)::iterator> remove_list;
	decltype(peaks)::iterator standard = peaks.begin();
	standard->mass = round(standard->mass);
	double round_mass;
	for (auto it = peaks.begin() + 1; it != peaks.end(); ++it)
	{
		round_mass = round(it->mass);
		if (round_mass != standard->mass)
		{
			standard = it;
			standard->mass = round_mass;
		}
		else
		{
			standard->intensity += it->intensity;
			remove_list.emplace_back(it);
		}
	}
	// remove merged peaks
	for (auto it = remove_list.rbegin(); it != remove_list.rend(); ++it)
		peaks.erase(*it);
	// re-normalize
	is_normalized = false;
	normalizeAndSort();
}

void Spectrum::setIntensityThreshold(double threshold, int min_peaks)
{
	std::vector<Peak> result;
	std::sort(peaks.begin(), peaks.end(), sort_peaks_by_intensity);
	double thresh = peaks.front().intensity * threshold;
	for (unsigned int i = 0; i < peaks.size(); ++i)
	{
		result.emplace_back(peaks[i]);
		if (result[i].intensity < thresh && i + 1 >= min_peaks)
			break;
	}
	peaks = result;
	if (is_sorted)
		std::sort(peaks.begin(), peaks.end(), sort_peaks_by_mass);
}

void Spectrum::getHighPeaks(int num)
{
	std::sort(peaks.begin(), peaks.end(), sort_peaks_by_intensity);
	peaks.erase(peaks.begin() + num, peaks.end());
	if (is_sorted)
		std::sort(peaks.begin(), peaks.end(), sort_peaks_by_mass);
}

void Spectrum::clean( double abs_mass_tol, double ppm_mass_tol ){

	//Ensure the initial spectrum is normalized and sorted by mass
	normalizeAndSort();

	//Filter peaks that are too close together
	// - Remove peaks within abs_mass_tol of a larger peak
	//	 and any peaks below an absolute intensity threshold
	double abs_intensity_thresh = 0.01;
	std::vector<bool> peak_flags( peaks.size(), true );
	
	//Forward Pass
	double prev_intensity = -1.0;
	double prev_mass = -100.0;
	std::vector<Peak>::iterator itp = peaks.begin();
	for( int idx=0; itp != peaks.end(); ++itp, idx++ ){
		if( itp->intensity < abs_intensity_thresh ) peak_flags[idx] = false;
		
		double mass_tol = getMassTol( abs_mass_tol, ppm_mass_tol, itp->mass );
		if( fabs( itp->mass - prev_mass ) < mass_tol &&
			itp->intensity < prev_intensity )
			peak_flags[idx] = false;
		prev_mass = itp->mass;
		prev_intensity = itp->intensity;
	}

	//Reverse Pass
	prev_intensity = -1.0;
	prev_mass = -100.0;
	std::vector<Peak>::reverse_iterator ritp = peaks.rbegin();
	for( int idx=peaks.size()-1; ritp != peaks.rend(); ++ritp, idx-- ){
		double mass_tol = getMassTol( abs_mass_tol, ppm_mass_tol, ritp->mass );
		if( fabs( ritp->mass - prev_mass ) < mass_tol &&
			ritp->intensity < prev_intensity )
			peak_flags[idx] = false;
		prev_mass = ritp->mass;
		prev_intensity = ritp->intensity;
	}

	//Alter the spectrum to include the selected peaks
	std::vector<Peak> peaks_copy(peaks);
	peaks.clear();
	itp = peaks_copy.begin();
	for( int idx=0; itp != peaks_copy.end(); ++itp, idx++ )
		if( peak_flags[idx] ) peaks.push_back( *itp );

	//Re-normalize
	normalizeAndSort();

}

void Spectrum::sortAndNormalizeAnnotations(){

	std::vector<Peak>::iterator it = peaks.begin();
	for( ; it != peaks.end();  ++it ){
		//Sort
		std::sort( it->annotations.begin(), it->annotations.end(), sort_annotations_by_score);
		
		//Normalize
		std::vector<annotation_t>::iterator itt = it->annotations.begin();
		double total = 0.0;
		for( ; itt != it->annotations.end(); ++itt ) total += itt->second;
		double norm = 1.0;
		if( total > 0.0) norm = it->intensity/(100.0*total);
		for( itt = it->annotations.begin(); itt != it->annotations.end(); ++itt ) 
			itt->second *= norm;
	}
}

void Spectrum::removePeaksWithNoFragment( std::vector<double> &frag_masses, double abs_tol, double ppm_tol ){
	
	//Remove any peaks more than mass_tol away from any fragment
	std::vector<Peak>::iterator itp = peaks.begin();
	for( ; itp != peaks.end(); ){
			
		double mass_tol = getMassTol( abs_tol, ppm_tol, itp->mass );
		bool found = false;
		std::vector<double>::iterator it = frag_masses.begin();
		for( ; it != frag_masses.end(); ++it ){
			if( fabs( *it - itp->mass ) < mass_tol ){ 
				found = true;
				break;
			}
		}
		if(!found) itp = peaks.erase( itp );
		else ++itp;
	}

	//Renormalise
	normalizeAndSort();

}