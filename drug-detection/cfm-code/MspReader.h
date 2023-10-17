/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# MspReader.h
#
# Description: 	Class for handling spectrum input in .msp format
#
# Copyright (c) 2014, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MSP_READER_H__
#define __MSP_READER_H__

#include "Spectrum.h"
#include <map>
#include <string>

class MspReadException: public std::exception{

	virtual const char* what() const throw(){
		return "Could not read msp file.";
	}
};

class MspIdException: public std::exception{

	virtual const char* what() const throw(){
		return "Could not find id in msp file.";
	}
};

class MspReader{
public:
	//Constructor
	MspReader(){};
	MspReader( const char *filename, const char *pre_id = "", bool normalize_and_sort = false, int quantise_dec_pl = -1){ readInMspFile( filename, pre_id, normalize_and_sort, quantise_dec_pl ); };

	//Function to fetch a spectrum with a given id from the msp library
	const std::vector<Spectrum> *fetchSpectrumForId( const char *id  );
	void writeLibraryToMspFile( const char *filename, int ionization_mode = DEFAULT_IONIZATION_MODE ) const;

private:
	void readInMspFile( const char *filename, const char *pre_id, bool normalize_and_sort, int quantise_dec_pl);

	//Library to read all the spectra into - accessible by id
	std::map< std::string, std::vector<Spectrum> > library;
	std::string mspfilename;
};

#endif // __MSP_READER_H__
