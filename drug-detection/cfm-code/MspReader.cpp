/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Config.cpp
#
# Description: 	Structs, functions, defaults for setting general 
#			    configuration data.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "MspReader.h"
#include "Spectrum.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>


const std::vector<Spectrum> *MspReader::fetchSpectrumForId( const char *id  ) { 
	
	if( library.find(id) == library.end() ){
		std::cout << "Could not find spectrum for " << id << " in MSP file " << mspfilename << std::endl;
		throw MspIdException();		
	}
	return &library[id]; 
};

void MspReader::readInMspFile( const char *filename, const char *pre_id, bool normalize_and_sort, int quantise_dec_pl ){

	std::string line;
	std::ifstream ifs ( filename, std::ifstream::in  );
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(";:");


	int numpeaks = 0;
	int lineno = 0;
	double mass = 0.0, intensity = 0.0;
	std::string id = "";
	int num_added = 0;

	//Read the config file into the paramater update config structure
	int energy = 0;
	if(!ifs) std::cout << "Could not open file " << filename << std::endl;
	while( ifs.good() ){

		getline( ifs, line );
		lineno++;
		if( line.size() < 3 ){ 
			numpeaks = 0; id = "";	//Empty line, so reset
			continue;	
		}
	
		if( numpeaks > 0 ){	
			if( line.substr( 0, 3 ) == "ID:" || line.substr( 0, 5 ) == "Name:"){
				//Note: this will only work if there are less peaks than stated
				std::cout << "Warning: Number of peaks read did not match stated number of peaks for " << id << std::endl;
				numpeaks = 0;
			}
			else{
				std::stringstream ss1(line);
				tokenizer tokens(line, sep);
				tokenizer::iterator tok_iter = tokens.begin();
				for ( ; tok_iter != tokens.end(); ++tok_iter){
					if( tok_iter->size() < 2 ) continue;
					std::stringstream ss1(*tok_iter);	
					ss1 >> mass >> intensity;
					library[id][energy].push_back( Peak(mass, intensity) );
					if( numpeaks == 1 ){ 
						if(normalize_and_sort) 
							library[id][energy].normalizeAndSort();
						if(quantise_dec_pl >= 0) 
							library[id][energy].quantisePeaksByMass(quantise_dec_pl);
					}
					numpeaks--;
				}
			}
		}

		if( boost::to_upper_copy(line.substr( 0, 3 )) == "ID:" ){
			std::stringstream ss1(line);
			std::string tag;
			ss1 >> tag >> id;
			id = pre_id + id;
		}

		if( line.substr( 0, 15 ) == "Comment: Energy" ){
			std::stringstream ss1(line);
			int tmp = line.size();
			energy = atoi( line.substr(15,line.size()-15).c_str() );
		}

		if( boost::to_upper_copy(line.substr( 0, 10 )) == "NUM PEAKS:" ){
			tokenizer tokens(line, sep);
			tokenizer::iterator tok_iter = tokens.begin();
			if( ++tok_iter != tokens.end() )
				numpeaks = atoi(tok_iter->c_str());
			if( id == "" ){
				std::cout << "Missing ID field for entry at line " << lineno << std::endl;
				throw MspReadException();
			}
			library[id];	//initialise the spectrum for this id
			if( library[id].size() < energy + 1)
				library[id].resize(energy + 1);
			num_added += 1;
		}
	}
	std::cout << "Found " << num_added << " molecules in MSP file " << filename << std::endl;
}


void MspReader::writeLibraryToMspFile( const char *filename, int ionization_mode ) const{

	std::ofstream out;
	out.open(filename, std::fstream::out | std::fstream::app);

	std::map< std::string, std::vector<Spectrum> >::const_iterator it = library.begin();
	for( ; it != library.end(); ++it ){
		std::vector<Spectrum>::const_iterator its = it->second.begin();
		for( int energy = 0; its != it->second.end(); ++its, energy++ )
			its->outputToMspStream( out, it->first, ionization_mode, energy );
	}
	out.close();

}