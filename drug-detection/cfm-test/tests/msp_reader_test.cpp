/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# msp_reader_test.h
#
# Description: 	Test code for MspReader
#
# Author: Felicity Allen
# Created: August 2014
#########################################################################*/

#include "msp_reader_test.h"
#include "MolData.h"
#include "MspReader.h"

#include <boost/filesystem.hpp>

MspReaderTest::MspReaderTest(){
	description = "Test MSP spectrum reader";
}

void MspReaderTest::runTest(){
	
	bool pass = true;
	double tol = 1e-5;

	MspReader msp = MspReader("tests/test_data/nist2011_cutdown.msp", "NIST2011_");

	std::cout << std::setprecision(8);

	//Test first spectrum in msp
	config_t cfg; initDefaultConfig(cfg);
	MolData mol1( "NIST2011_1", "InChI=1S/H2/h1H", &cfg );
	mol1.readInSpectraFromMSP( msp );

	if( mol1.getNumSpectra() != 1 ){
		std::cout << "Unexpected number of spectra for NIST2011_1: " << mol1.getNumSpectra() <<  std::endl;
		pass = false;
	}
	else{
		const Spectrum *spec = mol1.getSpectrum(0);
		if( spec->size() != 2 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 2 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 1.0) > tol || fabs(spec->getPeak(0)->intensity - 2.056903) > tol ||
				fabs(spec->getPeak(1)->mass - 2.0) > tol || fabs(spec->getPeak(1)->intensity - 97.9430962) > tol ){		
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				std::cout << "Expecting  1 20.98, 2 97.9431" << std::endl;
				std::cout << "Found " << spec->getPeak(0)->mass << " " << spec->getPeak(0)->intensity << ", " << spec->getPeak(1)->mass << " " << spec->getPeak(1)->intensity << std::endl;
				pass = false;		
			}
		}
	}
	
	//Test mid spectrum in msp
	MolData mol2( "NIST2011_71459", "InChI=1S/H2/h1H", &cfg );
	mol2.readInSpectraFromMSP( msp );

	if( mol2.getNumSpectra() != 1 ){
		std::cout << "Unexpected number of spectra for NIST2011_71459: " << mol2.getNumSpectra() <<  std::endl;
		pass = false;
	}
	else{
		const Spectrum *spec = mol2.getSpectrum(0);
		if( spec->size() != 166 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 166 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 33.0) > tol || fabs(spec->getPeak(0)->intensity - 0.010019287) > tol ||
				fabs(spec->getPeak(165)->mass - 469.0) > tol || fabs(spec->getPeak(165)->intensity - 0.0050096436) > tol ){
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				std::cout << "Expecting 33 0.010019287, ...., 469 0.0050096436" << std::endl;
				std::cout << "Found " << spec->getPeak(0)->mass << " " << spec->getPeak(0)->intensity << ",...., " << spec->getPeak(165)->mass << " " << spec->getPeak(165)->intensity << std::endl;
				pass = false;		
			}
		}
	}

	//Test last spectrum in msp
	MolData mol3( "NIST2011_212964", "InChI=1S/H2/h1H", &cfg );
	mol3.readInSpectraFromMSP( msp );

	if( mol3.getNumSpectra() != 1 ){
		std::cout << "Unexpected number of spectra for NIST2011_212964: " << mol3.getNumSpectra() <<  std::endl;
		pass = false;
	}
	else{
		const Spectrum *spec = mol3.getSpectrum(0);
		if( spec->size() != 30 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 30 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 192.0) > tol || fabs(spec->getPeak(0)->intensity -  0.87150087) > tol ||
				fabs(spec->getPeak(29)->mass - 1168.0) > tol || fabs(spec->getPeak(29)->intensity - 18.15324) > tol ){
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				std::cout << "Expecting  192 0.87150087, ...., 1168.0 18.15324" << std::endl;
				std::cout << "Found " << spec->getPeak(0)->mass << " " << spec->getPeak(0)->intensity << ",...., " << spec->getPeak(29)->mass << " " << spec->getPeak(29)->intensity << std::endl;
				pass = false;		
			}
		}
	}


	passed = pass;

}

MspReaderMultipleEnergiesTest::MspReaderMultipleEnergiesTest(){
	description = "Test MSP spectrum reader with multiple energy levels per molecule";
}

void MspReaderMultipleEnergiesTest::runTest(){
	
	bool pass = true;
	double tol = 1e-5;

	MspReader msp = MspReader("tests/test_data/three_energies.msp", "");

	std::cout << std::setprecision(8);

	//Test first spectrum in msp
	config_t cfg; initDefaultConfig(cfg);
	MolData mol1( "Test3", "InChI=1S/H2/h1H", &cfg );
	mol1.readInSpectraFromMSP( msp );

	if( mol1.getNumSpectra() != 3 ){
		std::cout << "Unexpected number of spectra for Test1: " << mol1.getNumSpectra() <<  std::endl;
		pass = false;
	}
	else{
		const Spectrum *spec = mol1.getSpectrum(0);
		if( spec->size() != 2 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 2 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 54.07127368) > tol || fabs(spec->getPeak(0)->intensity - 10.0) > tol ||
				fabs(spec->getPeak(1)->mass - 76.08692374) > tol || fabs(spec->getPeak(1)->intensity - 90.0) > tol ){		
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				pass = false;		
			}
		}
		spec = mol1.getSpectrum(1);
		if( spec->size() != 2 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 2 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 15.02292652) > tol || fabs(spec->getPeak(0)->intensity - 40.0) > tol ||
				fabs(spec->getPeak(1)->mass - 47.06037464) > tol || fabs(spec->getPeak(1)->intensity - 60.0) > tol ){		
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				pass = false;		
			}
		}
		spec = mol1.getSpectrum(2);
		if( spec->size() != 2 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 2 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 15.02292652) > tol || fabs(spec->getPeak(0)->intensity - 20.0) > tol ||
				fabs(spec->getPeak(1)->mass - 29.01342445) > tol || fabs(spec->getPeak(1)->intensity - 80.0) > tol ){		
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				pass = false;		
			}
		}
	}

	//Test last spectrum in msp
	MolData mol3( "Test5", "InChI=1S/H2/h1H", &cfg );
	mol3.readInSpectraFromMSP( msp );

	if( mol3.getNumSpectra() != 3 ){
		std::cout << "Unexpected number of spectra for Test5: " << mol3.getNumSpectra() <<  std::endl;
		pass = false;
	}
	else{
		const Spectrum *spec = mol3.getSpectrum(0);
		if( spec->size() != 2 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 2 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 241.1798211) > tol || fabs(spec->getPeak(0)->intensity -  50.0) > tol ||
				fabs(spec->getPeak(1)->mass - 315.2166005) > tol || fabs(spec->getPeak(1)->intensity - 50.0) > tol ){
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				pass = false;		
			}
		}
		spec = mol3.getSpectrum(1);
		if( spec->size() != 3 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 3 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 299.1853004) > tol || fabs(spec->getPeak(0)->intensity -  25.0) > tol ||
				fabs(spec->getPeak(2)->mass - 315.2166005) > tol || fabs(spec->getPeak(2)->intensity - 25.0) > tol ){
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				pass = false;		
			}
		}
		spec = mol3.getSpectrum(2);
		if( spec->size() != 1 ){
			std::cout << "Unexpected number of peaks in spectrum: expecting 1 but found " << spec->size() << std::endl;
			pass = false;
		}
		else{
			if( fabs(spec->getPeak(0)->mass - 29.03857658) > tol || fabs(spec->getPeak(0)->intensity -  100.0) > tol){
				std::cout << "Unexpected peaks in spectrum:" << std::endl;
				pass = false;		
			}
		}
	}

	passed = pass;
}