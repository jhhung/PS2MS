/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# spectrum_clean_test.h
#
# Description: 	Test code for MolData::cleanSpectra
#
# Author: Felicity Allen
# Created: January 2013
#########################################################################*/

#include "spectrum_test.h"
#include "MolData.h"

#include <boost/filesystem.hpp>

SpectrumQuantiseTest::SpectrumQuantiseTest(){
	description = "Test quantising of spectrum";
}

void SpectrumQuantiseTest::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	config_t cfg; initDefaultConfig(cfg);
	MolData moldata("Test4","C", &cfg);
	std::string specfile = "tests/test_data/test_pspec/Test3.txt";

	//3 Decimal places
	moldata.readInSpectraFromFile(specfile, true);
	moldata.quantisePredictedSpectra(3);
	const Spectrum *spec = moldata.getPredictedSpectrum(1);
	if( spec->size() != 7 ){
		std::cout << "Unexpected number of peaks in 3 decimal place spectrum: expected 7 but found " << spec->size() << std::endl;
		pass = false;
	}
	double exp_masses3[7] = {10.0, 50.99, 60.0, 76.008, 76.31, 76.431, 76.432};
	double exp_intensities3[7] = {20.0, 20.0, 35.0, 5.0, 5.0, 5.0, 10.0};
	for( int i = 0; i < 7; i++){ 
		if( fabs( spec->getPeak(i)->mass - exp_masses3[i]) > tol ){
			std::cout << "Unexpected peak mass for 3 decimal place spectrum: " << spec->getPeak(i)->mass << " vs " << exp_masses3[i] << std::endl;
			pass = false;
		}
		if( fabs( spec->getPeak(i)->intensity - exp_intensities3[i]) > tol ){
			std::cout << "Unexpected peak intensity for 3 decimal place spectrum: " << spec->getPeak(i)->intensity << " vs " << exp_intensities3[i] << std::endl;
			pass = false;
		}
	}

	//2 Decimal places
	moldata.readInSpectraFromFile(specfile, true);
	moldata.quantisePredictedSpectra(2);
	spec = moldata.getPredictedSpectrum(1);
	if( spec->size() != 6 ){
		std::cout << "Unexpected number of peaks in 2 decimal place spectrum: expected 6 but found " << spec->size() << std::endl;
		pass = false;	
	}
	double exp_masses2[6] = {10.0, 50.99, 60.0, 76.01, 76.31, 76.43};
	double exp_intensities2[6] = {20.0, 20.0, 35.0, 5.0, 5.0, 15.0};
	for( int i = 0; i < 6; i++){ 
		if( fabs( spec->getPeak(i)->mass - exp_masses2[i]) > tol ){
			std::cout << "Unexpected peak mass for 3 decimal place spectrum: " << spec->getPeak(i)->mass << " vs " << exp_masses2[i] << std::endl;
			pass = false;
		}
		if( fabs( spec->getPeak(i)->intensity - exp_intensities2[i]) > tol ){
			std::cout << "Unexpected peak intensity for 3 decimal place spectrum: " << spec->getPeak(i)->intensity << " vs " << exp_intensities2[i] << std::endl;
			pass = false;
		}
	}

	//0 Decimal places
	moldata.readInSpectraFromFile(specfile, true);
	moldata.quantisePredictedSpectra(0);
	spec = moldata.getPredictedSpectrum(1);
	if( spec->size() != 4 ){
		std::cout << "Unexpected number of peaks in 0 decimal place spectrum: expected 4 but found " << spec->size() << std::endl;
		pass = false;	
	}
	double exp_masses0[4] = {10.0, 51.0, 60.0, 76.0};
	double exp_intensities0[4] = {20.0, 20.0, 35.0, 25.0};
	for( int i = 0; i < 4; i++){ 
		if( fabs( spec->getPeak(i)->mass - exp_masses0[i]) > tol ){
			std::cout << "Unexpected peak mass for 3 decimal place spectrum: " << spec->getPeak(i)->mass << " vs " << exp_masses0[i] << std::endl;
			pass = false;
		}
		if( fabs( spec->getPeak(i)->intensity - exp_intensities0[i]) > tol ){
			std::cout << "Unexpected peak intensity for 3 decimal place spectrum: " << spec->getPeak(i)->intensity << " vs " << exp_intensities0[i] << std::endl;
			pass = false;
		}
	}

	passed = pass;
}


SpectrumCleanTest::SpectrumCleanTest(){
	description = "Test cleaning of spectrum";
}

void SpectrumCleanTest::runTest(){
	
	bool pass = true;
    double tol = 1e-10;

	//Test with exact matches
	config_t cfg; initDefaultConfig(cfg);
	MolData moldata("Test4","C", &cfg);
	std::string specfile = "tests/test_data/test_spec/Test4.txt";
	moldata.readInSpectraFromFile(specfile);
	moldata.cleanSpectra(0.1, 10.0);

	//Check that there is one peak over intensity 10.0 within 0.1 of 113.07 in each spectrum
	for( int energy = 0; energy < 3; energy++ ){
		const Spectrum *spec = moldata.getSpectrum(energy);
		Spectrum::const_iterator it = spec->begin();
		int count = 0;
		for( ; it != spec->end(); ++it )
			if( fabs(it->mass - 113.07) < 0.01 && it->intensity > 10.0 ) count++;
		if( count != 1 ){
			std::cout << "Incorrect count for 113.07 peak"	<< std::endl;
			pass = false;
		}
	}

	passed = pass;

}
