/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# comparator_test.cpp
#
# Description: Test code for Comparators.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#include "mpi.h"

#include "comparator_test.h"
#include "MolData.h"

#include <boost/filesystem.hpp>

bool testComparator( Comparator &cmp, double exp_scores[3], double max_score ){

    double tol = 1e-8;
	bool pass = true;

	std::vector<double> scores(3);
	//Test with exact matches
	config_t cfg; initDefaultConfig(cfg);
	MolData moldata("Test1","C", &cfg);
	std::string specfile = "tests/test_data/test_spec/Test1.txt";
	moldata.readInSpectraFromFile(specfile);
	std::string pspecfile = "tests/test_data/test_pspec/Test1.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i));
	if( fabs( scores[0] - max_score ) > tol || fabs( scores[1] - max_score ) > tol || fabs( scores[2] - max_score ) > tol ){
		std::cout << "Perfect Matches: Expecting " << max_score << " scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test with no matches
	moldata = MolData("Test2","C", &cfg);
	specfile = "tests/test_data/test_spec/Test2.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test2.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i) );
	if( fabs( scores[0] - 0.0 ) > tol || fabs( scores[1] - 0.0 ) > tol || fabs( scores[2] - 0.0 ) > tol ){
		std::cout << "No Match: Expecting 0.0 scores but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	//Test between
	moldata = MolData("Test3","C", &cfg);
	specfile = "tests/test_data/test_spec/Test3.txt";
	moldata.readInSpectraFromFile(specfile);
	pspecfile = "tests/test_data/test_pspec/Test3.txt";
	moldata.readInSpectraFromFile(pspecfile, true);
	for( int i = 0; i < 3; i++ ) 
		scores[i] = cmp.computeScore( moldata.getSpectrum(i), moldata.getPredictedSpectrum(i));
	if( fabs( scores[0] - exp_scores[0] ) > tol || fabs( scores[1] - exp_scores[1] ) > tol || fabs( scores[2] - exp_scores[2] ) > tol ){
		std::cout << "Mixed: Expecting " << std::setprecision(10) << exp_scores[0] << " " << exp_scores[1] << " " << exp_scores[2] << " but found " << scores[0] << " " << scores[1] << " " << scores[2] << std::endl;
		pass = false;
	}

	return pass;

}



ComparatorsTestWeightedRecall::ComparatorsTestWeightedRecall(){
	description = "Test computing of weighted recall score";
}

void ComparatorsTestWeightedRecall::runTest(){
	
	WeightedRecall cmp(10.0, 0.01);
	double exp_scores[3] = {100.0, 65.0, 50.0};
	bool pass = testComparator( cmp, exp_scores, 100.0 );
	passed = pass;

}

ComparatorsTestRecall::ComparatorsTestRecall(){
	description = "Test computing of recall score";
}

void ComparatorsTestRecall::runTest(){
	
	Recall cmp(10.0, 0.01);
	double exp_scores[3] = {100.0, 66.6666666666666, 50.0};
	bool pass = testComparator( cmp, exp_scores, 100.0 );
	passed = pass;

}


ComparatorsTestPrecision::ComparatorsTestPrecision(){
	description = "Test computing of precision score";
}

void ComparatorsTestPrecision::runTest(){
	
	Precision cmp(10.0, 0.01);
	double exp_scores[3] = {50.0, 57.142857142, 50.0};
	bool pass = testComparator( cmp, exp_scores, 100.0 );
	passed = pass;

}

ComparatorsTestWeightedPrecision::ComparatorsTestWeightedPrecision(){
	description = "Test computing of weighted precision score";
}

void ComparatorsTestWeightedPrecision::runTest(){
	
	WeightedPrecision cmp(10.0, 0.01);
	double exp_scores[3] = {57.142857142, 70.0, 55.5555555555555};
	bool pass = testComparator( cmp, exp_scores, 100.0 );
	passed = pass;
}

ComparatorsTestJaccard::ComparatorsTestJaccard(){
	description = "Test computing of jaccard score";
}

void ComparatorsTestJaccard::runTest(){
	
	Jaccard cmp(10.0, 0.01);
	double exp_scores[3] = {0.666666667, 0.61538461538461542, 0.5};
	bool pass = testComparator( cmp, exp_scores, 1.0 );
	passed = pass;

}
