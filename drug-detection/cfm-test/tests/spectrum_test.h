/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# spectrum_test.h
#
# Description: 	Test code for MolData::cleanSpectra
#
# Author: Felicity Allen
# Created: January 2013
#########################################################################*/

#ifndef __SPECTRUM_TEST_H__
#define __SPECTRUM_TEST_H__

#include "../test.h"

class SpectrumQuantiseTest : public Test {
public:
	SpectrumQuantiseTest();
	void runTest();
};

class SpectrumCleanTest : public Test {
public:
	SpectrumCleanTest();
	void runTest();
};

#endif // __SPECTRUM_TEST_H__