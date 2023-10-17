/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# em_test.h
#
# Description: 	Test code for EM.cpp
#
# Author: Felicity Allen
# Created: August 2013
#########################################################################*/

#ifndef __EM_TEST_H__
#define __EM_TEST_H__

#include "../test.h"

class EMTestSelfProduction : public Test {
public:
	EMTestSelfProduction();
	void runTest();
};

class EMTestSingleEnergySelfProduction : public Test {
public:
	EMTestSingleEnergySelfProduction();
	void runTest();
};

class EMTestNNSingleEnergySelfProduction : public Test {
public:
	EMTestNNSingleEnergySelfProduction();
	void runTest();
};

class EMTestSingleEnergyIsotopeSelfProduction : public Test {
public:
	EMTestSingleEnergyIsotopeSelfProduction();
	void runTest();
};

class EMTestLBFGSvsOriginalGradientAscent : public Test {
public:
	EMTestLBFGSvsOriginalGradientAscent();
	void runTest();
};

class EMTestMultiProcessor: public Test {
public:
	EMTestMultiProcessor();
	void runTest();
};

class EMTestMultiProcessorLBFGS: public Test {
public:
	EMTestMultiProcessorLBFGS();
	void runTest();
};

class EMTestBiasPreLearning: public Test {
public:
	EMTestBiasPreLearning();
	void runTest();
};

class EMTestMiniBatchSelection: public Test {
public:
	EMTestMiniBatchSelection();
	void runTest();
};

#endif // __EM_TEST_H__