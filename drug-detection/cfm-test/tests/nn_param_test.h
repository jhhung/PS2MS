/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# nn_param_test.h
#
# Description: 	Test code for NNParam.cpp
#
# Author: Felicity Allen
# Created: May 2015
#########################################################################*/

#ifndef __NN_PARAM_TEST_H__
#define __NN_PARAM_TEST_H__

#include "../test.h"
#include "NNParam.h"

class NNParamsTestComputeTransitionThetas : public Test {
public:
	NNParamsTestComputeTransitionThetas();
	void runTest();
};

class NNParamsTestSaveAndLoadFromFile : public Test {
public:
	NNParamsTestSaveAndLoadFromFile();
	void runTest();
};

class NNParamsTestComputeDeltas : public Test {
public:
	NNParamsTestComputeDeltas();
	void runTest();
};

class NNParamsTestComputeUnweightedGradients : public Test {
public:
	NNParamsTestComputeUnweightedGradients();
	void runTest();
};

class NNParamsTestComputeAndAccumulateGradient : public Test {
public:
	NNParamsTestComputeAndAccumulateGradient();
	void runTest();
};

class NNParamsTestBiasIndexes : public Test {
public:
	NNParamsTestBiasIndexes();
	void runTest();
};


#endif // __NN_PARAM_TEST_H__