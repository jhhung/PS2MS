/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# param_test.h
#
# Description: 	Test code for Param.cpp
#
# Author: Felicity Allen
# Created: January 2013
#########################################################################*/

#ifndef __PARAM_TEST_H__
#define __PARAM_TEST_H__

#include "../test.h"
#include "Param.h"

class ParamsTestComputeTransitionThetas : public Test {
public:
	ParamsTestComputeTransitionThetas();
	void runTest();
};

class ParamsTestComputeAndAccumulateGradient : public Test {
public:
	ParamsTestComputeAndAccumulateGradient();
	void runTest();
};

#endif // __PARAM_TEST_H__