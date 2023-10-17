/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# ipfp_test.h
#
# Description: 	Test code for IPFP.cpp
#
# Author: Felicity Allen
# Created: January 2013
#########################################################################*/

#ifndef __IPFP_TEST_H__
#define __IPFP_TEST_H__

#include "../test.h"
#include "IPFP.h"

class IPFPTestComputeBeliefsConverge : public Test {
public:
	IPFPTestComputeBeliefsConverge();
	void runTest();
};

class IPFPTestComputeBeliefsNonConverge : public Test {
public:
	IPFPTestComputeBeliefsNonConverge();
	void runTest();
};

class IPFPTestComputeBeliefsSharedMass : public Test {
public:
	IPFPTestComputeBeliefsSharedMass();
	void runTest();
};

class IPFPTestComputeBeliefsSingleEnergy : public Test {
public:
	IPFPTestComputeBeliefsSingleEnergy();
	void runTest();
};

#endif // __IPFP_TEST_H__