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

#ifndef __INFER_TEST_H__
#define __INFER_TEST_H__

#include "../test.h"
#include "Inference.h"

class InferenceTestCompareVsIPFPSingleEnergyCase : public Test {
public:
	InferenceTestCompareVsIPFPSingleEnergyCase();
	void runTest();
};

class InferenceTestCompareVsIPFPSharedMassCase : public Test {
public:
	InferenceTestCompareVsIPFPSharedMassCase();
	void runTest();
};

class InferenceTestSpectrumMessage : public Test {
public:
	InferenceTestSpectrumMessage();
	void runTest();
};

class InferenceTestSpectrumMessageNoisePeak : public Test {
public:
	InferenceTestSpectrumMessageNoisePeak();
	void runTest();
};

class InferenceTestSpectrumMessageWithIsotopes : public Test {
public:
	InferenceTestSpectrumMessageWithIsotopes();
	void runTest();
};

class InferenceTestSpectrumMessageWithIsotopesAndNoisePeak : public Test {
public:
	InferenceTestSpectrumMessageWithIsotopesAndNoisePeak();
	void runTest();
};

#endif // __INFER_TEST_H__