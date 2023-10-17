/*#########################################################################
# Mass Spec Prediction of HMDB Molecules
#
# inference_test.h
#
# Description: 	Test code for inference.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#ifndef __LIBDAI_INFERENCE_TEST_H__
#define __LIBDAI_INFERENCE_TEST_H__

#include "../test.h"

class InferenceTestvsIPFPDepth2 : public Test {
public:
	InferenceTestvsIPFPDepth2();
	void runTest();
};


class InferenceTestSimpleDepth3 : public Test {
public:
	InferenceTestSimpleDepth3();
	void runTest();
};

class InferenceTestComplexDepth3 : public Test {
public:
	InferenceTestComplexDepth3();
	void runTest();
};

class InferenceTestSimpleDepth2 : public Test {
public:
	InferenceTestSimpleDepth2();
	void runTest();
};

class InferenceTestComplexDepth2 : public Test {
public:
	InferenceTestComplexDepth2();
	void runTest();
};

class InferenceTestSimpleDepth1 : public Test {
public:
	InferenceTestSimpleDepth1();
	void runTest();
};

class InferenceTestComplexDepth1 : public Test {
public:
	InferenceTestComplexDepth1();
	void runTest();
};

class InferenceTestComputeTransitionProbs : public Test {
public:
	InferenceTestComputeTransitionProbs();
	void runTest();
};

#endif // __LIBDAI_INFERENCE_TEST_H__