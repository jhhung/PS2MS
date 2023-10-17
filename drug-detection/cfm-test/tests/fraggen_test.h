/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# fraggen_test.h
#
# Description: 	Test code for Fragment Tree Generation
#
# Author: Felicity Allen
# Created: August 2014
#########################################################################*/

#ifndef __FRAGGEN_TEST_H__
#define __FRAGGEN_TEST_H__

#include "../test.h"
#include "FragmentGraphGenerator.h"

class FragGenTestPositiveESI : public Test {
public:
	FragGenTestPositiveESI();
	void runTest();
};

class FragGenTestRingPositiveESI : public Test {
public:
	FragGenTestRingPositiveESI();
	void runTest();
};

class FragGenTestPositiveESISplitCharge : public Test {
public:
	FragGenTestPositiveESISplitCharge();
	void runTest();
};

class FragGenTestNegativeESI : public Test {
public:
	FragGenTestNegativeESI();
	void runTest();
};

class FragGenTestRingNegativeESI : public Test {
public:
	FragGenTestRingNegativeESI();
	void runTest();
};

class FragGenTestPositiveEI : public Test {
public:
	FragGenTestPositiveEI();
	void runTest();
};

class FragGenTestPositiveEIMultibreak : public Test {
public:
	FragGenTestPositiveEIMultibreak();
	void runTest();
};

class FragGenTestPositiveEITriple : public Test {
public:
	FragGenTestPositiveEITriple();
	void runTest();
};

class FragGenTestPositiveEIOxygenAromatic : public Test {
public:
	FragGenTestPositiveEIOxygenAromatic();
	void runTest();
};

class FragGenTestPositiveEIAlkane : public Test {
public:
	FragGenTestPositiveEIAlkane();
	void runTest();
};

class FragGenTestRingPositiveEI : public Test {
public:
	FragGenTestRingPositiveEI();
	void runTest();
};

class FragGenTestPositiveEISplitCharge : public Test {
public:
	FragGenTestPositiveEISplitCharge();
	void runTest();
};

class FragGenTestPositiveEINistExceptions : public Test {
public:
	FragGenTestPositiveEINistExceptions();
	void runTest();
};

class FragGenTestPositiveEIAndESIDegreeLpBonding : public Test {
public:
	FragGenTestPositiveEIAndESIDegreeLpBonding();
	void runTest();
};

class FragGenTestCasesFromGross : public Test {
public:
	FragGenTestCasesFromGross();
	void runTest();
};

class FragGenTestMaxElectronMovement : public Test {
public:
	FragGenTestMaxElectronMovement();
	void runTest();
};

class FragGenTestDisallowDetourTransitions : public Test {
public:
	FragGenTestDisallowDetourTransitions();
	void runTest();
};

class FragGenTestMaxRingBreaks : public Test {
public:
	FragGenTestMaxRingBreaks();
	void runTest();
};

class FVFragGraphSaveAndLoadState : public Test {
public:
	FVFragGraphSaveAndLoadState();
	void runTest();
};

#endif // __FRAGGEN_TEST_H__