/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# features_test.h
#
# Description: 	Test code for Features.cpp
#
# Author: Felicity Allen
# Created: November 2012
#########################################################################*/

#ifndef __FEATURE_TEST_H__
#define __FEATURE_TEST_H__

#include "../test.h"
#include "Features.h"

class FeaturesTestInit : public Test {
public:
	FeaturesTestInit();
	void runTest();
};

class FeaturesTestBreakAtomPair : public Test {
public:
	FeaturesTestBreakAtomPair();
	void runTest();
};

class FeaturesTestRootPairs : public Test {
public:
	FeaturesTestRootPairs();
	void runTest();
};

class FeaturesTestRootTriples : public Test {
public:
	FeaturesTestRootTriples();
	void runTest();
};

class FeaturesTestGasteigerCharges : public Test {
public:
	FeaturesTestGasteigerCharges();
	void runTest();
};

class FeaturesTestHydrogenMovement : public Test {
public:
	FeaturesTestHydrogenMovement();
	void runTest();
};

class FeaturesTestRadicalFeatures : public Test {
public:
	FeaturesTestRadicalFeatures();
	void runTest();
};

class FeaturesTestIonicFeatures : public Test {
public:
	FeaturesTestIonicFeatures();
	void runTest();
};

class FeaturesTestRingFeatures : public Test {
public:
	FeaturesTestRingFeatures();
	void runTest();
};

class FeaturesTestExtraRingFeatures : public Test {
public:
	FeaturesTestExtraRingFeatures();
	void runTest();
};

class FeaturesTestFunctionalGroups : public Test {
public:
	FeaturesTestFunctionalGroups();
	void runTest();
};

class FeaturesTestExtraFunctionalGroups : public Test {
public:
	FeaturesTestExtraFunctionalGroups();
	void runTest();
};

class FeaturesTestFunctionalGroupsRootOnly : public Test {
public:
	FeaturesTestFunctionalGroupsRootOnly();
	void runTest();
};

class FeaturesTestRootMMFFAtomType :  public Test {
public:
	FeaturesTestRootMMFFAtomType();
	void runTest();
};

class FeaturesTestNeighbourMMFFAtomType :  public Test {
public:
	FeaturesTestNeighbourMMFFAtomType();
	void runTest();
};

class FeaturesTestBrokenOrigBondType :  public Test {
public:
	FeaturesTestBrokenOrigBondType();
	void runTest();
};

class FeaturesTestNeighbourOrigBondType : public Test {
public:
	FeaturesTestNeighbourOrigBondType();
	void runTest();
};

class FeaturesTestRootAtom :  public Test {
public:
	FeaturesTestRootAtom();
	void runTest();
};

class FeaturesTestQuadraticFeatures : public Test {
public:
	FeaturesTestQuadraticFeatures();
	void runTest();
};

class FeaturesTestLength : public Test {
public:
	FeaturesTestLength();
	void runTest();
};

class FeaturesTestMetlinExample : public Test {
public:
	FeaturesTestMetlinExample();
	void runTest();
};

#endif // __FEATURE_TEST_H__