/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# comms_test.h
#
# Description: 	Test code for Comms.cpp
#
# Author: Felicity Allen
# Created: March 2013
#########################################################################*/

#ifndef __COMMS_TEST_H__
#define __COMMS_TEST_H__

#include "../test.h"
#include "Comms.h"

class CommsTestSetMasterUsedIdxs : public Test {
public:
	CommsTestSetMasterUsedIdxs();
	void runTest();
};

class CommsTestCollectQInMaster : public Test {
public:
	CommsTestCollectQInMaster();
	void runTest();
};

class CommsTestCollectGradsInMaster : public Test {
public:
	CommsTestCollectGradsInMaster();
	void runTest();
};

class CommsTestBroadcastParams : public Test {
public:
	CommsTestBroadcastParams();
	void runTest();
};

#endif // __COMMS_TEST_H__