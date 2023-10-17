/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# comparator_test.h
#
# Description: 	Test code for Comparators.cpp
#
# Author: Felicity Allen
# Created: January 2013
#########################################################################*/

#ifndef __COMPARATOR_TEST_H__
#define __COMPARATOR_TEST_H__

#include "../test.h"
#include "Comparators.h"

class ComparatorsTestRecall : public Test {
public:
	ComparatorsTestRecall();
	void runTest();
};

class ComparatorsTestWeightedRecall : public Test {
public:
	ComparatorsTestWeightedRecall();
	void runTest();
};


class ComparatorsTestPrecision : public Test {
public:
	ComparatorsTestPrecision();
	void runTest();
};

class ComparatorsTestWeightedPrecision : public Test {
public:
	ComparatorsTestWeightedPrecision();
	void runTest();
};

class ComparatorsTestJaccard : public Test {
public:
	ComparatorsTestJaccard();
	void runTest();
};


#endif // __COMPARATOR_TEST_H__