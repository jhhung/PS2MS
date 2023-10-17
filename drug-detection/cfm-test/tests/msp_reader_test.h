/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# msp_reader_test.h
#
# Description: 	Test code for MspReader
#
# Author: Felicity Allen
# Created: August 2014
#########################################################################*/

#ifndef __MSP_READER_TEST_H__
#define __MSP_READER_TEST_H__

#include "../test.h"

class MspReaderTest : public Test {
public:
	MspReaderTest();
	void runTest();
};

class MspReaderMultipleEnergiesTest : public Test {
public:
	MspReaderMultipleEnergiesTest();
	void runTest();
};

#endif // __MSP_READER_TEST_H__