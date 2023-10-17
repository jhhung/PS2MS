/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# message_test.h
#
# Description: 	Test code for Message.cpp
#
# Author: Felicity Allen
# Created: June 2013
#########################################################################*/

#ifndef __MESSAGE_TEST_H__
#define __MESSAGE_TEST_H__

#include "../test.h"
#include "Message.h"

class MessageTestBasics : public Test {
public:
	MessageTestBasics();
	void runTest();
};

class MessageTestIteration : public Test {
public:
	MessageTestIteration();
	void runTest();
};

class MessageTestCopyAssignment : public Test {
public:
	MessageTestCopyAssignment();
	void runTest();
};

#endif // __MESSAGE_TEST_H__