/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# test.h
#
# Description: 	Test code
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __TEST_H__
#define __TEST_H__

#include <string>

class Test{

public:
	Test();
	void run();
	virtual void runTest();
	bool passed;
	bool excpt_occurred;
	std::string description;
	virtual ~Test(){}
};

#endif // __TEST_H__