/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# FunctionalGroups.h
#
# Description: 	Pickle of the RDKit::FragCatParams object containing 
#				functional groups to be used in the FunctionalGroups features.
#
# Copyright (c) 2015, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FGRPS_H__
#define __FGRPS_H__

#include <string>

static const int NUM_FGRPS = 161;

extern const std::string FGRPS_PICKLE;

static const int NUM_EXTRA_FGRPS = 13;

extern const std::string EXTRA_FGRPS_PICKLE;

#endif // __FGRPS_H__