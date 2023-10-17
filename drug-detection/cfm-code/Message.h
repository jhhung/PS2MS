/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# message.h
#
# Description: 	Utility class for storing and passing around sparse
#               messages.
#				
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __MESSAGE_H__
#define __MESSAGE_H__

#include <boost/numeric/ublas/vector_sparse.hpp>
#include <limits>

static const double DBL_MAXIMUM = std::numeric_limits<double>::max();

//Log Domain Double (converts null entry to -DBL_MAX)
class logdbl {
  double value;
public:
  logdbl(double f=-DBL_MAXIMUM)
  {
    value = f;
  }

  logdbl& operator=(double f)
  {
    value = f;
    return *this;
  }

  operator double()   { return value; }
  operator double() const { return value; }
};

//Normalized sparse vector in the log domain
class Message{

public:
	Message(unsigned int size = 0);
	double getIdx(unsigned int i);
	void addToIdx(unsigned int i, double value);
	void multIdx(unsigned int i, double value);
	void reset(unsigned int size);
	unsigned int isEmpty();
	Message &operator=(const Message &rhs);
	void addWeightedMessage( Message &msg, double weight );
	typedef boost::numeric::ublas::mapped_vector<logdbl>::iterator iterator;
	typedef boost::numeric::ublas::mapped_vector<logdbl>::const_iterator const_iterator;	
	iterator begin(){ if( !normalized ) normalize(); return message.begin(); }
	iterator end(){ return message.end(); }
	unsigned int size(){ return message.size(); }
	void print();

private:
	double log_sum;
	boost::numeric::ublas::mapped_vector<logdbl> message;
	unsigned int empty;
	unsigned int normalized;
	void normalize();
};


#endif // __MESSAGE_H__
