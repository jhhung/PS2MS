/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# message.cpp
#
# Description: 	Utility class for storing and passing around sparse
#               messages.
#				
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Message.h"
#include "Util.h"

Message::Message(unsigned int size){
	reset(size);
}


double Message::getIdx(unsigned int i){
	
	if( !normalized ) normalize();
	logdbl dbl = message[i];
	return double(dbl);
}

void Message::addWeightedMessage( Message &msg, double weight ){
	Message::const_iterator it = msg.begin();
	for( ; it != msg.end(); it++ )
		addToIdx( it.index(), *it + log(weight) );
}

void Message::addToIdx(unsigned int i, double value){
	logdbl dbl = message[i];
	message[i] = logAdd( dbl, value );
	empty = 0;
	normalized = 0;
}

void Message::multIdx(unsigned int i, double value){
	logdbl dbl = message[i];
	if( dbl > -DBL_MAXIMUM ) message[i] = dbl + value;
	else message[i] = value;
	empty = 0;
	normalized = 0;
}

void Message::normalize(){
	log_sum = -DBL_MAXIMUM;
	iterator it = message.begin();
	for( ; it != message.end(); ++it )
		log_sum = logAdd( log_sum, *it );
	it = message.begin();
	for( ; it != message.end(); ++it )
		*it = *it - log_sum;
	log_sum = 0;
	normalized = 1;
}

void Message::reset(unsigned int size){
	empty = 1;
	normalized = 1;
	message.resize( size );
	message.clear();
	log_sum = -DBL_MAXIMUM;
}

unsigned int Message::isEmpty(){
	return empty;
}

Message &Message::operator=(const Message &rhs){

	log_sum = rhs.log_sum;
	empty = rhs.empty;
	message = rhs.message;
	normalized = rhs.normalized;
	return *this;
}

void Message::print(){
	if( !normalized ) normalize();
	iterator it = message.begin();
	for( ; it != message.end(); ++it )
		if( exp(*it) > 0.001 ) 
			std::cout << it.index() << ": " << exp(*it) << ", ";
	std::cout << std::endl;
}
