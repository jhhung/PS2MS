/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# message_test.cpp
#
# Description: Test code for message.cpp
#
# Author: Felicity Allen
# Created: June 2013
#########################################################################*/

#include "message_test.h"

MessageTestBasics::MessageTestBasics(){
	description = "Test message basics";
}

void MessageTestBasics::runTest(){
	
	bool pass = true;
	double tol = 1e-8;
	
	//Create a message
	Message m(150);
	m.addToIdx(100, 0.55);
	m.addToIdx(123, 1.65);
	m.addToIdx(12, -0.8);

	//Get the values back and check they are correctly normalized
	double tmp = m.getIdx(100);
	if( fabs( tmp - -1.4500684380957891 ) > tol ){
		std::cout << "Expecting -1.4500684380957891 at idx 100, but found " << tmp << std::endl;
		pass = false;
	}
	tmp = m.getIdx(123);
	if( fabs( tmp - -0.35006843809578925 ) > tol ){
		std::cout << "Expecting -0.35006843809578925 at idx 123, but found " << tmp << std::endl;
		pass = false;
	}
	tmp = m.getIdx(12);
	if( fabs( tmp - -2.8000684380957894 ) > tol ){
		std::cout << "Expecting -2.8000684380957894 at idx 12, but found " << tmp << std::endl;
		pass = false;
	}
	tmp = m.getIdx(5);	//Other index
	if( tmp > -DBL_MAXIMUM ){
		std::cout << "Expecting -DBL_MAXIMUM at idx 5, but found " << tmp << std::endl;
		pass = false;
	}

	passed = pass;
}

MessageTestIteration::MessageTestIteration(){
	description = "Test message iteration";
}

void MessageTestIteration::runTest(){
	
	bool pass = true;
	double tol = 1e-8;

	//Create that same message again
	Message m(150);
	m.addToIdx(100, 0.55);
	m.addToIdx(123, 1.65);
	m.addToIdx(12, -0.8);

	//Iterate through the components
	Message::const_iterator it = m.begin();
	if( it.index() != 12 || fabs( *it - -2.8000684380957894 ) > tol ){
		std::cout << "Expecting -2.8000684380957894 at idx 12, but found " << *it << " at " << it.index() << std::endl;
		pass = false;
	}
	++it;
	if( it.index() != 100 || fabs( *it - -1.4500684380957891 ) > tol ){
		std::cout << "Expecting -1.4500684380957891 at idx 100, but found " << *it << " at " << it.index() << std::endl;
		pass = false;
	}
	++it;
	if( it.index() != 123 || fabs( *it - -0.35006843809578925 ) > tol ){
		std::cout << "Expecting -0.35006843809578925 at idx 123, but found " << *it << " at " << it.index() << std::endl;
		pass = false;
	}

	passed = pass;

}

MessageTestCopyAssignment::MessageTestCopyAssignment(){
	description = "Test message copy assignment";
}

void MessageTestCopyAssignment::runTest(){
	
	bool pass = true;
	double tol = 1e-8;

	//Create that same message again
	Message m(150);
	m.addToIdx(100, 0.55);
	m.addToIdx(123, 1.65);
	m.addToIdx(12, -0.8);

	//Copy it
	Message m2 = m;

	//Change the old message
	m.addToIdx(5, 0.25);

	//Check index 5
	double tmp = m.getIdx(5);
	double tmp2 = m2.getIdx(5);
	if( fabs( tmp - tmp2 ) < tol || tmp2 > -DBL_MAXIMUM || fabs( tmp - -1.910282456760932) > tol){
		std::cout << "Unexpected values at index 5 after copy: " << tmp << " vs " << tmp2 << std::endl;
		pass = false;
	}

	//Check index 123
	tmp = m.getIdx(123);
	tmp2 = m2.getIdx(123);
	if( fabs( tmp - -0.510282456760932) > tol || fabs( tmp2 - -0.350068438095789) > tol ){
		std::cout << "Unexpected values at index 123 after copy: " << tmp << " vs " << tmp2 << std::endl;
		pass = false;
	}

	passed = pass;

}

//MessageTestMultiply::MessageTestMultiply(){
//	description = "Test message multiplication";
//}
//
//void MessageTestMultiply::runTest(){
//	
//	bool pass = true;
//	double tol = 1e-8;
//	
//
//	passed = pass;
//}