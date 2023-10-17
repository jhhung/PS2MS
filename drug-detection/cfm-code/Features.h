/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# Features.h
#
# Description: 	Code for computing features for fragmentations.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#ifndef __FEATURE_H__
#define __FEATURE_H__

#include "Util.h"
#include "FunctionalGroups.h"

#include <GraphMol/FragCatalog/FragCatParams.h>

#include <boost/ptr_container/ptr_vector.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

struct input_file_t;

//Exception to throw when the input feature configuration file is invalid 
class InvalidConfigException: public std::exception{

	virtual const char* what() const throw(){
		return "Invalid Feature Configuration File";
	}
};

class FeatureCalculationException: public std::exception{
private:
    std::string message_;
public:
	FeatureCalculationException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << "Error computing feature vector: " << message_ << std::endl;
		return message_.c_str();
	}
	~FeatureCalculationException() throw() {};
};

class FeatureHelperException: public std::exception{
private:
    std::string message_;
public:
	FeatureHelperException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << "Error in FeatureHelper: " << message_ << std::endl;
		return message_.c_str();
	}
	~FeatureHelperException() throw() {};
};

//Structure to hold a sparse computed feature vector
typedef unsigned int feature_t;
class FeatureVector{
public:
	FeatureVector(){fv_idx = 0;};
	void addFeature( double value );
	void addFeatureAtIdx( double value, unsigned int idx );
	unsigned int const getTotalLength() const {return fv_idx;};
	feature_t const getFeature( int idx ) const { return fv[idx];};
	std::vector<feature_t>::const_iterator getFeatureBegin() const { return fv.begin(); };
	std::vector<feature_t>::const_iterator getFeatureEnd() const { return fv.end(); };
	unsigned int const getNumSetFeatures() const { return fv.size();};

private:
	std::vector<feature_t> fv;
	unsigned int fv_idx;
};

//Base class to compute a feature - all features should inherit from this
class Feature{

public:
	virtual void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const = 0;
	unsigned int getSize() const { return size; };
	std::string getName() const { return name; };
	virtual ~Feature(){};

protected:
	unsigned int size;	
	std::string name;
	static const std::vector<std::string> &OKsymbols();
	static const std::vector<std::string> &OKSymbolsLess();
	void replaceUncommonWithX( std::string &symbol ) const;

};

//Class to compute a feature vector
class FeatureCalculator{

public:
	//Constructor: Initialise the calculator using a config file listing features
	FeatureCalculator( std::string &config_filename );

	//Constructor: Initialise the calculator using a list of feature names
	FeatureCalculator( std::vector<std::string> &feature_list );
	
	//Compute the expected number of total features
	unsigned int getNumFeatures();
	
	//Retrieve the list of feature names being used
	std::vector<std::string> getFeatureNames();

	//Retrieve a list of valid feature names (for testing)
	static const std::vector<std::string> getValidFeatureNames();

	//Compute the feature vector for the input ion and nl (with labeled Root atoms)
	// - NB: responsibility of caller to delete.
	FeatureVector *computeFV( const RootedROMolPtr *ion, const RootedROMolPtr *nl );

	bool includesFeature( const std::string &fname );

private:
	//List of feature classes ready to be used
	static const boost::ptr_vector<Feature> &featureCogs();

	//Indexes of feature classes that are selected for use
	std::vector<int> used_feature_idxs;

	//Helper function - Configure feature for use
	void configureFeature( std::string &name );
};

//*************************
//FEATURE IMPLEMENTATIONS:
//*************************

class BreakAtomPair : public Feature {
public:
	BreakAtomPair(){ size = 72; name = "BreakAtomPair"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class RootPathFeature : public Feature {
protected:
	typedef std::vector<std::string> path_t;
	void computeRootPaths(std::vector<path_t> &paths, const RootedROMolPtr *mol, int len, bool ring_break) const;
	void addRootPairFeatures(FeatureVector &fv, std::vector<path_t> &paths, int ring_break) const;
	void addRootTripleFeatures(FeatureVector &fv, std::vector<path_t> &paths, int ring_break) const;
private:
	void addPathsFromAtom( std::vector<path_t> &paths, const RDKit::Atom *atom, const romol_ptr_t mol, const RDKit::Atom *prev_atom, path_t &path_so_far, int len ) const;
};

class IonRootPairs : public RootPathFeature {
public:
	IonRootPairs(){ size = 145; name = "IonRootPairs"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootPairs : public RootPathFeature {
public:
	NLRootPairs(){ size = 145; name = "NLRootPairs"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonRootTriples : public RootPathFeature {
public:
	IonRootTriples(){ size = 865; name = "IonRootTriples"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootTriples : public RootPathFeature {
public:
	NLRootTriples(){ size = 865; name = "NLRootTriples"; }; 
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class GasteigerCharges : public RootPathFeature {
public:
	GasteigerCharges(){ size = 72; name = "GasteigerCharges"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl)  const;
private:
	int discretizeGasteigerCharge( double gc ) const;
};

class HydrogenMovement : public Feature {
public:
	HydrogenMovement(){ size = 10; name = "HydrogenMovement"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class HydrogenRemoval : public Feature {
public:
	HydrogenRemoval(){ size = 10; name = "HydrogenRemoval"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class BrokenOrigBondType : public Feature {
public:
	BrokenOrigBondType(){ size = 7; name = "BrokenOrigBondType"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NeighbourOrigBondTypes : public Feature {
public:
	NeighbourOrigBondTypes(){ size = 12; name = "NeighbourOrigBondTypes"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class RadicalFeatures : public Feature {
public:
	RadicalFeatures(){ size = 3; name = "RadicalFeatures"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonicFeatures: public Feature {
public:
	IonicFeatures(){ size = 5; name = "IonicFeatures"; };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class RingFeatures : public Feature {
public:
	RingFeatures(){ size = 12; name = "RingFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
private:
	//Helper function - compute the distance between two root 
	//atoms in a molecule (assumes ring break)
	int calcRootDistance(const RootedROMolPtr *mol)  const;
};

class ExtraRingFeatures : public Feature {
public:
	ExtraRingFeatures(){ size = 3; name = "ExtraRingFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class IonRootMMFFAtomType : public Feature {
public:
	IonRootMMFFAtomType(){ size = 100; name = "IonRootMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootMMFFAtomType : public Feature {
public:
	NLRootMMFFAtomType(){ size = 100; name = "NLRootMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class RootAtomFeature : public Feature {
protected:
	void computeRootAtomFeature( FeatureVector &fv, const RootedROMolPtr *mol, bool ring_break ) const;
};

class IonRootAtom : public RootAtomFeature {
public:
	IonRootAtom(){ size = 13; name = "IonRootAtom";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLRootAtom : public RootAtomFeature {
public:
	NLRootAtom(){ size = 13; name = "NLRootAtom";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};


class NeighbourMMFFFeature : public Feature {
protected:
	void addNeighbourAtomTypes( FeatureVector &fv, const RootedROMolPtr *mol, const RDKit::Atom *root, int offset ) const;

};

class IonNeighbourMMFFAtomType : public NeighbourMMFFFeature {
public:
	IonNeighbourMMFFAtomType(){ size = 101; name = "IonNeighbourMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class NLNeighbourMMFFAtomType : public NeighbourMMFFFeature {
public:
	NLNeighbourMMFFAtomType(){ size = 101; name = "NLNeighbourMMFFAtomType";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl )  const;
};

class FunctionalGroupFeature : public Feature {
protected:
	void addFunctionalGroupFeaturesFromAtom( std::vector<int> &tmp_full_fv, const RDKit::Atom *atom, const romol_ptr_t mol, const RDKit::Atom *prev_atom, int max_depth, int depth, bool extra ) const;
	void addFunctionalGroupFeatures( FeatureVector &fv, const RootedROMolPtr *mol, int max_dist, int is_ring_break, bool extra = false) const;
};

class IonFunctionalGroupFeatures : public FunctionalGroupFeature {
public:
	IonFunctionalGroupFeatures(){ size = (NUM_FGRPS+1)*2; name = "IonFunctionalGroupFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class NLFunctionalGroupFeaturesD2 : public FunctionalGroupFeature {
public:
	NLFunctionalGroupFeaturesD2(){ size = (NUM_FGRPS+1)*3; name = "NLFunctionalGroupFeaturesD2";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class IonFunctionalGroupFeaturesD2 : public FunctionalGroupFeature {
public:
	IonFunctionalGroupFeaturesD2(){ size = (NUM_FGRPS+1)*3; name = "IonFunctionalGroupFeaturesD2";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class NLFunctionalGroupFeatures : public FunctionalGroupFeature {
public:
	NLFunctionalGroupFeatures(){ size = (NUM_FGRPS+1)*2; name = "NLFunctionalGroupFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class IonFunctionalGroupRootOnlyFeatures : public FunctionalGroupFeature {
public:
	IonFunctionalGroupRootOnlyFeatures(){ size = NUM_FGRPS+1; name = "IonFunctionalGroupRootOnlyFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class NLFunctionalGroupRootOnlyFeatures : public FunctionalGroupFeature {
public:
	NLFunctionalGroupRootOnlyFeatures(){ size = NUM_FGRPS+1; name = "NLFunctionalGroupRootOnlyFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class IonExtraFunctionalGroupFeatures : public FunctionalGroupFeature {
public:
	IonExtraFunctionalGroupFeatures(){ size = (NUM_EXTRA_FGRPS+1)*2; name = "IonExtraFunctionalGroupFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class NLExtraFunctionalGroupFeatures : public FunctionalGroupFeature {
public:
	NLExtraFunctionalGroupFeatures(){ size = (NUM_EXTRA_FGRPS+1)*2; name = "NLExtraFunctionalGroupFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};


class QuadraticFeatures : public Feature {
public:
	QuadraticFeatures(){ size = 0; name = "QuadraticFeatures";  };
	void compute( FeatureVector &fv, const RootedROMolPtr *ion, const RootedROMolPtr *nl ) const;
};

class FeatureHelper{
public:
	FeatureHelper(){ exec_flags.resize(6); for(int i=0;i<6;i++) exec_flags[i] = 0; };
	FeatureHelper( FeatureCalculator *fc ){
		exec_flags.resize(6);
		exec_flags[0] = fc->includesFeature("GasteigerCharges");
		exec_flags[1] = fc->includesFeature("HydrogenMovement") || fc->includesFeature("HydrogenRemoval");
		exec_flags[2] = fc->includesFeature("IonRootMMFFAtomType") || fc->includesFeature("NLRootMMFFAtomType") ||
			fc->includesFeature("IonNeighbourMMFFAtomType") || fc->includesFeature("NLNeighbourMMFFAtomType");
		exec_flags[3] = fc->includesFeature("BrokenOrigBondType") || fc->includesFeature("NeighbourOrigBondTypes");
		exec_flags[4] = fc->includesFeature("IonFunctionalGroupFeatures") || fc->includesFeature("NLFunctionalGroupFeatures") ||
						fc->includesFeature("IonFunctionalGroupFeaturesD2") || fc->includesFeature("NLFunctionalGroupFeaturesD2") ||
						fc->includesFeature("IonFunctionalGroupRootOnlyFeatures") || fc->includesFeature("NLFunctionalGroupRootOnlyFeatures");
		exec_flags[5] = fc->includesFeature("IonExtraFunctionalGroupFeatures") || fc->includesFeature("NLExtraFunctionalGroupFeatures");
		if( exec_flags[4]){
			//fparams = new RDKit::FragCatParams( 0, 20, "cfmid_functional_groups.csv" );	
			//std::string fgrps_serial = fparams->Serialize();
			//std::ofstream of;
			//of.open("functional_groups_serial.txt");
			//of << fgrps_serial << std::endl;
			//of.close();
			fparams = new RDKit::FragCatParams( FGRPS_PICKLE );
			if( fparams->getNumFuncGroups() != NUM_FGRPS ) 
				throw FeatureHelperException("Mismatch in expected and found number of functional groups");
		}
		if( exec_flags[5] ){
			xfparams = new RDKit::FragCatParams( EXTRA_FGRPS_PICKLE );
			if( xfparams->getNumFuncGroups() != NUM_EXTRA_FGRPS ) 
				throw FeatureHelperException("Mismatch in expected and found number of extra functional groups");	
		
		}
	};
	~FeatureHelper(){
		if( exec_flags[4]) delete fparams;
		if( exec_flags[5]) delete xfparams;
	}

	void addLabels( RDKit::RWMol *rwmol){
		initialiseRoots(rwmol);
		labelAromatics(rwmol);
		if( exec_flags[0] ) labelGasteigers(rwmol); 
		if( exec_flags[1] ) labelOriginalMasses(rwmol); 
		if( exec_flags[2] ) labelMMFFAtomTypes(rwmol); 
		if( exec_flags[3] ) labelOriginalBondTypes(rwmol);
		if( exec_flags[4] ) labelFunctionalGroups(rwmol, false);
		if( exec_flags[5] ) labelFunctionalGroups(rwmol, true);
		labelAtomsWithLonePairs(rwmol);
	};
	bool getExecFlag( unsigned int idx ){ return exec_flags[idx]; };

private:

	std::vector<int> exec_flags;
	//Helper functions - used to create labels on atoms and bonds, 
	//that will be used in Feature Calculations and can't be computed once
	//a molecule is broken
	static void initialiseRoots( RDKit::RWMol *rwmol );
	static void labelGasteigers( RDKit::RWMol *rwmol );
	static void labelAromatics( RDKit::RWMol *rwmol );
	static void labelOriginalMasses( RDKit::RWMol *rwmol );
	static void labelMMFFAtomTypes( RDKit::RWMol *rwmol );
	static void labelAtomsWithLonePairs( RDKit::RWMol *rwmol );
	static void labelOriginalBondTypes( RDKit::RWMol *rwmol );
	void labelFunctionalGroups( RDKit::RWMol *rwmol, bool extra );	//Not static because it uses fparams.

	RDKit::FragCatParams *fparams;
	RDKit::FragCatParams *xfparams;
};


#endif // __FEATURE_H__