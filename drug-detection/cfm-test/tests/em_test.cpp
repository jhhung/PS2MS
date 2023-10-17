/*#########################################################################
#  Mass Spec Prediction and Identification of Metabolites
#
# em_test.cpp
#
# Description: Test code for EM.cpp
#
# Author: Felicity Allen
# Created: August 2013
#########################################################################*/

#include "mpi.h"
#include "time.h"

#include "em_test.h"

#include "EM.h"
#include "EM_NN.h"

bool compareSpectra(const Spectrum *orig_spec, const Spectrum *predicted_spec, double abs_mass_tol, double ppm_mass_tol, double intensity_tol  ){
	
	bool pass = true;
	Spectrum::const_iterator ito = orig_spec->begin();
	for( ; ito != orig_spec->end(); ++ito ){
			
		//Find a peak in the predicted spectrum with the same mass
		Spectrum::const_iterator itp = predicted_spec->begin();
		int found = 0;
		for( ; itp != predicted_spec->end(); ++itp ){
			double mass_tol = getMassTol( abs_mass_tol, ppm_mass_tol, itp->mass );
			if( fabs(itp->mass - ito->mass ) < mass_tol ){
					
				//Check the intensity values of the matching peaks
				if( fabs(itp->intensity - ito->intensity) > intensity_tol ){
					std::cout << "Mismatch in predicted peak intensity for mass " << ito->mass << ": Expecting " << ito->intensity << " but found " << itp->intensity << std::endl;
					pass = false;
				}
				found = 1; 
				break;
			}
		}
		if( !found ){
			std::cout << "Could not find matching predicted peak for mass " << ito->mass << std::endl;
			pass = false;
		}
	}
	return pass;
}

EMTestSelfProduction::EMTestSelfProduction(){
	description = "Test integrated CE-CFM EM ability to learn single spectrum";
}

void EMTestSelfProduction::runTest(){
	
	bool pass = true; 
	double intensity_tol = 2.0; //For intensity in range 0-100
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( cfg, param_cfg_file );
	cfg.lambda = 0.0000001;
	cfg.ga_converge_thresh = 0.00001;
	cfg.include_h_losses = true; 

	for( int config_state = 0; config_state <= 3; config_state++ ){

		if( config_state == 1 ){
			std::cout << "Testing Depth 1-2-3 Configuration" << std::endl;
			cfg.spectrum_depths[0] = 1; cfg.spectrum_depths[1] = 2;  cfg.spectrum_depths[2] = 3;
			cfg.model_depth = 3;
			initDerivedConfig(cfg);
		}else std::cout << "Testing Depth 2-4-6 Configuration" << std::endl;

		//Feature Calculator
		std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
		FeatureCalculator fc( feature_cfg_file );

		//Prepare some simple data
		std::vector<MolData> data;
		std::string id = "TestMol", smiles = "NCCCN";
		std::string spec_file = "tests/test_data/example_spectra.txt";
		data.push_back( MolData( id, smiles, 0, &cfg ) );
		data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
		data[0].readInSpectraFromFile( spec_file );

		//Run EM
		double best_Q = -100000.0;
		std::string param_filename = "tmp_param_output.log";
		for( int trial = 0; trial < 3; trial++ ){
			std::string status_file = "tmp_status_file.log";
			std::string tmp_file = "tmp.log";
			EM em( &cfg, &fc, status_file);
			double Q = em.run( data, 1, tmp_file  );
			if( Q > best_Q ){
				em.writeParamsToFile( param_filename );
				best_Q = Q;
			}
		}

		//Predict the output spectra
		Param param( param_filename );
		data[0].computePredictedSpectra( param );
		data[0].postprocessPredictedSpectra(100.0, 0, 100000);

		//Compare the original and predicted spectra - should be able to overfit
		//very close to the actual values since training on same (and only same) mol
		for( unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++ ){
			const Spectrum *orig_spec = data[0].getSpectrum(energy);
			const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
			pass &= compareSpectra(orig_spec, predicted_spec, cfg.abs_mass_tol, cfg.ppm_mass_tol, intensity_tol);
		}
	}

	passed = pass;

}

EMTestSingleEnergySelfProduction::EMTestSingleEnergySelfProduction(){
	description = "Test integrated SE-CFM EM ability to learn single spectrum";
}

void EMTestSingleEnergySelfProduction::runTest(){
	
	bool pass = true; 
	double intensity_tol = 2.0; //For intensity in range 0-100
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t orig_cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( orig_cfg, param_cfg_file );
	orig_cfg.lambda = 0.0000001;
	orig_cfg.use_single_energy_cfm = 1;
	orig_cfg.spectrum_depths[1] = 2;
	orig_cfg.spectrum_depths[2] = 2;
	orig_cfg.include_h_losses = true; 

	//Feature Calculator
	std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( feature_cfg_file );

	//Prepare some simple data
	std::vector<MolData> data;
	std::string id = "TestMol", smiles = "NCCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";
	data.push_back( MolData( id, smiles, 0, &orig_cfg ) );
	data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
	data[0].readInSpectraFromFile( spec_file );

	Param *final_params;
	for( int energy = 0; energy < 3; energy++ ){

		config_t cfg;
		initSingleEnergyConfig( cfg, orig_cfg, energy );

		//Run EM
		std::string status_file = "tmp_status_file.log";
		std::string tmp_file = "tmp.log";
		EM em( &cfg, &fc, status_file );
		double Q = em.run( data, 1, tmp_file );
		std::string param_filename = "tmp_param_output.log";
		em.writeParamsToFile( param_filename );

		if(energy == 0 ) final_params = new Param(param_filename);
		else{ 
			Param eparam(param_filename);
			final_params->appendNextEnergyParams( eparam, energy );
		}
	}

	//Predict the output spectra
	data[0].computePredictedSpectra( *final_params );
	data[0].postprocessPredictedSpectra( 100.0, 0, 1000 );

	//Compare the original and predicted spectra - should be able to overfit
	//very close to the actual values since training on same (and only same) mol
	for( unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++ ){
		const Spectrum *orig_spec = data[0].getSpectrum(energy);
		const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
		pass &= compareSpectra(orig_spec, predicted_spec, orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, intensity_tol);
		
	}
	passed = pass;
}


EMTestNNSingleEnergySelfProduction::EMTestNNSingleEnergySelfProduction(){
	description = "Test Neural Net SE-CFM EM ability to learn single spectrum from itself";
}

void EMTestNNSingleEnergySelfProduction::runTest(){
	
	bool pass = true; 
	double intensity_tol = 2.0; //For intensity in range 0-100
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t orig_cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( orig_cfg, param_cfg_file );
	orig_cfg.lambda = 0.0000001;
	orig_cfg.use_single_energy_cfm = 1;
	orig_cfg.spectrum_depths[1] = 2;
	orig_cfg.spectrum_depths[2] = 2;
	orig_cfg.include_h_losses = true; 
	std::vector<int> hlayer_numnodes(3), act_ids(3);
	hlayer_numnodes[0] = 4; hlayer_numnodes[1] = 2; hlayer_numnodes[2] = 1;
	act_ids[0] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[1] = RELU_NN_ACTIVATION_FUNCTION; 
	act_ids[2] = LINEAR_NN_ACTIVATION_FUNCTION;	//Final theta should be linear
	orig_cfg.theta_nn_hlayer_num_nodes = hlayer_numnodes;
	orig_cfg.theta_nn_layer_act_func_ids = act_ids;
	//orig_cfg.ga_converge_thresh = 0.0001;
	//orig_cfg.em_converge_thresh = 0.0001;

	//Feature Calculator
	std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( feature_cfg_file );

	//Prepare some simple data
	std::vector<MolData> data;
	std::string id = "TestMol", smiles = "NCCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";
	data.push_back( MolData( id, smiles, 0 , &orig_cfg ) );
	data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
	data[0].readInSpectraFromFile( spec_file );

	NNParam *final_params;
	for( int energy = 0; energy < 3; energy++ ){

		config_t cfg;
		initSingleEnergyConfig( cfg, orig_cfg, energy );
		std::string param_filename = "tmp_param_output.log";

		//Run EM (multiple times, and take best Q)
		double best_Q = -1000000.0;
		std::vector<double> Qs;
		for( int trial = 0; trial < 5; trial++){
			std::string status_file = "tmp_status_file.log";
			std::string tmp_file = "tmp.log";
			EM_NN em( &cfg, &fc, status_file );
			double Q = em.run( data, 1, tmp_file );
			Qs.push_back(Q);

			if( Q > best_Q ){
				em.writeParamsToFile( param_filename );
				best_Q = Q;
			}
		}
		std::vector<double>::iterator it = Qs.begin();
		for( ; it != Qs.end(); ++it ) std::cout << *it << " ";
		std::cout << " Best=" << best_Q << std::endl;

		if(energy == 0 ) final_params = new NNParam(param_filename);
		else{ 
			NNParam eparam(param_filename);
			final_params->appendNextEnergyParams( eparam, energy );
		}
	}

	//Predict the output spectra
	data[0].computePredictedSpectra( *final_params );
	data[0].postprocessPredictedSpectra( 100.0, 0, 1000 );

	//Compare the original and predicted spectra - should be able to overfit
	//very close to the actual values since training on same (and only same) mol
	for( unsigned int energy = 0; energy < data[0].getNumSpectra(); energy++ ){
		const Spectrum *orig_spec = data[0].getSpectrum(energy);
		const Spectrum *predicted_spec = data[0].getPredictedSpectrum(energy);
		pass &= compareSpectra(orig_spec, predicted_spec, orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, intensity_tol);
		
	}
	passed = pass;
}

class IsotopeTestMol : public MolData {
public:
	IsotopeTestMol(config_t *a_cfg) : MolData("Isotope Test Mol", "CCNCl", a_cfg){

		//Set the spectra
		spectra.resize(1);
		spectra[0].push_back( Peak(80.0261499201, 100.0000000000 ) );
		spectra[0].push_back( Peak(81.0287615888, 2.7114457737 ) );
		spectra[0].push_back( Peak(82.0232058112, 32.4230213142 ) );
		spectra[0].push_back( Peak(63.9948515201, 100.0000000000 ) );
		spectra[0].push_back( Peak(64.9967761924, 1.5320451066 ) );
		spectra[0].push_back( Peak(65.9919024752, 32.4042671242 ) );
		spectra[0].normalizeAndSort();
	}
};

EMTestSingleEnergyIsotopeSelfProduction::EMTestSingleEnergyIsotopeSelfProduction(){
	description = "Test integrated SE-CFM EM ability to learn single spectrum with isotopes";
}

void EMTestSingleEnergyIsotopeSelfProduction::runTest(){
	
	bool pass = true; 
	double intensity_tol = 2.0; //For intensity in range 0-100
    
	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t orig_cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( orig_cfg, param_cfg_file );
	orig_cfg.lambda = 0.0000001;
	orig_cfg.use_single_energy_cfm = 1;
	orig_cfg.spectrum_depths.resize(1);
	orig_cfg.spectrum_weights.resize(1);
	orig_cfg.include_isotopes = 1;
	orig_cfg.use_lbfgs_for_ga = 1;
	orig_cfg.ga_converge_thresh = 0.0001;

	//Feature Calculator
	std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( feature_cfg_file );

	//Prepare some simple data
	std::vector<MolData> data;
	data.push_back( IsotopeTestMol(&orig_cfg) );
	data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, true);

	Param *final_params;

	config_t cfg;
	initSingleEnergyConfig( cfg, orig_cfg, 0 );

	//Run EM
	std::string status_file = "tmp_status_file.log";
	std::string tmp_file = "tmp.log";
	EM em( &cfg, &fc, status_file );
	double Q = em.run( data, 1, tmp_file );
	std::string param_filename = "tmp_param_output.log";
	em.writeParamsToFile( param_filename );
	final_params = new Param(param_filename);

	//Predict the output spectra
	data[0].computePredictedSpectra( *final_params );
	data[0].postprocessPredictedSpectra( 100.0, 0, 1000 );

	//Compare the original and predicted spectra - should be able to overfit
	//very close to the actual values since training on same (and only same) mol
	const Spectrum *orig_spec = data[0].getSpectrum(0);
	const Spectrum *predicted_spec = data[0].getPredictedSpectrum(0);
	pass &= compareSpectra(orig_spec, predicted_spec, orig_cfg.abs_mass_tol, orig_cfg.ppm_mass_tol, intensity_tol);

	passed = pass;
}

EMTestLBFGSvsOriginalGradientAscent::EMTestLBFGSvsOriginalGradientAscent(){
	description = "Test LBFGS gradient ascent vs my original gradient ascent";
}

void EMTestLBFGSvsOriginalGradientAscent::runTest(){
	
	double tol = 0.01;
	bool pass = true; 
	config_t cfg;
	std::string cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( cfg, cfg_file );
	cfg.use_lbfgs_for_ga = 0;
	cfg.ga_converge_thresh = 0.0001;
	cfg.em_init_type = PARAM_FULL_ZERO_INIT;
	cfg.include_h_losses = true; 
	std::string fc_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( fc_file );

	std::string id1 = "TestMol1", id2 = "TestMol2", id3 = "TestMol3";
	std::string smiles1 = "NCCCN", smiles2 = "N=CC(OC)CN", smiles3 = "N(CCCC)CCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";

	std::vector<MolData> data;
	data.push_back( MolData( id1, smiles1, 0, &cfg ) );
	data.push_back( MolData( id2, smiles2, 0, &cfg ) );
	data.push_back( MolData( id3, smiles3, 0, &cfg ) );
	for( int i = 0 ; i < 3; i++ ){
		data[i].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
		data[i].readInSpectraFromFile( spec_file );
	}
	
	std::cout << "Running EM using original gradient ascent algorithm" << std::endl;
	std::string status_file = "tmp_status_file.log";
	std::string tmp_file = "tmp.log";
	EM em1( &cfg, &fc, status_file );
	double orig_Q = em1.run( data, 1, tmp_file );
	std::cout << "Original Q = " << orig_Q << std::endl;

	std::cout << "Running EM using LBFGS gradient ascent algorithm" << std::endl;
	cfg.use_lbfgs_for_ga = 1;
	cfg.ga_converge_thresh = 0.00001;
	EM em2( &cfg, &fc, status_file );
	double lbfgs_Q = em2.run( data, 1, tmp_file );

	std::cout << "Original Q: " << orig_Q << " LBFGS Q: " << lbfgs_Q << std::endl;
	if( orig_Q - lbfgs_Q > tol ){
		std::cout << "Mismatch and/or poorer Q values between running LBFGS vs original gradient ascent" << std::endl;
		pass = false;
	}

	passed = pass;
}

bool runMultiProcessorEMTest( config_t &cfg ){

	bool pass = true;
	double tol = 1e-4;

	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump <= 2 ){
		std::cout << "Error: Test intended for multiple processors" << std::endl;
		return false;
	}

	cfg.em_init_type = PARAM_FULL_ZERO_INIT;
	std::string fc_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( fc_file );


	std::string id1 = "TestMol1", id2 = "TestMol2", id3 = "TestMol3";
	std::string smiles1 = "NCCCN", smiles2 = "N=CC(OC)CN", smiles3 = "N(CCCC)CCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";

	//First Run three molecules on the Master Only
	std::cout << "Running on master only..." << std::endl;
	std::vector<MolData> data;
	if( mpi_rank == MASTER ){
		data.push_back( MolData( id1, smiles1, 0, &cfg ) );
		data.push_back( MolData( id2, smiles2, 0, &cfg ) );
		data.push_back( MolData( id3, smiles3, 0, &cfg ) );
		for( int i = 0 ; i < 3; i++ ){
			data[i].computeFragmentGraphAndReplaceMolsWithFVs(&fc, &cfg);
			data[i].readInSpectraFromFile( spec_file );
		}
	}
	std::string status_file = "tmp_status_file.log";
	std::string tmp_file = "tmp.log";
	EM em1( &cfg, &fc, status_file );
	double all_master_Q = em1.run( data, 1, tmp_file );
	std::cout << "Q = " << all_master_Q << std::endl;

	//Now run the same molecules on 3 separate processes
	std::cout << "Running on 3 processes..." << std::endl;
	data.clear();
	cfg.include_h_losses = true; 
	if( mpi_rank == MASTER ) data.push_back( MolData( id1, smiles1, 0, &cfg ) );
	else if( mpi_rank == 1 ) data.push_back( MolData( id2, smiles2, 0, &cfg ) );
	else if( mpi_rank == 2 ) data.push_back( MolData( id3, smiles3, 0, &cfg ) );
	data[0].computeFragmentGraph(&fc);
	data[0].computeFeatureVectors(&fc);
	data[0].readInSpectraFromFile( spec_file );
	EM em2( &cfg, &fc, status_file );
	double separated_Q = em2.run( data, 1, tmp_file );
	std::cout << "Q = " << separated_Q << std::endl;

	//Check that the output Q values are the same
	std::cout << "Q1: " << all_master_Q << " Q2: " << separated_Q << std::endl;
	if( fabs( all_master_Q - separated_Q ) > tol ){
		std::cout << "Mismatch in Q values between running on one vs many processes" << std::endl;
		pass = false;
	}
	return pass;
}


EMTestMultiProcessor::EMTestMultiProcessor(){
	description = "Test EM running on multiple processors";
}

void EMTestMultiProcessor::runTest(){
	
	bool pass = true;

	//Initialisation
	config_t cfg;
	std::string cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( cfg, cfg_file );
	cfg.use_lbfgs_for_ga = 0;
	runMultiProcessorEMTest( cfg );

	passed = pass;
}

EMTestMultiProcessorLBFGS::EMTestMultiProcessorLBFGS(){
	description = "Test EM running on multiple processors using LBFGS";
}

void EMTestMultiProcessorLBFGS::runTest(){
	
	bool pass = true;

	//Initialisation
	config_t cfg;
	std::string cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( cfg, cfg_file );
	cfg.use_lbfgs_for_ga = 1;
	pass = runMultiProcessorEMTest( cfg );

	passed = pass;
}


EMTestMiniBatchSelection::EMTestMiniBatchSelection(){
	description  = "Test selection of random mini-batches during gradient ascent";
}

void EMTestMiniBatchSelection::runTest(){

	bool pass = true;

	config_t cfg;
	initDefaultConfig(cfg);
	cfg.ga_minibatch_nth_size = 10; //Take 1 in 10
	initDerivedConfig(cfg);

	std::string status_filename = "dummy_status.log";
	std::vector<std::string> fnames;
	fnames.push_back("HydrogenRemoval");
	FeatureCalculator fc(fnames);

	std::vector<int> flags1(1000, 1);	//Select 100
	std::vector<int> flags2(1000, 1);
	//Turn off every second index (to test exclusion of validation molecules)
	for( int i = 0; i < 1000; i += 2 ){
		flags1[i] = 0; flags2[i] = 0; 
	}

	EM em( &cfg, &fc, status_filename );
	em.selectMiniBatch(flags1);
	em.selectMiniBatch(flags2);

	//Check that the random selections select the right number of molecules, 
	//and that the selected molecules are different in the two runs
	int num_on1 = std::count(flags1.begin(), flags1.end(), 1);
	int num_on2 = std::count(flags2.begin(), flags2.end(), 1);
	if( num_on1 != 100 || num_on2 != 100 ){
		std::cout << "Incorrect number of selected molecules in mini batch" << std::endl;
		pass = false;
	}

	std::vector<int>::iterator it1 = flags1.begin(), it2 = flags2.begin();
	int count_same = 0;
	for( ; it1 != flags1.end(); ++it1, ++it2 ){
		count_same += (*it1)*(*it2);
	}
	if( count_same == num_on1 || count_same == num_on2){
		std::cout << "Found all matching molecules between two sets" << std::endl;
		pass = false;		
	}

	passed = pass;

}


EMTestBiasPreLearning::EMTestBiasPreLearning(){
	description = "Test pre-learning of bias parameter before others";
}

void EMTestBiasPreLearning::runTest(){

	bool pass = true; 
	time_t before, after;

	int mpi_rank, mpi_nump;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nump );

	if( mpi_nump != 1 ){
		std::cout << "Error: Test intended for one processor only" << std::endl;
		passed = false;
		return;
	}

	//Config
	config_t orig_cfg;
	std::string param_cfg_file = "tests/test_data/example_param_config.txt";
	initConfig( orig_cfg, param_cfg_file );
	orig_cfg.lambda = 0.0000001;
	orig_cfg.use_single_energy_cfm = 1;
	orig_cfg.spectrum_depths[1] = 2;
	orig_cfg.spectrum_depths[2] = 2;

	//Feature Calculator
	std::string feature_cfg_file = "tests/test_data/example_feature_config_withquadratic.txt";
	FeatureCalculator fc( feature_cfg_file );

	//Prepare some simple data
	std::vector<MolData> data;
	std::string id = "TestMol", smiles = "NCCCN";
	std::string spec_file = "tests/test_data/example_spectra.txt";
	data.push_back( MolData( id, smiles, 0, &orig_cfg ) );
	data[0].computeFragmentGraphAndReplaceMolsWithFVs(&fc, false);
	data[0].readInSpectraFromFile( spec_file );


	for( int energy = 0; energy < 3; energy++ ){

		config_t cfg;
		initSingleEnergyConfig( cfg, orig_cfg, energy );
		cfg.update_bias_first = 0;
		cfg.em_init_type = PARAM_FULL_ZERO_INIT;

		std::cout << "Running with no separate bias update" << std::endl;
		before = time( NULL );

		//Run EM
		std::string status_file = "tmp_status_file.log";
		std::string tmp_file = "tmp.log";
		EM em( &cfg, &fc, status_file );
		double Q = em.run( data, 1, tmp_file );
		std::string param_filename = "tmp_param_output.log";
		em.writeParamsToFile( param_filename );
		Param *comb_params = new Param(param_filename);
		after = time( NULL );
		time_t tcomb = after - before;
		std::cout << "Time Elapsed = "<< tcomb << " seconds" << std::endl;


		std::cout << "Running with separate bias update" << std::endl;
		before = after;
		cfg.update_bias_first = 1;
		EM em2( &cfg, &fc, status_file );
		Q = em2.run( data, 1, tmp_file );
		em2.writeParamsToFile( param_filename );
		Param *sep_params = new Param(param_filename);
		after = time( NULL );
		time_t tsep = after - before;
		std::cout << "Time Elapsed = "<< tsep << " seconds" << std::endl;

		//Check that the resulting parameter is roughly the same and that it took less time...
		//TODO!


	}

	passed = pass;
}
