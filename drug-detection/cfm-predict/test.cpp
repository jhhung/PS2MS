/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# main.cpp
#
# Description: 	Predict the MS/MS spectra for a given structure using a
#				pre-trained CFM model.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Config.h"
#include "Param.h"
#include "Features.h"
#include "MolData.h"
#include "Comparators.h"
#include "thread_pool.hpp"

#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <boost/filesystem.hpp>

int main(int argc, char *argv[]);

#include <iostream>
#include <fstream>
#include <string>
#include <csignal>

#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>

class SpectrumPredictionException: public std::exception{
private:
    std::string message_;
public:
	SpectrumPredictionException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << message_ << std::endl;
		return message_.c_str();
	}
	~SpectrumPredictionException() throw() {};
};

class FileException: public std::exception{
private:
    std::string message_;
public:
	FileException(const std::string& message) throw() : message_(message) {};
	virtual const char* what() const throw(){
		std::cout << message_ << std::endl;
		return message_.c_str();
	}
	~FileException() throw() {};
};

Result result("output.msp");
ThreadPool thread_pool;
std::atomic<bool> quit(false);
void predict_and_record(boost::shared_ptr<MolData> data, config_t* cfg, NNParam* nn_param, OrigSteinDotProduct* comparator, std::string name, double prob_thresh_for_prune, int do_annotate, int apply_postprocessing, int suppress_exceptions)
{
	Param *param;
	try{
		//Calculate the pruned FragmentGraph
		LikelyFragmentGraphGenerator *fgen;
		fgen = new LikelyFragmentGraphGenerator(nn_param, cfg, prob_thresh_for_prune ); 
		data->computeLikelyFragmentGraphAndSetThetas(*fgen, prob_thresh_for_prune, do_annotate);

		//Predict the spectra (and post-process, use existing thetas)
		data->computePredictedSpectra( *param, apply_postprocessing, true );
	}
	catch( RDKit::MolSanitizeException e ){
		std::cout << "Could not sanitize input: " << data->getSmilesOrInchi() << std::endl;
		if(!suppress_exceptions) throw SpectrumPredictionException("RDKit could not sanitize input: " + data->getSmilesOrInchi());
		return;
	}
	catch( RDKit::SmilesParseException pe ){
		std::cout << "Could not parse input: " << data->getSmilesOrInchi() << std::endl;
		if(!suppress_exceptions) throw SpectrumPredictionException("RDKit could not parse input: " + data->getSmilesOrInchi());
		return;
	}
	catch( FragmentGraphGenerationException ge ){
		std::cout << "Could not compute fragmentation graph for input: " << data->getSmilesOrInchi() << std::endl;
		if(!suppress_exceptions) throw SpectrumPredictionException("Could not compute fragmentation graph for input: " + data->getSmilesOrInchi());
		return;
	}
	catch( IonizationException ie ){
		std::cout << "Could not ionize: " << data->getSmilesOrInchi() << std::endl;		
		if(!suppress_exceptions) throw IonizationException();
		return;
	}
	catch( FragmentGraphTimeoutException ie ){
		std::cout << "Unable to proceed: " << data->getSmilesOrInchi() << ", Molecule fragment graph computation exceeded timeout" << std::endl;		
		if(!suppress_exceptions) throw FragmentGraphTimeoutException();
		return;
	}
	result.write_file(data);
	result.add_sum(name, comparator->computeScore( data->getSpectrum(0), data->getPredictedSpectrum(0) ));
}

void got_signal(int signal)
{
    std::cerr << "Signal raised: " << signal << "\n";
    quit.store(true);
}

int main(int argc, char *argv[])
{
	bool to_stdout = true;
	int do_annotate = 0;
	int apply_postprocessing = 1;
	int suppress_exceptions = 1;
	std::string output_filename = "output.msp";
	std::string param_filename = "param_output.log";
	std::string config_filename = "param_config.txt";
	double prob_thresh_for_prune = 0.001;
	std::signal(SIGINT, got_signal);

	/*
	if (argc != 6 && argc != 2 && argc != 5 && argc != 3 && argc != 7 && argc != 8 && argc != 9)
	{
		std::cout << std::endl << std::endl;
		std::cout << std::endl << "Usage: cfm-predict.exe <input_smiles_or_inchi> <prob_thresh_for_prune> <param_filename> <config_filename> <include_annotations> <output_filename> <apply_post_processing>" << std::endl << std::endl << std::endl;
		std::cout << std::endl << "input_smiles_or_inchi_or_file:" << std::endl << "The smiles or inchi string of the structure whose spectra you want to predict, or a .txt file containing a list of <id smiles> pairs, one per line." << std::endl;
		std::cout << std::endl << "prob_thresh_for_prune (opt):" << std::endl << "The probability below which to prune unlikely fragmentations (default 0.001)" << std::endl;
		std::cout << std::endl << "param_filename (opt):" << std::endl << "The filename where the parameters of a trained cfm model can be found (if not given, assumes param_output.log in current directory)" << std::endl;
		std::cout << std::endl << "config_filename (opt):" << std::endl << "The filename where the configuration parameters of the cfm model can be found (if not given, assumes param_config.txt in current directory)" << std::endl;
		std::cout << std::endl << "include_annotations (opt):" << std::endl << "Whether to include fragment information in the output spectra (0 = NO (default), 1 = YES ). Note: ignored for msp/mgf output." << std::endl;
		std::cout << std::endl << "output_filename_or_dir (opt):" << std::endl << "The filename of the output spectra file to write to (if not given, prints to stdout), OR directory if multiple smiles inputs are given (else current directory) OR msp or mgf file." << std::endl;
		std::cout << std::endl << "apply_postprocessing (opt):" << std::endl << "Whether or not to post-process predicted spectra to take the top 80% of energy (at least 5 peaks), or the highest 30 peaks (whichever comes first) (0 = OFF, 1 = ON (default) )." << std::endl;
		std::cout << std::endl << "suppress_exception (opt):" << std::endl << "Suppress exceptions so that the program returns normally even when it fails to produce a result (0 = OFF (default), 1 = ON)." << std::endl;
		exit(1);
	}*/
	
	//Initialise model configuration
	config_t cfg;
	if( !boost::filesystem::exists( config_filename ) ){
		std::cout << "Could not find file: " <<  config_filename << std::endl;
		throw FileException("Could not find file: " + config_filename);
	}
	initConfig( cfg, config_filename );

	//Read in the parameters
	if( !boost::filesystem::exists( param_filename ) ){
		std::cout << "Could not find file: " <<  param_filename << std::endl;
		throw FileException("Could not find file: " + param_filename);
	}
	Param *param; NNParam *nn_param;
	if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
		nn_param = new NNParam( param_filename );
	else
		param = new Param(param_filename);

	std::ifstream msp_library("GMD.msp");
	std::ofstream outfile(output_filename);
	std::string line, upper_line, smiles, name;
	int num_peaks = 0, num_predicted = 0;
	bool ignore = false;
	boost::shared_ptr<MolData> data;
	OrigSteinDotProduct comparator(10.0, 0.01);
	double score, sum = 0;
	while (getline(msp_library, line))
	{
		if (line.size() < 3)
		{
			if (ignore)
				ignore = false;
			else if (data)
				thread_pool.submit(0, std::bind(predict_and_record, data, &cfg, nn_param, &comparator, name, prob_thresh_for_prune, do_annotate, apply_postprocessing, suppress_exceptions));
			smiles.clear();
			name.clear();
			data.reset();
			num_peaks = 0;
			continue;
		}
		if (!ignore)
		{
			upper_line = boost::to_upper_copy(line);
			if (upper_line.substr(0, 4) == "NAME")
				name = line.substr(6);
			else if (upper_line.substr(0, 13) == "SPECTRUM_TYPE")
			{
				if (upper_line.substr(15) != "MS1")
					ignore = true;
			}
			// else if (upper_line.substr(0, 15) == "INSTRUMENT_TYPE")
			// 	if (upper_line.find("EI") == std::string::npos)
			// 		ignore = true;
			// else if (upper_line.substr(0, 8) == "ION_MODE")
			// 	if (upper_line.substr(10) != "P")
			// 		ignore = true;
			else if (upper_line.substr(0, 8) == "COMMENTS")
			{
				std::vector<std::string> comments;
				std::string c = line.substr(10);
				boost::split(comments, c, boost::is_any_of("\""));
				for (std::size_t i = 1; i < comments.size(); i += 2)
				{
					if (boost::to_upper_copy(comments[i]).substr(0, 15) == "COMPUTED SMILES")
						smiles = comments[i].substr(16);
					else if (boost::to_upper_copy(comments[i]).substr(0, 6) == "SMILES")
						smiles = comments[i].substr(7);
				}
			}
			else if (upper_line.substr(0, 9) == "NUM PEAKS")
			{
				num_peaks = std::stoi(line.substr(11));
				data.reset(new MolData(name, smiles, 0, &cfg));
				data->readInSpectraFromMSPFileStream(msp_library, num_peaks);
			}
		}
	}
	unsigned int remain = thread_pool.get_task_num();
    std::cerr << "Task add finished, size = " << remain << std::endl;
    while (remain != 0 && !quit.load())
    {
        std::this_thread::sleep_for(std::chrono::seconds(10));
        remain = thread_pool.get_task_num();
        std::cerr << "Task remain: " << remain << "\n";
    }
    std::cerr << "Wait threads to complete remaining tasks ...\n";
    thread_pool.terminate_all_thread();
	std::cout << "Test cases: " << result.size() << std::endl;
	std::cout << "Avg. score: " << result.get_avg_score() << std::endl;
	return 0;
}