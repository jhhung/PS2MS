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
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <RDGeneral/FileParseException.h>
#include <filesystem>

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

class Result
{
  public:
    Result()
    : score_sum (0)
	, rank_sum  (0)
    , num       (0)
    {}
    
	void open(const std::string& filename)
	{
		outfile.open(filename);
	}

    void write_file(boost::shared_ptr<MolData> data, int rank, double score)
    {
        std::lock_guard<std::mutex> lock(mut);
        // data->writePredictedSpectraToMspFileStream(outfile);
		outfile << data->getId() << ":" << rank << "," << score << "\n";
    }

    void add_sum(std::string name, int rank, double score)
    {
        score_sum += score;
		rank_sum += rank;
        ++num;
        std::cout << name << ": " << rank << ", " << score << "\n";
    }

    int size()
    {
        return num;
    }

    double get_avg_score()
    {
        return score_sum / num;
    }

	double get_avg_rank()
    {
        return rank_sum / num;
    }

  private:
    mutable std::mutex mut;
    std::ofstream outfile;
    double score_sum;
	long long int rank_sum;
    int num;
};

Spectrum normalization(const Spectrum &target){
	Spectrum spec;
	//Compute the normalizer
	double max = 0.0;
	std::vector<Peak>::const_iterator itp = target.begin();
	for( ; itp != target.end(); ++itp )
		if (itp->intensity > max)
			max = itp->intensity;
	double norm = 1.0;
	if( max > 0.0 ) norm = 100.0/max;

	//Adjust the values
	for( itp = target.begin(); itp != target.end(); ++itp )
	{
		// std::cout << itp->intensity << "->" << itp->intensity * norm << std::endl;
		spec.push_back(Peak(itp->mass, itp->intensity * norm));
	}
	return spec;
}

std::pair<std::string, Spectrum> read_a_spectrum(std::ifstream& infile, bool id_is_rt = false)
{
	std::string line, upper_line, id;
	int num_peaks = 0;
	Spectrum spectrum;
	std::vector<std::string> split;
	while (getline(infile, line))
	{
		if (line.size() < 3)
		{
			if (num_peaks != 0)
				return std::make_pair(id, spectrum);
			continue;
		}

		if (num_peaks != 0)
		{
			boost::split(split, line, boost::is_any_of(" ;"));
			spectrum.push_back(Peak(std::stod(split[0]), std::stod(split[1])));
			continue;
		}
		
		upper_line = boost::to_upper_copy(line);
		if (!id_is_rt && upper_line.substr(0, 2) == "ID")
			id = line.substr(4);
		else if (id_is_rt && upper_line.substr(0, 2) == "RT")
			id = "RT=" + line.substr(4);
		else if (upper_line.substr(0, 9) == "NUM PEAKS")
			num_peaks = std::stoi(line.substr(11));
	}
	return std::make_pair(id, spectrum);
}

std::pair<int, double> compare_candidate(std::string& database_filename, std::string& name, Spectrum& answer, Spectrum& predicted, DotProduct* comparator)
{
	double score;
	std::string max_name;
	std::ifstream candidates(database_filename);
	std::pair<std::string, Spectrum> candidate;
	double threshold = comparator->computeScore(&answer, &predicted);
	int rank = 1;
	// read candidates and compare with target
	while (true)
	{
		candidate = read_a_spectrum(candidates);
		if (candidate.first.empty())
			break;
		candidate.second.roundPeaksToInteger(); // normalizeAndSort
		score = comparator->computeScore(&candidate.second, &predicted);
		if (name == candidate.first)
			continue;
		if (score > threshold)
			++rank;
	}
	return std::make_pair(rank, threshold);
}

Result result;
ThreadPool<void> thread_pool;
std::atomic<bool> quit(false);
std::mutex mut;
void predict_and_record(boost::shared_ptr<MolData> data, config_t* cfg, NNParam* nn_param, DotProduct* comparator, std::string name, double prob_thresh_for_prune, int do_annotate, int apply_postprocessing, int suppress_exceptions, std::string database_filename, std::ofstream& outfile)
{
	/*
	Param *param;
	LikelyFragmentGraphGenerator *fgen;
	try{
		//Calculate the pruned FragmentGraph
		fgen = new LikelyFragmentGraphGenerator(nn_param, cfg, prob_thresh_for_prune ); 
		data->computeLikelyFragmentGraphAndSetThetas(*fgen, prob_thresh_for_prune, do_annotate);

		//Predict the spectra (and post-process, use existing thetas)
		data->computePredictedSpectra( *param, apply_postprocessing, true );
		delete fgen;
	}
	catch( RDKit::MolSanitizeException e ){
		std::cout << "Could not sanitize input: " << data->getSmilesOrInchi() << std::endl;
		if(!suppress_exceptions) throw SpectrumPredictionException("RDKit could not sanitize input: " + data->getSmilesOrInchi());
		delete fgen;
		return;
	}
	catch( RDKit::SmilesParseException pe ){
		std::cout << "Could not parse input: " << data->getSmilesOrInchi() << std::endl;
		if(!suppress_exceptions) throw SpectrumPredictionException("RDKit could not parse input: " + data->getSmilesOrInchi());
		delete fgen;
		return;
	}
	catch( FragmentGraphGenerationException ge ){
		std::cout << "Could not compute fragmentation graph for input: " << data->getSmilesOrInchi() << std::endl;
		if(!suppress_exceptions) throw SpectrumPredictionException("Could not compute fragmentation graph for input: " + data->getSmilesOrInchi());
		delete fgen;
		return;
	}
	catch( IonizationException ie ){
		std::cout << "Could not ionize: " << data->getSmilesOrInchi() << std::endl;		
		if(!suppress_exceptions) throw IonizationException();
		delete fgen;
		return;
	}
	catch( FragmentGraphTimeoutException ie ){
		std::cout << "Unable to proceed: " << data->getSmilesOrInchi() << ", Molecule fragment graph computation exceeded timeout" << std::endl;		
		if(!suppress_exceptions) throw FragmentGraphTimeoutException();
		delete fgen;
		return;
	}
	catch (std::exception& e)
	{
		std::cout << "Exception occured: " << e.what() << std::endl;
		delete fgen;
		return;
	}
	result.write_file(data);
	*/
	Spectrum /*s = *(data->getSpectrum(0)),*/ t;
	try
	{
		LikelyFragmentGraphGenerator *fgen = new LikelyFragmentGraphGenerator(nn_param, cfg, prob_thresh_for_prune );
		data->computeLikelyFragmentGraphAndSetThetas(*fgen, prob_thresh_for_prune, do_annotate);
		const FragmentGraph* fg = data->getFragmentGraph();
		// std::vector<bool> is_children(fg->getNumFragments(), true);
		// for (unsigned int i = 0; i < fg->getNumTransitions(); ++i)
		// 	is_children[fg->getTransitionAtIdx(i)->getFromId()] = false;
		// for (unsigned int i = 0; i < is_children.size(); ++i)
		// 	if (is_children[i])
		// 		s.push_back(Peak())
		for (unsigned int i = 0; i < fg->getNumFragments(); ++i)
			t.push_back(Peak(fg->getFragmentAtIdx(i)->getMass(), 0));
		delete fgen;
	}
	catch (...) {}
	// s.roundPeaksToInteger();
	t.roundPeaksToInteger();
	std::lock_guard<std::mutex> lock(mut);
	// s.outputToMspStream(outfile, data->getId(), cfg->ionization_mode, 0);
	t.outputToMspStream(outfile, data->getId(), cfg->ionization_mode, 0);
	// auto r = compare_candidate(database_filename, name, s, t, comparator);
	// result.add_sum(name, std::get<0>(r), std::get<1>(r)); //modify
	// result.write_file(data, std::get<0>(r), std::get<1>(r));
}

bool find_answer(boost::shared_ptr<MolData> data, std::ifstream& database)
{
	bool ignore = false;
	std::string line, upper_line, name = data->getId();
	while (getline(database, line))
	{
		if (line.size() < 3)
		{
			if (ignore)
				ignore = false;
			continue;
		}
		if (!ignore)
		{
			upper_line = boost::to_upper_copy(line);
			if (upper_line.substr(0, 2) == "ID")
			{
				if (line.substr(4) != name)
					ignore = true;
			}
			else if (upper_line.substr(0, 9) == "NUM PEAKS")
			{
				int num_peaks = std::stoi(line.substr(11));
				data->readInSpectraFromMSPFileStream(database, num_peaks);
				return true;
			}
		}
	}
	return false;
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
	int apply_postprocessing = 0;
	int suppress_exceptions = 1;
	std::string output_filename = argv[2];
	std::string param_filename = "param_output.log";
	std::string config_filename = "param_config.txt";
	std::string database_filename = argv[1];
	double prob_thresh_for_prune = 0.001;
	std::signal(SIGINT, got_signal);
	// result.open(argv[2]);

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
	if( !std::filesystem::exists( config_filename ) ){
		std::cout << "Could not find file: " <<  config_filename << std::endl;
		throw FileException("Could not find file: " + config_filename);
	}
	initConfig( cfg, config_filename );

	//Read in the parameters
	if( !std::filesystem::exists( param_filename ) ){
		std::cout << "Could not find file: " <<  param_filename << std::endl;
		throw FileException("Could not find file: " + param_filename);
	}
	Param *param; NNParam *nn_param;
	if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
		nn_param = new NNParam( param_filename );
	else
		param = new Param(param_filename);

	std::ifstream database_library(database_filename);
	std::ofstream outfile(output_filename);
	std::string line, upper_line, smiles, name, sdf_name;
	int num_peaks = 0, num_predicted = 0;
	bool ignore = false;
	RDKit::ROMOL_SPTR mol;
	boost::shared_ptr<MolData> data;
	// OrigSteinDotProduct comparator(1000.0, 1.0); // 10.0, 0.01
	DotProduct comparator(1000.0, 1.0);
	double score, sum = 0, ignore_num = 0;
	while (getline(database_library, line) && !quit.load())
	{
		std::cout << ".." << line << std::endl;
		data.reset(new MolData(line, line, 0, &cfg));
		thread_pool.submit(0, std::bind(predict_and_record, data, &cfg, nn_param, &comparator, line, prob_thresh_for_prune, do_annotate, apply_postprocessing, suppress_exceptions, database_filename, std::ref(outfile)));
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
	std::cout << "Can't be parsed: " << ignore_num << std::endl;
	// std::cout << "Avg. score: " << result.get_avg_score() << std::endl;
	// std::cout << "Avg. rank: " << result.get_avg_rank() << std::endl;
	return 0;
}
