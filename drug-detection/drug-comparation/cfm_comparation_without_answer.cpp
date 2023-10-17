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
#include "Spectrum.h"
#include "thread_pool.hpp"

#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <RDGeneral/FileParseException.h>
#include <boost/filesystem.hpp>

int main(int argc, char *argv[]);

#include <iostream>
#include <fstream>
#include <string>
#include <csignal>
#include <cmath>
#include <algorithm>
#include <bitset>
#include <set>
#include <sys/stat.h>

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

        void write_file(boost::shared_ptr<MolData> data, std::vector<std::pair<double, std::string>>& similar_list)
        {
            std::lock_guard<std::mutex> lock(result_mut);
            int cnt = 0;
            outfile << data->getId() << "\n";
            for (auto& item : similar_list) {
                outfile << item.first << ": " << item.second << "; ";
                ++cnt;
                if (cnt%5==0)
                    outfile << "\n";
            } 
            if (similar_list.size()%5!=0)
              outfile << "\n";
            outfile << "\n";
        }

        int size()
        {
            return num;
        }
    
    private:
        mutable std::mutex result_mut;
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
        spec.push_back(Peak(itp->mass, itp->intensity * norm));
    }
    return spec;
}

void filter_peak(Spectrum &target, double hight_threshold = 0.05, unsigned int num_threshold = 10)
{
    std::vector<Peak> peaks(*target.getPeaks());

    int intensity_sum = 0;
    for (const auto& peak : peaks)
        intensity_sum += peak.intensity;
    hight_threshold *= intensity_sum;

    int end_idx = (num_threshold < peaks.size()) ? num_threshold : peaks.size();
    std::sort(peaks.begin(), peaks.end(), [](const Peak& left, const Peak& right){return left.intensity > right.intensity;});
    while (end_idx < peaks.size() && peaks[end_idx].intensity > hight_threshold)
        ++end_idx;

    target.clear();
    for (int i = 0; i < end_idx; ++i)
        target.push_back(peaks[i]);
}

std::tuple<std::string, Spectrum, std::string> read_a_spectrum(std::ifstream& infile, bool id_is_rt = false)
{
    std::string line, upper_line, id, fp;
    int num_peaks = 0;
    Spectrum spectrum;
    std::vector<std::string> split;
    while (getline(infile, line))
    {
        if (line.size() < 3)
        {
            if (num_peaks != 0)
            {
                return std::make_tuple(id, spectrum, fp);
            }
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
        else if (upper_line.substr(0, 2) == "FP")
            fp = line.substr(4);
        else if (upper_line.substr(0, 9) == "NUM PEAKS")
            num_peaks = std::stoi(line.substr(11));
    }
    return std::make_tuple(id, spectrum, fp);
}

double fp_similarity(std::string fp_str, std::string predicted_fp_str)
{
  int fp_size = (fp_str.size() >= predicted_fp_str.size()) ? fp_str.size() : predicted_fp_str.size() ;
  const int MAX_SIZE = 4096;
	fp_str.resize(MAX_SIZE, '0');
	predicted_fp_str.resize(MAX_SIZE, '0');
	std::bitset<MAX_SIZE> fp(fp_str);
	std::bitset<MAX_SIZE> predicted_fp(predicted_fp_str);
  std::bitset<MAX_SIZE> ans = ~(fp ^ predicted_fp);
  std::size_t ret = ans.count() - (MAX_SIZE-fp_size);
  return static_cast<double>(ret) / static_cast<double>(fp_size);
}

double jaccard_similarity(std::string fp_str, std::string predicted_fp_str)
{
    const int MAX_SIZE = 4096;
    fp_str.resize(MAX_SIZE, '0');
    predicted_fp_str.resize(MAX_SIZE, '0');
    std::bitset<MAX_SIZE> fp(fp_str);
    std::bitset<MAX_SIZE> predicted_fp(predicted_fp_str);
    std::size_t intersection = (fp & predicted_fp).count();
    return static_cast<double>(intersection) / (fp.count() + predicted_fp.count() - intersection);
}

void find_peaks(Spectrum &spec, std::set<int> &top_peaks, int &top_n) 
{
    spec.sortedByIntensity();
    for (int i=0;i<top_n;++i) {
        if (i>=spec.size()) break;
        top_peaks.insert(spec.getPeak(i)->mass);
    }
    spec.sortedByMass();
}

double jaccard_similarity(Spectrum &a, Spectrum &b, int &top_n)
{
    std::set<int> a_top_peaks;
    std::set<int> b_top_peaks;
    find_peaks(a, a_top_peaks, top_n);
    find_peaks(b, b_top_peaks, top_n);

    std::set<int> intersection;
    set_intersection(a_top_peaks.begin(), a_top_peaks.end(),
            b_top_peaks.begin(), b_top_peaks.end(),
            std::inserter(intersection, intersection.begin()));

    double size_inter = intersection.size();
    double size_a = a_top_peaks.size();
    double size_b = b_top_peaks.size();

    double jaccard = size_inter/(size_a+size_b-size_inter);
    return jaccard;
}

std::mutex sort_mut;

std::tuple<std::vector<std::pair<double, std::string>>> compare_candidate
    (std::string database_filename, std::string& name, Spectrum& predicted, std::string& predicted_fp, DotProduct* comparator
    , bool restrict_mw, int &top_n)
{
    double score;
    std::vector<std::pair<double, std::string>> similar_list;
    RDKit::RWMOL_SPTR mol(RDKit::SmilesToMol(name));
    double mw = RDKit::Descriptors::calcExactMW(*mol);
    bool include_fp = (predicted_fp.size() != 0);
    bool spectrum_jaccard = (top_n!=0);
    double thres = 0;
    //double spectrum_score, cosine_score, jaccard_score, fp_score;

    // read candidates and compare with target
    std::ifstream candidates(database_filename);
    RDKit::RWMOL_SPTR candidate_mol;
    while (true)
    {
        auto [candidate_id, candidate_spec, candidate_fp] = read_a_spectrum(candidates);
        if (candidate_id.empty())
            break;
        candidate_mol.reset(RDKit::SmilesToMol(candidate_id));
        if (restrict_mw && std::abs(mw - RDKit::Descriptors::calcExactMW(*candidate_mol)) > 1)
            continue;

        candidate_spec.roundPeaksToInteger(); // normalizeAndSort
        double candidate_score;
        if (spectrum_jaccard) {
            candidate_score = pow((comparator->computeScore(&candidate_spec, &predicted) * jaccard_similarity(candidate_spec, predicted, top_n)), (1.0/2.0));
        }
        else
            candidate_score = comparator->computeScore(&candidate_spec, &predicted);
        if (candidate_fp.size() != 0 && predicted_fp.size() != 0) {
            score = (0.7*candidate_score + 0.3*jaccard_similarity(candidate_fp, predicted_fp));
        }
        else
            score = candidate_score;


        int rank_n = 100;
        if (similar_list.size() == rank_n && score>thres) // brian modify
        {
            similar_list.emplace_back(std::make_pair(score, candidate_id));
            std::sort(similar_list.begin(), similar_list.end(), [](auto& left, auto& right) {return left.first > right.first; });
            similar_list.pop_back();
            thres = similar_list[rank_n-1].first;
        } else if (similar_list.size()<rank_n) {
            similar_list.emplace_back(std::make_pair(score, candidate_id));
            std::sort(similar_list.begin(), similar_list.end(), [](auto& left, auto& right) {return left.first > right.first; });
        }
    } 
    return std::make_tuple(similar_list);
}

Result result;
ThreadPool<void> thread_pool;
std::atomic<bool> quit(false);
std::mutex mut;
void predict_and_record(boost::shared_ptr<MolData> data, config_t* cfg, DotProduct* comparator, std::string name, 
  std::string database_filename, bool restrict_mw, int &top_n)
{
    Spectrum t = *(data->getPredictedSpectrum(0));
    t.roundPeaksToInteger();

    auto r = compare_candidate(database_filename, name, t, data->predicted_fp, comparator, restrict_mw, top_n);

    result.write_file(data, std::get<0>(r));
}

void got_signal(int signal)
{
    std::cerr << "Signal raised: " << signal << "\n";
    quit.store(true);
}

int main(int argc, char *argv[])
{
    std::string config_filename = "param_config.txt";
    std::string msp_filename = argv[2]; // analyte; unknown spectrum; predicted fingerprint
    std::string database_filename = argv[1]; // library; enumerated smiles; predicted spectrum; calculated fingerprint
    std::signal(SIGINT, got_signal);
    result.open(argv[3]);
    bool restrict_mw = false;
    if (argc < 5 || std::strcmp(argv[4], "false") == 0)
        restrict_mw = false;
    else if (std::strcmp(argv[4], "true") == 0)
        restrict_mw = true;
    bool spectrum_jaccard = false;
    int top_n = 0;
    if( argc < 6 || std::atoi(argv[5])==0)
        spectrum_jaccard = false;
    else {
        spectrum_jaccard = true;
        top_n = std::atoi(argv[5]);
    }

    //Initialise model configuration
    config_t cfg;
    if( !boost::filesystem::exists( config_filename ) ){
        std::cout << "Could not find file: " <<  config_filename << std::endl;
        throw FileException("Could not find file: " + config_filename);
    }
    initConfig( cfg, config_filename );

    std::ifstream msp_library(msp_filename);
    std::string line, upper_line, smiles, name, sdf_name, fp;
    int num_peaks = 0, num_predicted = 0, mass_weight = 0;
    bool ignore = false;
    RDKit::ROMOL_SPTR mol;
    boost::shared_ptr<MolData> data;
    DotProduct comparator(0.0, 1.1);
    double score, sum = 0, ignore_num = 0;
    while (getline(msp_library, line) && !quit.load())
    {
        if (line.size() < 3)
        {
            if (ignore)
                ignore = false;
            else if (data)
            {
              thread_pool.submit(0, std::bind(predict_and_record, data, &cfg, &comparator, name, database_filename, restrict_mw, top_n));
            }
            smiles.clear();
            name.clear();
            fp.clear();
            data.reset();
            database_filename.clear();
            num_peaks = 0;
            continue;
        }
        if (!ignore)
        {
            upper_line = boost::to_upper_copy(line);
            if (upper_line.substr(0, 2) == "ID")
                name = line.substr(4);
            if (upper_line.substr(0, 2)=="MW") {
                mass_weight = int(std::stod(line.substr(4)));
            }
            if (name == "CC")
            {
                ignore = true;
                for (int i = 0; i < 3; ++i)
                    getline(msp_library, line);
            }
            else if (upper_line.substr(0, 2) == "FP")
            {
                fp = line.substr(4);
            }
            else if (upper_line.substr(0, 9) == "NUM PEAKS")
            {
                try
                {
                    num_peaks = std::stoi(line.substr(11));
                }
                catch(const std::exception& e)
                {
                    std::cerr << e.what() << ": " << name << '\n';
                    exit(0);
                }
                data.reset(new MolData(name, name, 0, &cfg));
                if (!fp.empty())
                    data->predicted_fp = fp;
                data->readInSpectraFromMSPFileStream(msp_library, num_peaks, true);
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
    std::cout << "Can't be parsed: " << ignore_num << std::endl;
    return 0;
}
