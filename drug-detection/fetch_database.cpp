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

#include "MolData.h"
#include "Comparators.h"
#include "thread_pool.hpp"

#include <GraphMol/SanitException.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <filesystem>

int main(int argc, char *argv[]);

#include <iostream>
#include <fstream>
#include <string>
#include <csignal>

#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>

std::mutex print_mut;
void print_result(std::pair<std::string, std::vector<std::pair<std::string, double>>> result)
{
	std::lock_guard<std::mutex> lock(print_mut);
	std::cerr << "Possible candidate of " << result.first << ":\n";
	for (auto& item : result.second)
		std::cerr << item.first << ": " << item.second << "\n";
	std::cerr << "\n";
}

ThreadPool<std::pair<std::string, std::vector<std::pair<std::string, double>>>> thread_pool(print_result);
std::atomic<bool> quit(false);

void got_signal(int signal)
{
    std::cerr << "Signal raised: " << signal << "\n";
    quit.store(true);
}

std::pair<std::string, Spectrum> read_a_spectrum(std::ifstream& infile, bool include_rt = false)
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
		if (upper_line.substr(0, 4) == "NAME")
			id = line.substr(6);
		else if (include_rt && upper_line.substr(0, 2) == "RT")
			id += ", RT=" + line.substr(4);
		else if (upper_line.substr(0, 9) == "NUM PEAKS")
			num_peaks = std::stoi(line.substr(11));
	}
	return std::make_pair(id, spectrum);
}

std::pair<std::string, std::vector<std::pair<std::string, double>>> 
compare_candidate(std::string& database_filename, std::pair<std::string, Spectrum> target, std::string& path, OrigSteinDotProduct& comparator, double threshold)
{
	double score;
    if (!std::filesystem::exists(path))
        std::filesystem::create_directory(path);
	std::ifstream candidates(database_filename);
	std::ofstream outfile(path + "/" + target.first + ".txt");
	std::pair<std::string, Spectrum> candidate;
	std::vector<std::pair<std::string, double>> result;
	// read candidates and compare with target
	target.second.normalizeAndSort();
	while (true)
	{
		candidate = read_a_spectrum(candidates);
		if (candidate.first.empty())
			break;
		candidate.second.roundPeaksToInteger(); // normalizeAndSort
		score = comparator.computeScore(&candidate.second, &target.second);
		// std::cerr << score << std::endl;
		if (score > threshold)
        {
			result.emplace_back(std::make_pair(candidate.first, score));
            outfile << candidate.first << "," << score << "\n";
        }
	}
    // std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs){ return lhs.second < rhs.second; });
	return std::make_pair(target.first, result);
}

int main(int argc, char *argv[])
{
	double threshold = 0;
	std::string sample_filename = argv[1];
	std::string database_filename = argv[2];
    std::string path = argv[3];
	std::signal(SIGINT, got_signal);

	std::ifstream samples(sample_filename);
	OrigSteinDotProduct comparator(1000.0, 1.0); // 10.0, 0.01
	std::pair<std::string, Spectrum> target;
	while (true)
	{
		target = read_a_spectrum(samples, true);
		if (target.first.empty())
			break;
		thread_pool.submit(0, std::bind(compare_candidate, std::ref(database_filename), target, std::ref(path), std::ref(comparator), threshold));
		// compare_candidate(database_filename, target, comparator, threshold);
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
	return 0;
}
