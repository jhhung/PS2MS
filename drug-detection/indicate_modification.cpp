#include <iostream>
#include <unordered_set>
#include <map>
#include <fstream>
#include <iomanip>
#include <csignal>
#include <algorithm>

#include "thread_pool.hpp"
#include "DruglikeFilter.hpp"
#include "modifications.hpp"

#include "Config.h"
#include "Param.h"
#include "Features.h"
#include "MolData.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SanitException.h>
using namespace RDKit; 

#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <gperftools/profiler.h>

std::vector<std::string> exceptions;
constexpr int MOD_DEPTH = 3;
constexpr int WRITE_FREQ = 100000;
std::atomic<bool> quit(false);
ThreadPool<std::pair<int, std::string>> thread_pool;
config_t cfg;
Param *param; NNParam *nn_param;
std::vector<ModificationType> mod_list = {
    make_haloalkane("F"),
    make_haloalkane("Cl"),
    make_haloalkane("Br"),
    make_haloalkane("I")
};

using FILTER_SPTR = boost::shared_ptr<Filter>;
using TripleIntVec = std::vector<std::vector<std::vector<uint8_t>>>;
using H_INDEX_LIST_SPTR = boost::shared_ptr<TripleIntVec>;

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

struct Result
{
  public:
    Result() : size(0)
    {
        if (std::filesystem::exists("result.txt")
            || std::filesystem::exists("filtered_smiles.txt"))
        {
            std::cerr << "[ERROR] Result file still exists.\n";
            // exit(-1); // modified
        }
        outfile.open("result.txt");
    }

    ~Result()
    {
        output_result();
        outfile.close();
    }

    bool find(std::string& key)
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return candidate.find(key) != candidate.end()
                || filtered_mol.find(key) != filtered_mol.end();
    }

    unsigned int get_candidate_size()
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return candidate.size();
    }

    unsigned int get_filtered_size()
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return filtered_mol.size();
    }

    void append_smiles(RWMOL_SPTR mol, std::string& h_smiles, std::string& reason)
    {
        std::unique_lock<std::shared_mutex> lock(mutex);
        // if (candidate.find(h_smiles) != candidate.end()
        //     || filtered_mol.find(h_smiles) != filtered_mol.end())
        {
            if (reason.empty())
            {
                candidate.emplace(h_smiles);
                std::cout << "candidate" << std::endl;
                outfile << MolToSmiles(*MolOps::removeHs(ROMol(*mol))) << "\n";
            }
            else
                filtered_mol.emplace(h_smiles, reason);
        }
        ++size;
        if (size == WRITE_FREQ)
        {
            output_result();
            size = 0;
        }
    }

    // void write_file(MolData* data)
    void write_file(const Spectrum& spec, const std::string& name)
    {
        std::lock_guard<std::mutex> lock(file_mutex);
        // data->writePredictedSpectraToMspFileStream(out_msp);
        spec.outputToMspStream(outfile, name, 0, 0);
    }
    std::ofstream outfile;

  private:
    void output_result()
    {
        // std::ofstream result("result.txt");//, std::ios::app);
        // for (auto& c : candidate)
        //     outfile << c << "\n";
        // result.close();

        std::ofstream outfile("filtered_smiles.txt");//, std::ios::app);
        for (auto& [key, value] : filtered_mol)
            outfile << std::setw(150) << std::left << key << "(" << value << ")\n";
        outfile.close();
        // filtered_mol.clear();
    }

    unsigned int size;
    mutable std::shared_mutex mutex;
    mutable std::mutex file_mutex;
    std::unordered_set<std::string> candidate;
    std::map<std::string, std::string> filtered_mol;
    // std::ofstream out_msp;
};

std::vector<std::vector<uint8_t>> enumerate_h(RWMOL_SPTR mol, std::vector<int>& except_h, int h_num)
{
    std::vector<Atom*> all_h;
    all_h.reserve(mol->getNumAtoms() - except_h.size());
    for (auto& atom : mol->atoms())
        if (atom->getSymbol() == "H" 
            && std::find(except_h.begin(), except_h.end(), atom->getIdx()) == except_h.end())
            all_h.emplace_back(atom);
    std::vector<std::vector<uint8_t>> index_list;
    if (all_h.size() < h_num)
        return index_list;
    
    // get permulation
    std::vector<uint8_t> combination(h_num);
    std::vector<uint8_t> permutation(h_num);
    for (int j = 0; j < h_num; ++j)
        combination[j] = j + 1;
    do {
        index_list.push_back(combination);
    } while (std::next_permutation(combination.begin(), combination.end()));
    std::vector<uint8_t>::iterator first = combination.begin();
    std::vector<uint8_t>::iterator last = combination.end();
    unsigned int n = all_h.size(), r = h_num;
    while((*first) != (int)(n-r+1))
    {
        std::vector<uint8_t>::iterator mt = last;

        while (*(--mt) == n-(last-mt)+1);
        (*mt)++;
        while (++mt != last) *mt = *(mt-1)+1;

        permutation.assign(first, last);
        std::sort(permutation.begin(), permutation.end());
        do {
            index_list.push_back(permutation);
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    
    // map index back to atom_index
    for (auto &item : index_list)
    {   
        for (std::size_t i = 0; i < item.size(); ++i)
            item[i] = all_h[item[i] - 1]->getIdx();
    }
    return index_list;
}

void add_h(RWMOL_SPTR mol)
{
    if (mol->needsUpdatePropertyCache())
        mol->updatePropertyCache();
    MolOps::addHs(*mol);
}

RWMOL_SPTR do_modification(
    std::vector<uint8_t>& h_list, RWMOL_SPTR parent, const std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>& mod)
{
    RWMOL_SPTR result(new RWMol(*parent, true));
    mod(result, h_list);
    // if (result->needsUpdatePropertyCache())
    //     result->updatePropertyCache();
    add_h(result);
    // MolOps::removeHs(*result);
    // add_h(result);
    return result;
}

void recursive_modification(RWMOL_SPTR mol, std::string mol_smiles, std::vector<int>& except_h, Result& result_pool, int depth);

std::pair<int, std::string> replace_h_and_filt(RWMOL_SPTR mol, H_INDEX_LIST_SPTR h_index_list, std::vector<int>& except_h, Result& result_pool, int mod_num, int depth)
{
    ModificationType& mod = mod_list[mod_num];
    std::pair<int, std::string> r = std::make_pair(depth, std::get<0>(mod));
    // std::vector<std::vector<uint8_t>> h_list = enumerate_h(mol, std::get<1>(mod));
    std::vector<std::vector<uint8_t>>& h_list = h_index_list->operator[](std::get<1>(mod));
    unsigned int h_size = h_list.size();
    // unsigned int h_print_interval = h_size / 10 == 0 ? 1 : (unsigned int)(h_size / 10);
    std::string reason;
    for (unsigned int j = 0; j < h_size; ++j)
    {
        for (const auto& mod_func : std::get<2>(mod))
        {
            RWMOL_SPTR result;
            try
            {
                result = do_modification(h_list[j], mol, mod_func);
            }
            catch (...)
            {
                std::cerr << "excetpion occur: " << MolToSmiles(*result) << std::endl;
                continue;
            }
            std::string h_smiles = MolToSmiles(*result);
            if (h_smiles.find('.') != std::string::npos)
            {
                if (std::get<0>(mod).find("bisthiosemicarbazone") == std::string::npos &&
                    std::get<0>(mod).find("methine") == std::string::npos)
                {
                    std::cerr << "found lonely H-H: " << std::get<0>(mod) << std::endl;
                }
                continue;
            }
            if (!result_pool.find(h_smiles))
            {
                // do filter
                // result->setProp("filtered", false);
                // result->setProp("failedfilter", "");
                // reason = (*filter)(result);
                // if (j % h_print_interval == 0)
                // {
                //     std::string indent(2 * 2 * (MOD_DEPTH - depth) + 1, ' ');
                //     std::cerr << indent << h_smiles << ": " << (reason.empty() ? "Pass!" : reason)
                //             << " (" << j << " / " << h_size << " -> " << j * 100 / h_size << "%)\n";
                // }
                // std::cout << h_smiles << ": " << (*filter)(result) << std::endl;
                // std::string reason;
                // result->getProp("failedfilter", reason);
                result_pool.append_smiles(result, h_smiles, reason);
                // add_h(result);
            }
            // add_h(result);
            recursive_modification(result, h_smiles, except_h, result_pool, depth - 1);
            if (quit.load())
                return r;
        }
    }
    return r;
}

void recursive_modification(RWMOL_SPTR mol, std::string mol_smiles, std::vector<int>& except_h, Result& result_pool, int depth)
{
    if (depth == 0)
        return;
    
    // TODO: reduce memory, move enumerate_h to replace_h_and_filt
    // TODO: reduce runtime, declare h_index_list as shared_ptr, and pass it to replace_h_and_filt
    // (X) enumerate all H still use too much memory
    int max_mod = 0;
    for (auto& item : mod_list)
        if (std::get<1>(item) > max_mod)
            max_mod = std::get<1>(item);
    H_INDEX_LIST_SPTR h_index_list(new TripleIntVec(1));
    h_index_list->reserve(max_mod + 1);
    for (int i = 1; i <= max_mod; ++i)
        h_index_list->emplace_back(enumerate_h(mol, except_h, i));

    unsigned int mod_size = mod_list.size(), remaining;
    // unsigned int mod_print_interval = (unsigned int)(mod_size / 10);
    // std::vector<std::vector<int>> h_list;
    for (unsigned int i = 0; i < mod_size; ++i)
    {
        // if (i % mod_print_interval == 0)
        // {
        //     std::string indent(2 * 2 * (MOD_DEPTH - depth), ' ');
        //     std::cerr << indent << "Mol = " << mol_smiles << ", Mod = " << std::get<0>(mod_list[i])
        //             << " (" << i << " / " << mod_size << " -> " << i * 100 / mod_size << "%)\n";
        // }
        // std::cerr << "Add Modification: " << std::get<0>(mod_list[i]) << " ...\n\n\n";
        thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, std::ref(except_h), std::ref(result_pool), i, depth));
        remaining = thread_pool.get_task_num();
        if (depth == MOD_DEPTH && remaining > 100)
        {
            std::this_thread::sleep_for(std::chrono::seconds(10));
            remaining = thread_pool.get_task_num();
            std::cerr << "Task remain: " << remaining << "\n";
        }
    }
}

void got_signal(int signal)
{
    std::cerr << "Signal raised: " << signal << "\n";
    quit.store(true);
}

int main(int argc, char **argv)
{
    // RDLogger.DisableLog('rdApp.*')
    std::signal(SIGINT, got_signal);

    const std::string target = argv[1];
    std::string h_str = argv[argc - 1];
    RWMOL_SPTR mol(SmilesToMol(target));
    MolOps::addHs(*mol);
    print_molecule(mol);
    Result result;
    // for (int i = 2; i < argc - 1; ++i)
    //     mod_list.emplace_back(mod_func(argv[i]))
    std::vector<int> except_h;
    for (int i = 0; i < mol->getNumAtoms(); ++i)
        except_h.emplace_back(i);
    size_t pos = 0;
    std::string delimiter(",");
    while ((pos = h_str.find(delimiter)) != std::string::npos)
    {
        except_h.erase(except_h.begin() + std::stoi(h_str.substr(0, pos)));
        h_str.erase(0, pos + delimiter.length());
    }
    except_h.erase(except_h.begin() + std::stoi(h_str));
    // ProfilerStart("filter.prof");
    recursive_modification(mol, target, except_h, result, MOD_DEPTH);
    // MolOps::removeHs(*mol);
    // predict_spectrum(mol, result, suppress_exceptions);

    std::this_thread::sleep_for(std::chrono::seconds(5));
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
    // ProfilerStop();
    std::cerr << "\nTotal " << result.get_candidate_size() << " candidate\n";
    std::cerr << "Total " << result.get_filtered_size() << " filtered\n";
    return 0;
}