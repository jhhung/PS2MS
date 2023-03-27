#include <iostream>
#include <unordered_set>
#include <map>
#include <fstream>
#include <iomanip>
#include <csignal>
#include <ctime>
#include <cmath>
#include <thread>

#include "thread_pool.hpp"
#include "DruglikeFilter.hpp"
#include "modifications.hpp"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SanitException.h>
using namespace RDKit; 

#include <filesystem>
#include <boost/algorithm/string.hpp>

std::vector<std::string> exceptions;
constexpr int MOD_DEPTH = 1;
constexpr int WRITE_FREQ = 400000000;
std::atomic<bool> quit(false);
ThreadPool<std::pair<int, std::string>> thread_pool;
time_t start_time;
int mod_start_idx, mod_end_idx;
std::ofstream logging;
std::string PATH_PREFIX = "/tmp/";

using FILTER_SPTR = boost::shared_ptr<Filter>;
using TripleIntVec = std::vector<std::vector<std::vector<uint8_t>>>;
using H_INDEX_LIST_SPTR = boost::shared_ptr<TripleIntVec>;

class ReasonFilterCount
{
  public:
    ReasonFilterCount(const std::string &num)
    {
      record_filename = PATH_PREFIX + "reason_count_record_" + num + ".txt";
      if (std::filesystem::exists(record_filename))
          std::cerr << "[ERROR] Reason count file still exists.\n";
      std::filesystem::remove(record_filename);
    }

    ~ReasonFilterCount()
    {
      writeFile();
    }

    void addReasonCount(const std::string& reason)
    {
      if (reason_counter.count(reason))
        ++reason_counter[reason];
      else
        reason_counter[reason] = 1;
    }

    void writeFile()
    {
      std::ofstream outfile(record_filename, std::ios::app);
      for (auto iter = reason_counter.begin(); iter!=reason_counter.end(); ++iter)
      {
        outfile << iter->first << "\t" << iter->second << "\n";
      }
      outfile.close();
      reason_counter.clear();
    }

  private:
    std::unordered_map<std::string, unsigned long long> reason_counter;
    std::string record_filename;
};

class FunctionalGroupFilterCount
{
  public:
    FunctionalGroupFilterCount(const std::string &num)
    {
      record_filename = PATH_PREFIX + "fg_count_record_" + num + ".txt";
      if (std::filesystem::exists(record_filename))
          std::cerr << "[ERROR] Functuinal group count file still exists.\n";
      std::filesystem::remove(record_filename);
    }

    ~FunctionalGroupFilterCount()
    {
      writeFile();
    }

    void addFGCount(const std::string& functional_group)
    {
      if (fg_counter.count(functional_group))
        ++fg_counter[functional_group];
      else
        fg_counter[functional_group] = 1;
    }

    void writeFile()
    {
      std::ofstream outfile(record_filename, std::ios::app);
      for (auto iter = fg_counter.begin(); iter!=fg_counter.end(); ++iter)
      {
        outfile << iter->first << "\t" << iter->second << "\n";
      }
      outfile.close();
      fg_counter.clear();
    }

  private:
    std::unordered_map<std::string, unsigned long long> fg_counter;
    std::string record_filename;
};

struct Result
{
  public:
    Result(std::string num) : candidate_size(0), filtered_size(0)
    {
        candidate_filename = PATH_PREFIX + "candidate_smiles_" + num + ".txt";
        filtered_filename = PATH_PREFIX + "filtered_smiles_" + num + ".txt";
        if (std::filesystem::exists(candidate_filename)
            || std::filesystem::exists(filtered_filename))
        {
            std::cerr << "[ERROR] Result file still exists.\n";
            // exit(-1); // modified
        }
        std::filesystem::remove(candidate_filename);
        std::filesystem::remove(filtered_filename);
    }

    ~Result()
    {
        write_file();
    }

    bool find(std::string& key)
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return candidate.find(key) != candidate.end()
                || filtered.find(key) != filtered.end();
    }

    std::uint64_t get_candidate_size()
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return candidate_size + candidate.size();
    }

    std::uint64_t get_filtered_size()
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return filtered_size + filtered.size();
    }

    std::uint64_t get_total_size()
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return candidate_size + filtered_size + candidate.size() + filtered.size();
    }

    void add_filtered(const std::string& h_smiles, const std::string& reason)
    {
        std::unique_lock<std::shared_mutex> lock(mutex);
        filtered.emplace(h_smiles, reason);
        if ((candidate_size + filtered_size) % WRITE_FREQ == 0)
            write_file();
    }

    void add_candidate(RWMOL_SPTR mol, const std::string& h_smiles)
    {
        std::unique_lock<std::shared_mutex> lock(mutex);
        candidate.emplace(h_smiles);
        if ((candidate_size + filtered_size) % WRITE_FREQ == 0)
            write_file();
    }

    void write_file()
    {
        std::ofstream outfile(candidate_filename, std::ios::app);
        for (const auto& smi : candidate)
            outfile << smi << "\n";
        outfile.close();
        // outfile.open(filtered_filename, std::ios::app);
        // for (auto& [key, value] : filtered)
        //     // outfile << std::setw(150) << std::left << key << "(" << value << ")\n";
        //     outfile << key << " " << value << "\n";
        // outfile.close();
        candidate_size += candidate.size();
        filtered_size += filtered.size();
        candidate.clear();
        filtered.clear();
    }

  private:
    std::uint64_t candidate_size;
    std::uint64_t filtered_size;
    std::string candidate_filename;
    std::string filtered_filename;
    mutable std::shared_mutex mutex;
    std::unordered_set<std::string> candidate;
    std::unordered_map<std::string, std::string> filtered;
};

void dfs(RWMOL_SPTR mol, uint8_t start_id, std::set<uint8_t>& remaining)
{
    remaining.erase(start_id);
    if (remaining.empty())
        return;
    for (const auto &nbri : boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(start_id))))
    {
        uint8_t nbr_id = (*mol)[nbri]->getIdx();
        if (remaining.find(nbr_id) != remaining.end())
            dfs(mol, nbr_id, remaining);
        if (remaining.empty())
            return;
    }
}

bool connected_component(RWMOL_SPTR mol, std::vector<uint8_t>& h_list)
{
    switch (h_list.size())
    {
        case 1:
            return true;
        case 2:
        {
            int target1 = get_target(mol, h_list[0]);
            int target2 = get_target(mol, h_list[1]);
            return target1 == target2 || mol->getBondBetweenAtoms(target1, target2) != nullptr;
        }
        default:
        {
            std::set<uint8_t> remaining;
            for (const auto& item : h_list)
                remaining.emplace(get_target(mol, item));
            dfs(mol, *remaining.begin(), remaining);
            if (remaining.empty())
                return true;
            else
                return false;
        }
    }
}

std::vector<std::vector<uint8_t>> enumerate_h(RWMOL_SPTR mol, int h_num)
{
    std::vector<Atom*> all_h;
    all_h.reserve(mol->getNumAtoms());
    for (auto& atom : mol->atoms())
        if (atom->getSymbol() == "H")
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
    
    // map index back to atom_index and filter too far H
    std::vector<std::vector<uint8_t>> index_list_rm_far;
    for (auto &item : index_list)
    {   
        for (std::size_t i = 0; i < item.size(); ++i)
            item[i] = all_h[item[i] - 1]->getIdx();
        if (connected_component(mol, item))
            index_list_rm_far.emplace_back(item);
    }
    return index_list_rm_far;
}

//brian
std::vector<ModificationType> mod_list = get_modification_list();
//std::vector<ModificationType> mod_list = get_modification_test_list();
std::vector<ModificationType> simple_mod_list = get_simple_modification_list();
// std::vector<ModificationType> mod_list;

void add_h(RWMOL_SPTR mol)
{
    if (mol->needsUpdatePropertyCache())
        mol->updatePropertyCache();
    MolOps::addHs(*mol);
}

void remove_h(RWMOL_SPTR mol)
{
    if (mol->needsUpdatePropertyCache())
        mol->updatePropertyCache();
    MolOps::removeHs(*mol);
}

RWMOL_SPTR do_modification(
    std::vector<uint8_t>& h_list, RWMOL_SPTR parent, const std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>& mod)
{
    RWMOL_SPTR result = boost::make_shared<RWMol>(*parent, true);
    mod(result, h_list);
    // if (result->needsUpdatePropertyCache())
    //     result->updatePropertyCache();
    // add_h(result);
    // MolOps::removeHs(*result);
    // add_h(result);
    return result;
}

void recursive_modification(RWMOL_SPTR mol, std::string mol_smiles, FILTER_SPTR filter, Result& result_pool, ReasonFilterCount& reason_count, FunctionalGroupFilterCount &fg_count, int depth, int suppress_exceptions);
std::mutex print_mutex;
std::pair<int, std::string> replace_h_and_filt(RWMOL_SPTR mol, H_INDEX_LIST_SPTR h_index_list, FILTER_SPTR filter, Result& result_pool, ReasonFilterCount& reason_count, FunctionalGroupFilterCount& fg_count, ModificationType& mod, int depth, int suppress_exceptions)
{
    // ModificationType& mod = mod_list[mod_num];
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
                // add_h(result);
                remove_h(result);
            }
            catch (...)
            {
                // std::lock_guard<std::mutex> lock(print_mutex);
                // std::cerr << "//////////////////////////////////////////////////" << std::endl;
                // std::cerr << depth << " " << std::get<0>(mod) << " ";
                // for (auto& item : h_list[j])
                //     std::cerr << (int)item << " ";
                // std::cerr << std::endl;
                // print_molecule(result);
                // std::cerr << "//////////////////////////////////////////////////" << std::endl;
                continue;
            }
            // std::string h_smiles = MolToSmiles(*result);
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
                reason = (*filter)(result);
                // if (j % h_print_interval == 0)
                // {
                //     std::string indent(2 * 2 * (MOD_DEPTH - depth) + 1, ' ');
                //     std::cerr << indent << h_smiles << ": " << (reason.empty() ? "Pass!" : reason)
                //             << " (" << j << " / " << h_size << " -> " << j * 100 / h_size << "%)\n";
                // }
                // std::cout << h_smiles << ": " << (*filter)(result) << std::endl;
                // std::string reason;
                // result->getProp("failedfilter", reason);
                // if (!result_pool.find(h_smiles)) // check again to prevent "some" race condition
                // {
                if (reason.empty())
                {
                    result_pool.add_candidate(result, h_smiles);
                    if (depth != 1)
                    {
                        add_h(result);
                        recursive_modification(result, h_smiles, filter, result_pool, reason_count, fg_count, depth - 1, suppress_exceptions);
                    }
                }
                else
                {
                    auto fg = std::get<0>(mod);
                    result_pool.add_filtered(h_smiles, reason);
                    reason_count.addReasonCount(reason);
                    fg_count.addFGCount(fg);
                }
            }
            if (quit.load())
                return r;
        }
    }
    return r;
}

void recursive_modification(RWMOL_SPTR mol, std::string mol_smiles, FILTER_SPTR filter, Result& result_pool, ReasonFilterCount& reason_count, FunctionalGroupFilterCount& fg_count, int depth, int suppress_exceptions)
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
    
    H_INDEX_LIST_SPTR h_index_list = boost::make_shared<TripleIntVec>(1);
    h_index_list->reserve(max_mod + 1);
    for (int i = 1; i <= max_mod; ++i)
        h_index_list->emplace_back(enumerate_h(mol, i));

    unsigned int mod_size = mod_list.size(), profile_count = 0, remaining;
    //unsigned int mod_size = simple_mod_list.size(), profile_count = 0, remaining;
    // unsigned int mod_print_interval = (unsigned int)(mod_size / 10);
    // std::vector<std::vector<int>> h_list;
    if (depth == MOD_DEPTH)
    {
        for (unsigned int i = mod_start_idx; i < mod_end_idx; ++i)
        {
            thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(reason_count), std::ref(fg_count), std::ref(mod_list[i]), depth, suppress_exceptions));
            //thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(simple_mod_list[i]), depth, suppress_exceptions));
            std::this_thread::sleep_for(std::chrono::seconds(1));
            while ((remaining = thread_pool.get_task_num()) > 1.0e5 && !quit.load())
            {
                logging << "Stop submit task at mod " << i + 2 << "/" << mod_size << ", Task remain: " << remaining << ", Candidates: " << result_pool.get_candidate_size() << ", Filtered: " << result_pool.get_filtered_size() << ", Speed: " << result_pool.get_total_size() / difftime(time(NULL), start_time) << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(60));
            }
        }
    }
    else
    {/*
	if (depth == MOD_DEPTH - 1)
        	for (unsigned int i = 0; i < mod_list.size(); ++i)
            		thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(mod_list[i]), depth, suppress_exceptions));
        else
		for (unsigned int i = 0; i < simple_mod_list.size(); ++i)
            		thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(simple_mod_list[i]), depth, suppress_exceptions));
        */
        for (unsigned int i = 0; i < mod_list.size(); ++i)
        	thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(reason_count), std::ref(fg_count), std::ref(mod_list[i]), depth, suppress_exceptions));
        /*
	for (unsigned int i = 0; i < simple_mod_list.size(); ++i)
        	thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(simple_mod_list[i]), depth, suppress_exceptions));
	*/
    }
}

void got_signal(int signal)
{
    logging << "Signal raised: " << signal << "\n";
    quit.store(true);
}

int main(int argc, char **argv)
{
    // RDLogger.DisableLog('rdApp.*')
    std::signal(SIGINT, got_signal);
    start_time = time(NULL);
    if (argc > 3)
        // PATH_PREFIX += argv[3];
        PATH_PREFIX = argv[3];
    else {
        PATH_PREFIX += "output_amph_depth1/";
        //PATH_PREFIX += "output_cath/";
        //PATH_PREFIX += "output_non_poison_depth_2/";
    }
    std::filesystem::create_directories(PATH_PREFIX);
    logging.open(PATH_PREFIX + argv[2] + ".log");

    //Initialise model configuration
    int suppress_exceptions = 1;
    // std::string param_filename = "param_output.log";
	// std::string config_filename = "param_config.txt";
	// if( !std::filesystem::exists( config_filename ) ){
	// 	std::cout << "Could not find file: " <<  config_filename << std::endl;
	// 	throw FileException("Could not find file: " + config_filename);
	// }
	// initConfig( cfg, config_filename );

	//Read in the parameters
	// if( !std::filesystem::exists( param_filename ) ){
	// 	std::cout << "Could not find file: " <<  param_filename << std::endl;
	// 	throw FileException("Could not find file: " + param_filename);
	// }
	// if( cfg.theta_function == NEURAL_NET_THETA_FUNCTION )
	// 	nn_param = new NNParam( param_filename );
	// else
	// 	param = new Param(param_filename);

    /*
    std::vector<ModificationType> constraint_mod;
    for (auto& item : simple_mod_list)
        if (std::get<1>(item) <= 2)
            constraint_mod.emplace_back(item);
    simple_mod_list = constraint_mod;
    logging << "# of Functional Group: " << simple_mod_list.size() << std::endl;
    
    int stride = std::ceil(simple_mod_list.size() / std::stoi(argv[1]));
    mod_start_idx = std::stoi(argv[2]) * stride;
    mod_end_idx = std::stoi(argv[2]) == std::stoi(argv[1]) - 1 ? simple_mod_list.size() : (std::stoi(argv[2]) + 1) * stride;
    logging << "Job: " << argv[2] << " (" << mod_start_idx << " - " << mod_end_idx << ")\n";
    */

    std::vector<ModificationType> constraint_mod;
    for (auto& item : mod_list)
        if (std::get<1>(item) <= 2)
            constraint_mod.emplace_back(item);
    mod_list = constraint_mod;
    logging << "# of Functional Group: " << mod_list.size() << std::endl;
    
    // for (auto& item : mod_list2)
    //     if (std::get<1>(item) <= 3)
    //         mod_list.emplace_back(item);

    int stride = std::ceil(mod_list.size() / std::stoi(argv[1]));
    mod_start_idx = std::stoi(argv[2]) * stride;
    mod_end_idx = std::stoi(argv[2]) == std::stoi(argv[1]) - 1 ? mod_list.size() : (std::stoi(argv[2]) + 1) * stride;
    logging << "Job: " << argv[2] << " (" << mod_start_idx << " - " << mod_end_idx << ")\n";

    const std::string amphetamine = "NC(CC1=CC=CC=C1)C";
    //const std::string cathinone = "O=C(c1ccccc1)[C@@H](N)C";
    //const std::string non_poison_mol = "CCCCN(CCCC)C(=O)Nc1ccccc1";
    RWMOL_SPTR mol(SmilesToMol(amphetamine));
    //RWMOL_SPTR mol(SmilesToMol(cathinone));
    //RWMOL_SPTR mol(SmilesToMol(non_poison_mol));
    MolOps::addHs(*mol);
    print_molecule(mol);
    FILTER_SPTR filter = boost::make_shared<DruglikeFilter>();
    Result result(argv[2]);
    ReasonFilterCount reason_count(argv[2]);
    FunctionalGroupFilterCount fg_count(argv[2]);
    // ProfilerStart("build_database.prof");
    // HeapProfilerStart("memory_profile/build_database.hprof");
    recursive_modification(mol, amphetamine, filter, result, reason_count, fg_count, MOD_DEPTH, suppress_exceptions);
    //recursive_modification(mol, cathinone, filter, result, MOD_DEPTH, suppress_exceptions);
    //recursive_modification(mol, non_poison_mol, filter, result, MOD_DEPTH, suppress_exceptions);
    // MolOps::removeHs(*mol);
    // predict_spectrum(mol, result, suppress_exceptions);
    std::this_thread::sleep_for(std::chrono::seconds(5));
    unsigned int remain = thread_pool.get_task_num(), profile_times = 0;
    logging << "Task add finished, size = " << remain << std::endl;
    while (remain != 0 && !quit.load())
    {
        std::this_thread::sleep_for(std::chrono::seconds(60));
        remain = thread_pool.get_task_num();
        logging << "Task remain: " << remain << ", Candidates: " << result.get_candidate_size() << ", Filtered: " << result.get_filtered_size() << ", Speed: " << result.get_total_size() / difftime(time(NULL), start_time) << std::endl;
        // if (++profile_times % 4320 == 0)
        //     HeapProfilerDump("12 hours");
    }
    logging << "Wait threads to complete remaining tasks ...\n";
    thread_pool.terminate_all_thread();
    // ProfilerStop();
    // HeapProfilerStop();
    result.write_file();
    reason_count.writeFile();
    fg_count.writeFile();
    logging << "Write completed\n";
    logging << "\nTotal " << result.get_candidate_size() << " candidate\n";
    logging << "Total " << result.get_filtered_size() << " filtered\n";
    logging << "Average iteration per second: " << result.get_total_size() / difftime(time(NULL), start_time) << "\n";
    return 0;
}
