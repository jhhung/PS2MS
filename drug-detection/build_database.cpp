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
constexpr int MOD_DEPTH = 4;
constexpr int WRITE_FREQ = 400000000;
std::atomic<bool> quit(false);
ThreadPool<std::pair<int, std::string>> thread_pool;
time_t start_time;
int mod_start_idx, mod_end_idx;
std::ofstream logging;
std::string PATH_PREFIX = "/mnt/ec/mammoth/blender/drug/data/";

using FILTER_SPTR = boost::shared_ptr<Filter>;
using TripleIntVec = std::vector<std::vector<std::vector<uint8_t>>>;
using H_INDEX_LIST_SPTR = boost::shared_ptr<TripleIntVec>;

class ReasonFilterCount
{
  public:
    ReasonFilterCount(const std::string &num)
    {
      record_filename = PATH_PREFIX + "reason_count_record_" + num + ".txt";
    //   if (std::filesystem::exists(record_filename))
    //       std::cerr << "[ERROR] Reason count file still exists.\n";
    //   std::filesystem::remove(record_filename);
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

        // ascending permution of r-elements subset in range of 1~h_num
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

std::vector<std::vector<ModificationType>> mod_list = {
    // depth 3.5
    {
        make_alkyl(0, 1),
        make_alkyl(1, 1),
        make_alkyl(3, 1),
    },
    // depth 3
    {
        make_dialkyl(0),
        make_dialkyl(1, 1),
        make_dialkyl(3, 1),
        make_acetoxy(),
        make_acetyl(),
        make_acryloyl(),
        make_acyloin(1),
        make_acylsilane(3),
        make_alcohol(),
        make_aldehyde(),
        make_alkene(1),
        make_alkyl(0, 1),
        make_alkyl(1),
        make_alkyl(2),
        make_alkyl(3),
        make_alkyl(4),
        make_alkyl(5),
        make_alkyl(1, 1),
        make_alkyl(2, 1),
        make_alkyl(3, 1),
        make_allyl(),
        make_amide(2),
        make_amine(2),
        make_aminophosphine_r2(0,2),
        make_aminophosphine_r2(1,1),
        make_aminophosphine_r2(1,2),
        make_aminophosphine_r3(2),
        make_basic_aluminium_r1(),
        make_basic_aluminium_r2(1),
        make_benzene(),
        make_benzyl(),
        make_beta_lactone("O", 3),
        make_beta_lactone("S", 3),
        make_boronic_acid(),
        make_carbene(1),
        make_carbonate_ester(1),
        make_carbonyl(1),
        make_carboximidate(),
        make_cyanohydrin(1),
        make_cyanomethyl(),
        make_cyclopropyl(),
        make_delta_lactone("O", 6),
        make_dioxolane(),
        make_dithioacetal(1),
        make_enol(1),
        make_enol_ether(3),
        make_enone(3),
        make_episulfide(3),
        make_epoxide(3),
        make_ester(1),
        make_ether(1),
        make_gamma_lactone("O",4),
        make_haloalkane("Br"),
        make_haloalkane("Cl"),
        make_haloalkane("F"),
        make_haloalkane("I"),
        make_hemithioacetal(1),
        make_iminium(3),
        make_ketone(1),
        make_ketyl(2),
        make_sulfide_ether(),
        make_sulfide_ether(1),
        make_methanedithiol(2),
        make_methylidene(),
        make_methylmorpholine(),
        make_nitrene(),
        make_nitronate(3),
        make_ortho_quinone_methide(),
        make_para_quinone_methide(),
        make_phosphaalkene(2),
        make_phosphinate(1),
        make_phosphinate(2),
        make_phosphine(2),
        make_phosphine_oxide(2),
        make_phosphinite(1),
        make_phosphinite(2),
        make_phosphinous_acids(1),
        make_phosphonite(1),
        make_phosphonite(2),
        make_phosphonium(3),
        make_phosphorane(4),
        make_propenyl(),
        make_selenenic_acid(),
        make_selenol(),
        make_silyl_ether(2),
        make_silyl_ether(3),
        make_tellurol(),
        make_test4(0),
        make_thioacetal(1),
        make_thioacetal(2),
        make_thioester(1),
        make_thioketal(2),
        make_tosyl(),
        make_vinyl(),
        make_vinylene(1),
        make_vinylidene(3),
        make_ynone(1)
    },
    // depth 2
    {
        make_thiocarbamate(false, 2),
        make_acetoxy(),
        make_acetyl(),
        make_acryloyl(),
        make_acylal(1),
        make_acylal(2),
        make_acyloin(1),
        make_acylsilane(2),
        make_acylsilane(3),
        make_alcohol(),
        make_aldehyde(),
        make_alkene(1),
        make_alkyl(0, 1),
        make_alkyl(1),
        make_alkyl(2),
        make_alkyl(3),
        make_alkyl(4),
        make_alkyl(5),
        make_alkyl(1, 1),
        make_alkyl(2, 1),
        make_alkyl(3, 1),
        make_alkyl(4, 1),
        make_alkyne(1),
        make_allyl(),
        make_amide(2),
        make_aminal(2),
        make_amine(2),
        make_aminophosphine_r2(0,2),
        make_aminophosphine_r2(1,1),
        make_aminophosphine_r2(1,2),
        make_aminophosphine_r3(1),
        make_aminophosphine_r3(2),
        make_azine(2),
        make_basic_aluminium_r1(),
        make_basic_aluminium_r2(),
        make_basic_aluminium_r2(1),
        make_benzene(),
        make_benzyl(1),
        make_beta_lactone("O", 3),
        make_beta_lactone("S", 3),
        make_boronic_acid(),
        make_carbamate(2),
        make_carbene(0),
        make_carbene(1),
        make_carbonate_ester(1),
        make_carbonyl(1),
        make_carboximidate(),
        make_cyanohydrin(1),
        make_cyanomethyl(),
        make_cyclopropyl(),
        make_delta_lactone("O", 6),
        make_dioxolane(),
        make_dithioacetal(1),
        make_enol(1),
        make_enol_ether(3),
        make_enone(2),
        make_enone(3),
        make_episulfide(3),
        make_epoxide(3),
        make_ester(),
        make_ester(1),
        make_ether(1),
        make_gamma_lactone("O",4),
        make_geminal_diol(),
        make_haloalkane("Br"),
        make_haloalkane("Cl"),
        make_haloalkane("F"),
        make_haloalkane("I"),
        make_hemithioacetal(0),
        make_hemithioacetal(1),
        make_hydrazone(1),
        make_hydroxamic_acid(1),
        make_imide(2),
        make_iminium(3),
        make_ketone(1),
        make_ketyl(1),
        make_ketyl(2),
        make_sulfide_ether(),
        make_sulfide_ether(1),
        make_methanedithiol(2),
        make_methylene_bridge(),
        make_methylidene(),
        make_methylmorpholine(),
        make_nitrene(),
        make_nitro(),
        make_nitronate(2),
        make_nitronate(3),
        make_ortho_quinone_methide(),
        make_para_quinone_methide(),
        make_phosphaalkene(2),
        make_phosphaalkyne(),
        make_phosphinate(1),
        make_phosphinate(2),
        make_phosphine(1),
        make_phosphine(2),
        make_phosphine_oxide(2),
        make_phosphinite(1),
        make_phosphinite(2),
        make_phosphinous_acids(1),
        make_phosphite_ester(1),
        make_phosphite_ester(2),
        make_phosphonite(1),
        make_phosphonite(2),
        make_phosphonium(3),
        make_phosphoramidate(3),
        make_phosphoramidite(3),
        make_phosphorane(4),
        make_phosphoryl(),
        make_propenyl(),
        make_selenenic_acid(),
        make_selenol(),
        make_selenonic_acid(),
        make_silyl_ether(2),
        make_silyl_ether(3),
        make_sulfonamide(2),
        make_sulfonanilide(1),
        make_tellurol(),
        make_test4(0),
        make_thioacetal(1),
        make_thioacetal(2),
        make_thioamide(2),
        make_thioester(0),
        make_thioester(1),
        make_thioketal(2),
        make_thiourea(true,3),
        make_tosyl(),
        make_urea(true,3),
        make_vinyl(),
        make_vinylene(1),
        make_vinylidene(3),
        make_ynone(1)
    },
    // depth 1
    {
        make_alkyl(0, 1),
        make_alkyl(1),
        make_alkyl(2),
        make_alkyl(3),
        make_alkyl(4),
        make_alkyl(5),
        make_alkyl(1, 1),
        make_alkyl(2, 1),
        make_alkyl(3, 1),
        make_alkyl(4, 1),
        make_acetal(),
        make_acetal(1),
        make_acetal(2),
        make_acetal(3),
        make_acetoxy(),
        make_acetyl(),
        make_acetylide("Li"),
        make_acetylide("Na"),
        make_acetylide("K"),
        make_acetylide("Mg"),
        make_acetylide("Ca"),
        make_acetylide("Fe"),
        make_acetylide("Al"),
        make_acid_anhydride(),
        make_acid_anhydride(1),
        make_acryloyl(),
        make_acyl_azide(),
        make_acyl_chloride(),
        make_acyl_halide("F"),
        make_acyl_halide("Cl"),
        make_acyl_halide("Br"),
        make_acyl_halide("I"),
        make_acylal(),
        make_acylal(1),
        make_acylal(2),
        make_acylhydrazine(),
        make_acylhydrazine(1),
        make_acylhydrazine(2),
        make_acylhydrazine(3),
        make_acyloin(),
        make_acyloin(1),
        make_acylsilane(),
        make_acylsilane(1),
        make_acylsilane(2),
        make_acylsilane(3),
        make_acylurea(),
        make_alcohol(),
        make_aldehyde(),
        make_aldimine(),
        make_aldimine(1),
        make_alkene(),
        make_alkene(1),
        make_alkoxide(),
        make_alkyl_nitrites(),
        make_alkyne(),
        make_alkyne(1),
        make_allyl(),
        make_amide(),
        make_amide(1),
        make_amide(2),
        make_amidine(),
        make_amidrazone(),
        make_amidrazone(true),
        make_aminal(),
        make_aminal(1),
        make_aminal(2),
        make_aminal(3),
        make_amine(),
        make_amine(1),
        make_amine(2),
        make_amine_oxide(),
        make_amine_oxide(1),
        make_amine_oxide(2),
        make_aminophosphine_r3(),
        make_aminophosphine_r3(1),
        make_aminophosphine_r3(2),
        make_aminophosphine_r2(0, 0),
        make_aminophosphine_r2(0, 1),
        make_aminophosphine_r2(0, 2),
        make_aminophosphine_r2(1, 0),
        make_aminophosphine_r2(1, 1),
        make_aminophosphine_r2(1, 2),
        make_aminophosphine_r2(2, 0),
        make_aminophosphine_r2(2, 1),
        make_aminophosphine_r1(0, 0),
        make_aminophosphine_r1(0, 1),
        make_aminophosphine_r1(0, 2),
        make_aminophosphine_r1(0, 3),
        make_aminophosphine_r1(0, 4),
        make_aminophosphine_r1(1, 0),
        make_aminophosphine_r1(1, 1),
        make_aminophosphine_r1(1, 2),
        make_aminophosphine_r1(1, 3),
        make_aminophosphine_r0(0),
        make_aminophosphine_r0(1),
        make_aminophosphine_r0(2),
        make_aminophosphine_r0(3),
        make_aminophosphine_r0(4),
        make_aminophosphine_r0(5),
        make_aminoxyl(),
        make_aminoxyl(1),
        make_azide(),
        make_azide(false),
        make_azine(),
        make_azine(1),
        make_azine(2),
        make_azine(3),
        make_aziridine(),
        make_aziridine(1),
        make_aziridine(2),
        make_aziridine(3),
        make_aziridine(4),
        make_azo(),
        make_azo(1),
        make_azole_3_hands_1_atom("N", 2),
        make_azole_3_hands_2_atoms("N"),
        make_azole_3_hands_2_atoms("N", 1),
        make_azole_2_hands_1_atom("O", 1),
        make_azoxy(),
        make_azoxy(1),
        make_basic_aluminium_r2(),
        make_basic_aluminium_r2(1),
        make_basic_aluminium_r1(),
        make_benzylidene_acetal(),
        make_benzylidene_acetal(1),
        make_bisthiosemicarbazone(),
        make_bisthiosemicarbazone(1),
        make_bisthiosemicarbazone(2),
        make_bisthiosemicarbazone(3),
        make_biuret(),
        make_biuret(1),
        make_biuret(2),
        make_boronic_acid(),
        make_carbamate(),
        make_carbamate(1),
        make_carbamate(2),
        make_carbamoyl_chloride(),
        make_carbamoyl_chloride(1),
        make_carbazide(),
        make_carbazide(1),
        make_carbene(),
        make_carbene(1),
        make_carbodiimide(),
        make_carbodiimide(1),
        make_carbonate_ester(),
        make_carbonate_ester(1),
        make_carbonyl(),
        make_carbonyl(1),
        make_carboximidate(),
        make_carboximidate(1),
        make_carboximidate(2),
        make_carboxylic_acid(),
        make_chloroformate(),
        make_cyanate(),
        make_cyanate(false),
        make_cyanate_ester(),
        make_cyanimide(),
        make_cyanimide(1),
        make_cyanohydrin(),
        make_cyanohydrin(1),
        make_cyanomethyl(),
        make_cyclopropyl(),
        make_diazo(),
        make_diazo(1),
        make_diazo(0, false),
        make_diazo(1, false),
        make_diazonium("F"),
        make_diazonium("Cl"),
        make_diazonium("Br"),
        make_diazonium("I"),
        make_dicarbonate(),
        make_dicarbonate(1),
        make_diketopiperazine(),
        make_diketopiperazine(1),
        make_diketopiperazine(0, 1),
        make_diketopiperazine(1, 1),
        make_diketopiperazine(0, 2),
        make_diketopiperazine(1, 2),
        make_dioxazolone(),
        make_dioxirane(),
        make_dioxirane(1),
        make_diphenyltriazene(),
        make_diphenyltriazene(1),
        make_disulfide(),
        make_disulfide(1),
        make_dithiocarbamate(),
        make_dithiocarbamate(1),
        make_dithiocarbamate(2),
        make_dithiol(),
        make_enamine(),
        make_enamine(1),
        make_enamine(2),
        make_enamine(3),
        make_enamine(4),
        make_vicinal_diol(),
        make_geminal_diol(),
        make_enediyne(),
        make_enediyne(1),
        make_enediyne(2),
        make_enediyne(3),
        make_enol(),
        make_enol(1),
        make_enol_ether(),
        make_enol_ether(1),
        make_enol_ether(2),
        make_enol_ether(3),
        make_enone(),
        make_enone(1),
        make_enone(2),
        make_enone(3),
        make_episulfide(),
        make_episulfide(1),
        make_episulfide(2),
        make_episulfide(3),
        make_epoxide(),
        make_epoxide(1),
        make_epoxide(2),
        make_epoxide(3),
        make_ester(),
        make_ester(1),
        make_ether(),
        make_ether(1),
        make_fluorosulfonate(),
        make_haloalkane("F"),
        make_haloalkane("Cl"),
        make_haloalkane("Br"),
        make_haloalkane("I"),
        make_halohydrin("F"),
        make_halohydrin("Cl"),
        make_halohydrin("Br"),
        make_halohydrin("I"),
        make_halohydrin("F", 1),
        make_halohydrin("Cl", 1),
        make_halohydrin("Br", 1),
        make_halohydrin("I", 1),
        make_haloketone("F"),
        make_haloketone("Cl"),
        make_haloketone("Br"),
        make_haloketone("I"),
        make_haloketone("F", 1),
        make_haloketone("Cl", 1),
        make_haloketone("Br", 1),
        make_haloketone("I", 1),
        make_haloketone("F", 2),
        make_haloketone("Cl", 2),
        make_haloketone("Br", 2),
        make_haloketone("I", 2),
        make_hemithioacetal(),
        make_hemithioacetal(1),
        make_carbohydrazides(),
        make_carbohydrazides(1),
        make_carbohydrazides(2),
        make_carbohydrazides(3),
        make_sulfonohydrazides(),
        make_sulfonohydrazides(1),
        make_sulfonohydrazides(2),
        make_sulfonohydrazides(3),
        make_phosphonic_dihydrazides(),
        make_phosphonic_dihydrazides(1),
        make_phosphonic_dihydrazides(2),
        make_phosphonic_dihydrazides(3),
        make_phosphonic_dihydrazides(4),
        make_phosphonic_dihydrazides(5),
        make_phosphonic_dihydrazides(6),
        make_hydrazone(),
        make_hydrazone(1),
        make_hydroperoxide(),
        make_hydroxamic_acid(),
        make_hydroxamic_acid(1),
        make_hydroxylamine(),
        make_hydroxylamine(1),
        make_imide(),
        make_imide(1),
        make_imide(2),
        make_imidic_acid(),
        make_imidic_acid(1),
        make_imidoyl_chloride(),
        make_imidoyl_chloride(1),
        make_imine(),
        make_imine(1),
        make_imine(2),
        make_iminium(),
        make_iminium(1),
        make_iminium(2),
        make_iminium(3),
        make_isocyanate(),
        make_isocyanide(),
        make_isodiazene(),
        make_isodiazene(1),
        make_isodiazomethane(),
        make_isodiazomethane(1),
        make_isothiocyanate(),
        make_isothiouronium(),
        make_ketene(),
        make_ketene(1),
        make_ketenimine(),
        make_ketenimine(1),
        make_ketenimine(2),
        make_ketone(),
        make_ketone(1),
        make_ketyl(),
        make_ketyl(1),
        make_ketyl(2),
        make_beta_lactone("O", 3),
        make_gamma_lactone("O", 4),
        make_delta_lactone("O", 6),
        make_beta_lactone("N"),
        make_beta_lactone("N", 1),
        make_beta_lactone("N", 2),
        make_beta_lactone("N", 3),
        make_gamma_lactone("N", 4),
        make_alpha_lactone("S"),
        make_alpha_lactone("S", 1),
        make_beta_lactone("S"),
        make_beta_lactone("S", 3),
        make_gamma_lactone("S", 4),
        make_gamma_lactone("S", 5),
        make_delta_lactone("S", 6),
        make_methanedithiol(),
        make_methanedithiol(1),
        make_methanedithiol(2),
        make_methanedithiol(3),
        make_methine(),
        make_methine(1),
        make_methine(2),
        make_methylene_bridge(),
        make_methylenedioxy(),
        make_methylenedioxy(1),
        make_methylidene(),
        make_nitrate_ester(),
        make_nitrene(),
        make_nitrile_ylide(),
        make_nitrile_ylide(1),
        make_nitrile_ylide(2),
        make_nitrilimine(),
        make_nitrilimine(1),
        make_nitro(),
        make_nitroalkene(),
        make_nitroalkene(1),
        make_nitroalkene(2),
        make_nitroamine(),
        make_nitroamine(1),
        make_nitrolic_acid(),
        make_nitronate(),
        make_nitronate(1),
        make_nitronate(2),
        make_nitronate(3),
        make_nitrone(),
        make_nitrone(1),
        make_nitrone(2),
        make_nitrosamine(),
        make_nitrosamine(1),
        make_nitroso(),
        make_s_nitrosothiol(),
        make_organic_acid_anhydride(),
        make_organic_acid_anhydride(1),
        make_organic_peroxide(),
        make_organic_peroxide(1),
        make_orthoester(),
        make_orthoester(1),
        make_orthoester(2),
        make_orthoester(3),
        make_oxaziridine(),
        make_oxaziridine(1),
        make_oxaziridine(2),
        make_phosphaalkene(),
        make_phosphaalkene(1),
        make_phosphaalkene(2),
        make_phosphaalkyne(),
        make_phosphate(),
        make_phosphate(1),
        make_phosphate(2),
        make_phosphinate(),
        make_phosphinate(1),
        make_phosphinate(2),
        make_phosphine(),
        make_phosphine(1),
        make_phosphine(2),
        make_phosphine_imide(),
        make_phosphine_imide(1),
        make_phosphine_imide(2),
        make_phosphine_imide(3),
        make_phosphine_oxide(),
        make_phosphine_oxide(1),
        make_phosphine_oxide(2),
        make_phosphinite(),
        make_phosphinite(1),
        make_phosphinite(2), 
        make_phosphinous_acids(),
        make_phosphinous_acids(1),
        make_phosphite_ester(),
        make_phosphite_ester(1),
        make_phosphite_ester(2),
        make_phosphonate(),
        make_phosphonate(1),
        make_phosphonate(2),
        make_phosphonite(),
        make_phosphonite(1),
        make_phosphonite(2),
        make_phosphonium(),
        make_phosphonium(1),
        make_phosphonium(2),
        make_phosphonium(3),
        make_phosphorodiamidate(),
        make_phosphorodiamidate(1),
        make_phosphorodiamidate(2),
        make_phosphorodiamidate(3),
        make_phosphorodiamidate(4),
        make_phosphoramidate(),
        make_phosphoramidate(1),
        make_phosphoramidate(2), 
        make_phosphoramidate(3), 
        make_phosphoramides(),
        make_phosphoramides(1),
        make_phosphoramides(2),
        make_phosphoramides(3),
        make_phosphoramides(4),
        make_phosphoramides(5),
        make_phosphoramidite(),
        make_phosphoramidite(1),
        make_phosphoramidite(2),
        make_phosphoramidite(3),
        make_phosphorane(),
        make_phosphorane(1),
        make_phosphorane(2),
        make_phosphorane(3),
        make_phosphorane(4),
        make_phosphorochloridate(),
        make_phosphorochloridate(1),
        make_phosphochloridite(),
        make_phosphodichloridite(),
        make_phosphodichloridite(1),
        make_phosphoryl(),
        make_propenyl(),
        make_para_quinone_methide(),
        make_ortho_quinone_methide(),
        make_reductone(),
        make_reductone(1),
        make_schiff_base(),
        make_schiff_base(1),
        make_schiff_base(2),
        make_selenenic_acid(),
        make_selenol(),
        make_selenonic_acid(),
        make_selone(),
        make_selone(1),
        make_semicarbazide(),
        make_semicarbazide(1),
        make_semicarbazide(2),
        make_semicarbazide(3),
        make_semicarbazide(4),
        make_semicarbazone(),
        make_semicarbazone(1),
        make_semicarbazone(2),
        make_semicarbazone(3),
        make_semicarbazone(4), 
        make_silyl_enol_ether(),
        make_silyl_enol_ether(1),
        make_silyl_enol_ether(2),
        make_silyl_enol_ether(3),
        make_silyl_enol_ether(4),
        make_silyl_enol_ether(5),
        make_silyl_ether(),
        make_silyl_ether(1),
        make_silyl_ether(2),
        make_silyl_ether(3), 
        make_sulfamoyl_fluoride(),
        make_sulfamoyl_fluoride(1),
        make_sulfenamide(),
        make_sulfenamide(1),
        make_sulfenamide(2),
        make_sulfenic_acid(),
        make_sulfenyl_chloride(),
        make_sulfide(),
        make_sulfide(1),
        make_sulfilimine(),
        make_sulfilimine(1),
        make_sulfilimine(2),
        make_sulfinamide(),
        make_sulfinamide(false),
        make_sulfinamide(true, 1),
        make_sulfinamide(false, 1),
        make_sulfinamide(true, 2),
        make_sulfinamide(false, 2),
        make_sulfinic_acid(),
        make_sulfite_ester(),
        make_sulfite_ester(1),
        make_sulfonamide(),
        make_sulfonamide(1),
        make_sulfonamide(2),
        make_sulfonanilide(),
        make_sulfonanilide(1),
        make_sulfonate(),
        make_sulfonate(1),
        make_sulfone(),
        make_sulfone(1),
        make_sulfonic_acid(),
        make_sulfonyl_halide("F"),
        make_sulfonyl_halide("Cl"),
        make_sulfonyl_halide("Br"),
        make_sulfonyl_halide("I"),
        make_sulfoxide(),
        make_sulfoxide(1),
        make_telluroketone(),
        make_telluroketone(1),
        make_tellurol(),
        make_thiadiazoles(),
        make_thiadiazoles(1),
        make_thiadiazoles(0, 1),
        make_thiadiazoles(1, 1),
        make_thiadiazoles(0, 2),
        make_thiadiazoles(1, 2),
        make_thiadiazoles(0, 3),
        make_thiadiazoles(0, 3),
        make_thial(),
        make_thioacetal(),
        make_thioacetal(1),
        make_thioacetal(2),
        make_dithioacetal(),
        make_dithioacetal(1),
        make_dithioacetal(2),
        make_thioacyl_chloride(),
        make_thioamide(),
        make_thioamide(1),
        make_thioamide(2),
        make_thiocarbamate(),
        make_thiocarbamate(false),
        make_thiocarbamate(true, 1),
        make_thiocarbamate(false, 1),
        make_thiocarbamate(true, 2),
        make_thiocarbamate(false, 2),
        make_thiocarboxylic_acid(),
        make_thiocarboxylic_acid(false),
        make_thiocyanate(),
        make_thioester(),
        make_thioester(1),
        make_thioketal(),
        make_thioketal(1),
        make_thioketal(2),
        make_thioketal(3),
        make_thioketene(),
        make_thioketene(1),
        make_thioketone(),
        make_thioketone(1),
        make_thiol(),
        make_thiophosphate(),
        make_thiophosphate(1),
        make_thiophosphate(2),
        make_thiourea(),
        make_thiourea(false),
        make_thiourea(true, 1),
        make_thiourea(false, 1),
        make_thiourea(true, 2),
        make_thiourea(false, 2),
        make_thiourea(true, 3),
        make_thiourea(false, 3),
        make_tosyl(), 
        make_tosylate(),
        make_tosylhydrazone(),
        make_tosylhydrazone(1),
        make_triazenes(),
        make_triazenes(1), 
        make_triazenes(2), 
        make_triuret(), 
        make_triuret(1), 
        make_triuret(2),
        make_triuret(3),
        make_triuret(4),
        make_triuret(5),
        make_urea(),
        make_urea(false),
        make_urea(true, 1),
        make_urea(false, 1),
        make_urea(true, 2),
        make_urea(false, 2),
        make_urea(true, 3),
        make_urea(false, 3),
        make_vanillyl(),
        make_vinyl(),
        make_vinylene(),
        make_vinylene(1),
        make_vinylidene(),
        make_vinylidene(1),
        make_vinylidene(2),
        make_vinylidene(3),
        make_xanthate("K"),
        make_xanthate("Na"),
        make_xanthate_ester(),
        make_xanthate_ester(1),
        make_ynolate(),
        make_ynolate(1),
        make_ynolate(2),
        make_ynolate(3),
        make_ynone(),
        make_ynone(1),    
        make_benzyl(1),
        make_benzene(0),
        make_thiol(),
        make_thiol(1),
        make_sulfide_ether(),
        make_sulfide_ether(1),
        make_methylmorpholine(),
        make_dioxolane(),
        make_test4(),
        make_Imide(),
        make_naphthalene(),
        make_ben1()
    }
};


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
    return result;
}

void recursive_modification(RWMOL_SPTR mol, std::string mol_smiles, FILTER_SPTR filter, Result& result_pool, ReasonFilterCount& reason_count, FunctionalGroupFilterCount &fg_count, int depth, int suppress_exceptions);
std::mutex print_mutex;
std::pair<int, std::string> replace_h_and_filt(RWMOL_SPTR mol, H_INDEX_LIST_SPTR h_index_list, FILTER_SPTR filter, Result& result_pool, ReasonFilterCount& reason_count, FunctionalGroupFilterCount& fg_count, ModificationType& mod, int depth, int suppress_exceptions)
{
    std::pair<int, std::string> r = std::make_pair(depth, std::get<0>(mod));
    std::vector<std::vector<uint8_t>>& h_list = h_index_list->operator[](std::get<1>(mod));
    unsigned int h_size = h_list.size();
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
                reason = (*filter)(result);
                if (reason.empty())
                {
                    result_pool.add_candidate(result, h_smiles);
                    if (depth >= 3 || std::get<0>(mod) == "dialkyl (num=0, h_num=0)" || std::get<0>(mod) == "dialkyl (num=1, h_num=1)" || std::get<0>(mod) == "dialkyl (num=3, h_num=1)")
                    {
                        // std::cout << h_smiles << std::endl;
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
    {
        return;
    }
    int max_mod = 0;
    for (auto& item : mod_list[depth -1])
        if (std::get<1>(item) > max_mod)
            max_mod = std::get<1>(item);
    
    H_INDEX_LIST_SPTR h_index_list = boost::make_shared<TripleIntVec>(1);
    h_index_list->reserve(max_mod + 1);
    for (int i = 1; i <= max_mod; ++i)
        h_index_list->emplace_back(enumerate_h(mol, i));

    unsigned int mod_size = mod_list[depth -1].size(), profile_count = 0, remaining;
    if (depth == MOD_DEPTH)
    {
        for (unsigned int i = mod_start_idx; i < mod_end_idx; ++i)
        {
            thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(reason_count), std::ref(fg_count), std::ref(mod_list[depth -1][i]), depth, suppress_exceptions));
            //thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(simple_mod_list[i]), depth, suppress_exceptions));
            std::this_thread::sleep_for(std::chrono::seconds(1));
            while ((remaining = thread_pool.get_task_num()) > 1.0e5 && !quit.load())
            {
                logging << "Stop submit task at mod " << i + 2 << "/" << mod_size << ", Task remain: " << remaining << ", Candidates: " << result_pool.get_candidate_size() << ", Filtered: " << result_pool.get_filtered_size() << ", Speed: " << result_pool.get_total_size() / difftime(time(NULL), start_time) << std::endl;
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < mod_list[depth -1].size(); ++i)
        	thread_pool.submit(depth, std::bind(replace_h_and_filt, mol, h_index_list, filter, std::ref(result_pool), std::ref(reason_count), std::ref(fg_count), std::ref(mod_list[depth -1][i]), depth, suppress_exceptions));
    }
}

void got_signal(int signal)
{
    logging << "Signal raised: " << signal << "\n";
    quit.store(true);
}

int main(int argc, char **argv)
{
    auto bg = std::chrono::high_resolution_clock::now(); //change
    std::signal(SIGINT, got_signal);
    start_time = time(NULL);
    if (argc > 3)
        PATH_PREFIX = argv[3];
    else {
        PATH_PREFIX += "output_cath/";
    }
    std::filesystem::create_directories(PATH_PREFIX);
    logging.open(PATH_PREFIX + argv[2] + ".log");

    //Initialise model configuration
    int suppress_exceptions = 1;

    std::vector<ModificationType> constraint_mod;
    for (auto& item : mod_list[MOD_DEPTH - 1])
        if (std::get<1>(item) <= 2)
            constraint_mod.emplace_back(item);
    mod_list[MOD_DEPTH - 1] = constraint_mod;
    logging << "# of Functional Group: " << mod_list[MOD_DEPTH - 1].size() << std::endl;

    int stride = std::ceil(mod_list[MOD_DEPTH - 1].size() / std::stoi(argv[1]));
    mod_start_idx = std::stoi(argv[2]) * stride;
    mod_end_idx = std::stoi(argv[2]) == std::stoi(argv[1]) - 1 ? mod_list[MOD_DEPTH - 1].size() : (std::stoi(argv[2]) + 1) * stride;
    logging << "Job: " << argv[2] << " (" << mod_start_idx << " - " << mod_end_idx << ")\n";

    const std::string cathinone = "O=C(c1ccccc1)[C@@H](N)C";
    RWMOL_SPTR mol(SmilesToMol(cathinone));
    MolOps::addHs(*mol);
    FILTER_SPTR filter = boost::make_shared<DruglikeFilter>();
    Result result(argv[2]);
    ReasonFilterCount reason_count(argv[2]);
    FunctionalGroupFilterCount fg_count(argv[2]);
    recursive_modification(mol, cathinone, filter, result, reason_count, fg_count, MOD_DEPTH, suppress_exceptions);
    std::this_thread::sleep_for(std::chrono::seconds(5));
    unsigned int remain = thread_pool.get_task_num(), profile_times = 0;
    logging << "Task add finished, size = " << remain << std::endl;
    while (remain != 0 && !quit.load())
    {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        remain = thread_pool.get_task_num();
        logging << "Task remain: " << remain << ", Candidates: " << result.get_candidate_size() << ", Filtered: " << result.get_filtered_size() << ", Speed: " << result.get_total_size() / difftime(time(NULL), start_time) << std::endl;
    }
    logging << "Wait threads to complete remaining tasks ...\n";
    thread_pool.terminate_all_thread();
    result.write_file();
    reason_count.writeFile();
    fg_count.writeFile();
    logging << "Write completed\n";
    logging << "\nTotal " << result.get_candidate_size() << " candidate\n";
    logging << "Total " << result.get_filtered_size() << " filtered\n";
    logging << "Average iteration per second: " << result.get_total_size() / difftime(time(NULL), start_time) << "\n";
    
    return 0;
}