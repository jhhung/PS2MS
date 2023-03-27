#include <iostream>
#include <set>
#include <string>
#include <map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <functional>
#include <atomic>
#include <csignal>
#include "DruglikeFilter.hpp"

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
using namespace RDKit; 

#include <boost/range/iterator_range.hpp>
#include <filesystem>
#include <gperftools/profiler.h>

std::vector<std::string> exceptions;
constexpr int MOD_DEPTH = 2;
constexpr int WRITE_FREQ = 100000;
std::atomic<bool> quit(false);

using ModificationType = 
    std::tuple<std::string, int, std::function<void(RWMOL_SPTR, std::vector<int>&)>>;
using FILTER_SPTR = boost::shared_ptr<Filter>;

struct Result
{
    Result() : size(0)
    {
        std::filesystem::remove("result_smiles.txt");
        std::filesystem::remove("filtered_smiles.txt");
    }

    ~Result()
    {
        output_result();
    }

    void append(std::string& smiles, std::string& reason)
    {
        if (reason != "")
            filtered_mol.emplace(smiles, reason);
        else
            candidate.emplace(smiles);
        ++size;
        if (size == WRITE_FREQ)
        {
            output_result();
            size = 0;
        }
    }

    void output_result()
    {
        std::ofstream outfile("result_smiles.txt");//, std::ios::app);
        for (auto& item : candidate)
            outfile << item << " " << item << "\n";
        outfile.close();
        // candidate.clear();

        outfile.open("filtered_smiles.txt");//, std::ios::app);
        for (auto& [key, value] : filtered_mol)
            outfile << std::setw(150) << std::left << key << "(" << value << ")\n";
        outfile.close();
        // filtered_mol.clear();
    }

    unsigned int size;
    std::set<std::string> candidate;
    std::map<std::string, std::string> filtered_mol;
};

std::vector<std::vector<int>> enumerate_h(RWMOL_SPTR mol, int h_num)
{
    std::vector<Atom*> all_h;
    all_h.reserve(mol->getNumAtoms());
    for (auto& atom : mol->atoms())
        if (atom->getSymbol() == "H")
            all_h.emplace_back(atom);
    
    // get permulation
    std::vector<std::vector<int>> index_list;
    std::vector<int> combination(h_num);
    std::vector<int> permutation(h_num);
    for (int j = 0; j < h_num; ++j)
        combination[j] = j + 1;
    do {
        index_list.push_back(combination);
    } while (std::next_permutation(combination.begin(), combination.end()));
    std::vector<int>::iterator first = combination.begin();
    std::vector<int>::iterator last = combination.end();
    unsigned int n = all_h.size(), r = h_num;
    while((*first) != (int)(n-r+1))
    {
        std::vector<int>::iterator mt = last;

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

void add_single_bond(RWMOL_SPTR mol, int src, int dst)
{
    Bond* target = mol->getBondBetweenAtoms(src, dst);
    if (target == nullptr)
        mol->addBond(src, dst, Bond::BondType::SINGLE);
    else
    {
        mol->removeBond(src, dst);
        switch ((int)target->getBondTypeAsDouble())
        {
            case 1:
                mol->addBond(src, dst, Bond::BondType::DOUBLE);
                break;
            case 2:
                mol->addBond(src, dst, Bond::BondType::TRIPLE);
                break;
            case 3:
                mol->addBond(src, dst, Bond::BondType::QUADRUPLE);
                break;
            case 4:
                mol->addBond(src, dst, Bond::BondType::QUINTUPLE);
                break;
            case 5:
                mol->addBond(src, dst, Bond::BondType::HEXTUPLE);
                break;
            default:
                std::cerr << "[ERROR] Unrecognized bond type: " << target->getBondTypeAsDouble() << "\n";
                exit(-1);
        }
    }
}

void remove_atoms(RWMOL_SPTR mol, std::vector<int>::iterator begin, std::vector<int>::iterator end)
{
    std::sort(begin, end, [] (int a, int b) {return a > b;});
    for (auto& it = begin; it != end; ++it)
        mol->removeAtom(*it);
}

inline int get_target(RWMOL_SPTR mol, int idx)
{
    return (*mol)[*(mol->getAtomNeighbors(mol->getAtomWithIdx(idx)).first)]->getIdx();
}

ModificationType make_acetal(int h_num = 0)
{
    char str[50];
    sprintf(str, "acetal (h_num=%d)", h_num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [h_num = h_num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("O"));
            int right_up_o, center_c;
            if (h_num < 3)
            {
                mol->replaceAtom(h_list[1], new Atom("O"));
                right_up_o = h_list[1];
                if (h_num < 2)
                {
                    mol->replaceAtom(h_list[2], new Atom("C"));
                    center_c = h_list[2];
                    if (h_num < 1)
                        add_single_bond(mol, center_c, get_target(mol, h_list[3]));
                }
                else
                    center_c = mol->addAtom(new Atom("C"));
            }
            else
            {
                center_c = mol->addAtom(new Atom("C"));
                right_up_o = mol->addAtom(new Atom("O"));
                int right_up_c = mol->addAtom(new Atom("C"));
                mol->addBond(right_up_o, right_up_c, Bond::BondType::SINGLE);
            }
            mol->addBond(center_c, h_list[0], Bond::BondType::SINGLE);
            mol->addBond(center_c, right_up_o, Bond::BondType::SINGLE);
            if (h_num == 0)
                mol->removeAtom(h_list[3]);
        };
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acetoxy()
{
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("O"));
            int center_c = mol->addAtom(new Atom("C"));
            mol->addBond(center_c, h_list[0], Bond::BondType::SINGLE);
            mol->addBond(center_c, mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(center_c, mol->addAtom(new Atom("C")), Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string("acetoxy ()"), 1, func);
}

ModificationType make_acetyl()
{
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], mol->addAtom(new Atom("C")), Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string("acetyl ()"), 1, func);
}

ModificationType make_acetylide(const char* matel)
{
    char str[50];
    sprintf(str, "acetylide (matel=%s)", matel);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [matel = matel] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int right_c = mol->addAtom(new Atom("C"));
            mol->addBond(h_list[0], right_c, Bond::BondType::TRIPLE);
            mol->addBond(right_c, mol->addAtom(new Atom(matel)), Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_acryloyl()
{
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int left_c = mol->addAtom(new Atom("C"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], left_c, Bond::BondType::SINGLE);
            mol->addBond(left_c, mol->addAtom(new Atom("C")), Bond::BondType::DOUBLE);
        };
    return std::make_tuple(std::string("acryloyl ()"), 1, func);
}

ModificationType make_acyl_azide()
{
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int right_n = mol->addAtom(new Atom("N"));
            Atom* atom = new Atom("N");
            atom->setFormalCharge(1);
            int right_n2 = mol->addAtom(atom);
            atom = new Atom("N");
            atom->setFormalCharge(-1);
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], right_n, Bond::BondType::SINGLE);
            mol->addBond(right_n, right_n2, Bond::BondType::DOUBLE);
            mol->addBond(right_n2, mol->addAtom(atom), Bond::BondType::DOUBLE);
        };
    return std::make_tuple(std::string("acyl azide ()"), 1, func);
}

ModificationType make_acyl(int h_num = 0)
{
    char str[50];
    sprintf(str, "acyl (h_num=%d)", h_num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [h_num = h_num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            if (h_num == 0)
                add_single_bond(mol, h_list[0], get_target(mol, h_list[1]));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            if (h_num == 0)
                mol->removeAtom(h_list[1]);
        };
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acyl_halide(const char* halide)
{
    char str[50];
    sprintf(str, "acyl halide (halide=%s)", halide);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [halide = halide] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], mol->addAtom(new Atom(halide)), Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_acylal(int h_num = 0)
{
    char str[50];
    sprintf(str, "acylal (h_num=%d)", h_num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [h_num = h_num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            // left
            mol->replaceAtom(h_list[0], new Atom("C"));
            int left_middle_o = mol->addAtom(new Atom("O"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], left_middle_o, Bond::BondType::SINGLE);
            // middle
            int middle_c;
            if (h_num < 2)
            {
                mol->replaceAtom(h_list[1], new Atom("C"));
                middle_c = h_list[1];
            }
            else
                middle_c = mol->addAtom(new Atom("C"));
            int right_middle_o = mol->addAtom(new Atom("O"));
            mol->addBond(middle_c, left_middle_o, Bond::BondType::SINGLE);
            mol->addBond(middle_c, right_middle_o, Bond::BondType::SINGLE);
            // right
            int right_c;
            if (h_num < 1)
            {
                mol->replaceAtom(h_list[2], new Atom("C"));
                right_c = h_list[2];
            }
            else
                right_c = mol->addAtom(new Atom("C"));
            mol->addBond(right_c, mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(right_c, right_middle_o, Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_acylhydrazine(int h_num = 0)
{
    char str[50];
    sprintf(str, "acylhydrazine (h_num=%d)", h_num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [h_num = h_num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int middle_n, right_n;
            if (h_num < 3)
            {
                mol->replaceAtom(h_list[1], new Atom("N"));
                middle_n = h_list[1];
                if (h_num < 2)
                {
                    mol->replaceAtom(h_list[2], new Atom("N"));
                    right_n = h_list[2];
                    if (h_num < 1)
                        add_single_bond(mol, right_n, get_target(mol, h_list[3]));
                }
                else
                    right_n = mol->addAtom(new Atom("N"));
            }
            else
            {
                middle_n = mol->addAtom(new Atom("N"));
                right_n = mol->addAtom(new Atom("N"));
            }
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], middle_n, Bond::BondType::SINGLE);
            mol->addBond(middle_n, right_n, Bond::BondType::SINGLE);
            if (h_num == 0)
                mol->removeAtom(h_list[3]);
        };
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acyloin(int h_num = 0)
{
    char str[50];
    sprintf(str, "acyloin (h_num=%d)", h_num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [h_num = h_num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int right_c;
            if (h_num == 0)
            {
                mol->replaceAtom(h_list[1], new Atom("C"));
                right_c = h_list[1];
            }
            else
                right_c = mol->addAtom(new Atom("C"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], right_c, Bond::BondType::SINGLE);
            mol->addBond(right_c, mol->addAtom(new Atom("O")), Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_acylsilane(int h_num = 0)
{
    char str[50];
    sprintf(str, "acylsilane (h_num=%d)", h_num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [h_num = h_num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int si;
            if (h_num < 3)
            {
                mol->replaceAtom(h_list[1], new Atom("Si"));
                si = h_list[1];
                if (h_num < 2)
                {
                    add_single_bond(mol, si, get_target(mol, h_list[2]));
                    if (h_num < 1)
                        add_single_bond(mol, si, get_target(mol, h_list[3]));
                }
                else
                    si = mol->addAtom(new Atom("Si"));
            }
            else
                si = mol->addAtom(new Atom("Si"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], si, Bond::BondType::SINGLE);
            if (h_num < 1)
                remove_atoms(mol, h_list.begin() + 2, h_list.end());
            else if (h_num < 2)
                mol->removeAtom(h_list[2]);
        };
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acylurea()
{
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int middle_n = mol->addAtom(new Atom("N"));
            int right_c = mol->addAtom(new Atom("C"));
            mol->addBond(h_list[0], mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(h_list[0], middle_n, Bond::BondType::SINGLE);
            mol->addBond(middle_n, right_c, Bond::BondType::SINGLE);
            mol->addBond(right_c, mol->addAtom(new Atom("O")), Bond::BondType::DOUBLE);
            mol->addBond(right_c, mol->addAtom(new Atom("N")), Bond::BondType::SINGLE);
        };
    return std::make_tuple(std::string("acylurea ()"), 1, func);
}

ModificationType make_alkyl(int num)
{
    char str[50];
    sprintf(str, "alkyl (num=%d)", num);
    
    std::function<void(RWMOL_SPTR, std::vector<int>&)> func = 
        [num = num] (RWMOL_SPTR mol, std::vector<int>& h_list)
        {
            mol->replaceAtom(h_list[0], new Atom("C"));
            int new_c;
            int last_c = h_list[0];
            for (int i = 0; i < num - 1; ++i)
            {
                new_c = mol->addAtom(new Atom("C"));
                mol->addBond(last_c, new_c, Bond::BondType::SINGLE);
            }
        };
    return std::make_tuple(std::string(str), 1, func);
}

std::vector<ModificationType> get_modification_list()
{
    return {
        make_alkyl(1),
        make_alkyl(2),
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
        make_acryloyl(),
        make_acyl_azide(),
        make_acyl(),
        make_acyl(1),
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
        make_acylurea()
    };
}

void add_h(RWMOL_SPTR mol)
{
    if (mol->needsUpdatePropertyCache())
        mol->updatePropertyCache();
    MolOps::addHs(*mol);
}

RWMOL_SPTR do_modification(
    std::vector<int>& h_list, RWMOL_SPTR parent, const std::function<void(RWMOL_SPTR, std::vector<int>&)>& mod)
{
    RWMOL_SPTR result(new RWMol(*parent, true));
    mod(result, h_list);
    // if (result->needsUpdatePropertyCache())
    //     result->updatePropertyCache();
    // MolOps::removeHs(*result);
    add_h(result);
    return result;
}

void recursive_modification(RWMOL_SPTR mol, std::string mol_smiles, FILTER_SPTR filter, Result& result_pool, int depth)
{
    if (depth == 0)
        return;
    
    std::vector<ModificationType> mod_list = get_modification_list();
    int max_mod = 0;
    for (auto& item : mod_list)
        if (std::get<1>(item) > max_mod)
            max_mod = std::get<1>(item);
    std::vector<std::vector<std::vector<int>>> h_index_list(1);
    h_index_list.reserve(max_mod + 1);
    for (int i = 1; i <= max_mod; ++i)
        h_index_list.emplace_back(enumerate_h(mol, i));

    unsigned int mod_size = mod_list.size();
    unsigned int h_size;
    unsigned int mod_print_interval = (unsigned int)(mod_size / 10);
    unsigned int h_print_interval;
    std::vector<std::vector<int>> h_list;
    std::string filt_result;
    for (unsigned int i = 0; i < mod_size; ++i)
    {
        if (i % mod_print_interval == 0)
        {
            std::string indent(2 * 2 * (MOD_DEPTH - depth), ' ');
            std::cerr << indent << "Mol = " << mol_smiles << ", Mod = " << std::get<0>(mod_list[i])
                      << " (" << i << " / " << mod_size << " -> " << i * 100 / mod_size << "%)\n";
        }
        h_list = h_index_list[std::get<1>(mod_list[i])];
        h_size = h_list.size();
        h_print_interval = h_size / 10 == 0 ? 1 : (unsigned int)(h_size / 10);
        for (unsigned int j = 0; j < h_size; ++j)
        {
            RWMOL_SPTR result = do_modification(h_list[j], mol, std::get<2>(mod_list[i]));
            std::string smiles = MolToSmiles(*result);
            if (result_pool.filtered_mol.find(smiles) == result_pool.filtered_mol.end()
                && result_pool.candidate.find(smiles) == result_pool.candidate.end())
            {
                // do filter
                result->setProp("filtered", false);
                result->setProp("failedfilter", "");
                filt_result = (*filter)(result);
                if (j % h_print_interval == 0)
                {
                    std::string indent(2 * 2 * (MOD_DEPTH - depth) + 1, ' ');
                    std::cerr << indent << smiles << ": " << (filt_result.empty() ? "Pass!" : filt_result)
                            << " (" << j << " / " << h_size << " -> " << j * 100 / h_size << "%)\n";
                }
                // std::cout << smiles << ": " << (*filter)(result) << std::endl;
                std::string reason;
                result->getProp("failedfilter", reason);
                result_pool.append(smiles, reason);
                add_h(result);
            }
            // add_h(result); // filter will suppress Hs
            recursive_modification(result, smiles, filter, result_pool, depth - 1);
            if (quit.load())
                return;
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

    const std::string amphetamine = "NC(CC1=CC=CC=C1)C";
    RWMOL_SPTR mol(SmilesToMol(amphetamine));
    MolOps::addHs(*mol);
    print_molecule(mol);
    FILTER_SPTR filter(new DruglikeFilter());
    Result result;
    ProfilerStart("filter.prof");
    recursive_modification(mol, amphetamine, filter, result, MOD_DEPTH);
    ProfilerStop();
    std::cerr << "Total " << result.candidate.size() << " candidate\n";
    std::cerr << "Total " << result.filtered_mol.size() << " filtered\n";
    return 0;
}