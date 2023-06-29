#include <vector>
#include <string>
#include <tuple>
#include <functional>
#include <algorithm>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
using namespace RDKit; 

#include <boost/range/iterator_range.hpp>
#include <boost/shared_ptr.hpp>

using ModificationType = 
    std::tuple<std::string, int, std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>>>;

void add_bond(RWMOL_SPTR mol, int src, int dst, int type = 1, Bond::BondDir direction = Bond::BondDir::NONE)
{
    Bond* target = mol->getBondBetweenAtoms(src, dst);
    if (target == nullptr)
        mol->getBondWithIdx(mol->addBond(src, dst, static_cast<Bond::BondType>(type)) - 1)->setBondDir(direction);
    else
    {
        double origin_type = target->getBondTypeAsDouble();
        mol->removeBond(src, dst);
        if (origin_type + type <= 6)
            mol->getBondWithIdx(mol->addBond(src, dst, static_cast<Bond::BondType>((int)origin_type + type)) - 1)->setBondDir(direction);
        else
        {
            std::cerr << "[ERROR] Unrecognized bond type: " << origin_type << "\n";
            exit(-1);
        }
    }
}

void remove_atoms(RWMOL_SPTR mol, std::vector<uint8_t>::iterator begin, std::vector<uint8_t>::iterator end)
{
    std::vector<int> copy(begin, end);
    std::sort(copy.begin(), copy.end(), [] (int a, int b) {return a > b;});
    for (auto it = copy.begin(); it != copy.end(); ++it)
        mol->removeAtom(*it);
}

inline int get_target(RWMOL_SPTR mol, int idx)
{
    return (*mol)[*(mol->getAtomNeighbors(mol->getAtomWithIdx(idx)).first)]->getIdx();
}

std::vector<uint8_t> calculate_link_point(RWMOL_SPTR mol, std::vector<int>& pattern, std::vector<uint8_t>& h_list)
{
    Atom atom_h("H");
    // don't directly modify pattern, or it'll be changed when next time you call this lambda
    std::vector<uint8_t> link_point(pattern.size());
    auto h_list_it = h_list.begin();
    for (int i = 0; i < pattern.size(); ++i)
    {
        if (pattern[i] == -1)
        {
            link_point[i] = mol->addAtom(&atom_h);
            mol->addAtom(&atom_h);
            add_bond(mol, link_point[i], link_point[i] + 1);
        }
        else
            link_point[i] = *(h_list_it++);
    }
    return link_point;
}

ModificationType make_acetal(int h_num = 0)
{
    // TODO: 目前每個h接的點位的位置固定，要改成接的點位也要改變
    // 先製造h_num個H-H，將這些h與要替換的h做排列，lambda再一視同仁的接上每個R，就可以得到接點shffule的目的
    char str[50];
    sprintf(str, "acetal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                mol->replaceAtom(link_point[1], &atom_o);
                // right_up_o = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // center_c = link_point[2];
                add_bond(mol, link_point[2], link_point[0]);
                add_bond(mol, link_point[2], link_point[1]);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acetoxy()
{
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            int center_c = mol->addAtom(&atom_c);
            add_bond(mol, center_c, h_list[0]);
            add_bond(mol, center_c, mol->addAtom(&atom_o), 2);
            add_bond(mol, center_c, mol->addAtom(&atom_c));
        }
    };
    return std::make_tuple(std::string("acetoxy ()"), 1, func);
}

ModificationType make_acetyl()
{
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_c));
        }
    };
    return std::make_tuple(std::string("acetyl ()"), 1, func);
}

ModificationType make_acetylide(const char* matel)
{
    char str[50];
    sprintf(str, "acetylide (matel=%s)", matel);
    
    Atom atom_c("C");
    Atom atom_matel(matel);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_matel = atom_matel] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int right_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], right_c, 3);
            add_bond(mol, right_c, mol->addAtom(&atom_matel));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_acid_anhydride(int h_num = 0)
{
    char str[50];
    sprintf(str, "acid anhydride (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o), 2);
                int middle_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], middle_o);
                add_bond(mol, link_point[1], middle_o);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_acryloyl()
{
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int left_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], left_c);
            add_bond(mol, left_c, mol->addAtom(&atom_c), 2);
        }
    };
    return std::make_tuple(std::string("acryloyl ()"), 1, func);
}

ModificationType make_acyl_azide()
{
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_n2("N");
    atom_n2.setFormalCharge(1);
    Atom atom_n3("N");
    atom_n3.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c, atom_n = atom_n, atom_n2 = atom_n2, atom_n3 = atom_n3] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int right_n = mol->addAtom(&atom_n);
            int right_n2 = mol->addAtom(&atom_n2);
            int right_n3 = mol->addAtom(&atom_n3);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], right_n);
            add_bond(mol, right_n, right_n2, 2);
            add_bond(mol, right_n2, right_n3, 2);
        }
    };
    return std::make_tuple(std::string("acyl azide ()"), 1, func);
}

ModificationType make_acyl_chloride()
{
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_cl));
        }
    };
    return std::make_tuple(std::string("acyl chloride ()"), 1, func);
}

ModificationType make_acyl_halide(const char* halide)
{
    char str[50];
    sprintf(str, "acyl halide (halide=%s)", halide);
    
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_halide(halide);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c, atom_halide = atom_halide] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_halide));
            // std::cerr <<  mol->getAtomWithIdx(h_list[0])->getSymbol() << std::endl;
            // std::cerr <<  mol->getAtomWithIdx(23)->getSymbol() << std::endl;
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_acylal(int h_num = 0)
{
    char str[50];
    sprintf(str, "acylal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                // left
                mol->replaceAtom(link_point[0], &atom_c);
                int left_middle_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], left_middle_o);
                // middle
                mol->replaceAtom(link_point[1], &atom_c);
                // middle_c = link_point[1];
                int right_middle_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[1], left_middle_o);
                add_bond(mol, link_point[1], right_middle_o);
                // right
                mol->replaceAtom(link_point[2], &atom_c);
                // right_c = link_point[2];
                add_bond(mol, link_point[2], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[2], right_middle_o);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_acylhydrazine(int h_num = 0)
{
    char str[50];
    sprintf(str, "acylhydrazine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_n);
                // middle_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // right_n = link_point[2];
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], link_point[2]);
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acyloin(int h_num = 0)
{
    char str[50];
    sprintf(str, "acyloin (h_num=%d)", h_num);
    
    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);

    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_acylsilane(int h_num = 0)
{
    char str[50];
    sprintf(str, "acylsilane (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_si("Si");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c, atom_si = atom_si] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_si);
                // si = link_point[1];
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_acylurea()
{
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int middle_n = mol->addAtom(&atom_n);
            int right_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], middle_n);
            add_bond(mol, middle_n, right_c);
            add_bond(mol, right_c, mol->addAtom(&atom_o), 2);
            add_bond(mol, right_c, mol->addAtom(&atom_n));
        }
    };
    return std::make_tuple(std::string("acylurea ()"), 1, func);
}

///////
ModificationType make_alcohol()
{
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            add_bond(mol, h_list[0], mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("alcohol ()"), 1, func);
}

ModificationType make_aldehyde()
{
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("aldehyde ()"), 1, func);
}

ModificationType make_aldimine(int h_num = 0)
{
    char str[50];
    sprintf(str, "aldimine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_n);
                add_bond(mol, link_point[0], link_point[1], 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_alkene(int h_num = 0)
{
    char str[50];
    sprintf(str, "alkene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_c);
                add_bond(mol, link_point[0], link_point[1], 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

// ModificationType make_alkene()
// {
//     Atom atom_c("C");
//     std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
//         [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
//         {
//             int target_idx = get_target(mol, h_list[1]);
//             if (get_target(mol, h_list[0]) == target_idx)
//             {
//                 mol->replaceAtom(h_list[0], &atom_c);
//                 add_bond(mol, h_list[0], target_idx);
//                 mol->removeAtom(h_list[1]);
//             }
//         }
//     };
//     return std::make_tuple(std::string("alkene ()"), 2, func);
// }

ModificationType make_alkoxide(int h_num = 0)
{
    char str[50];
    sprintf(str, "alkoxide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

// ModificationType make_alkoxide()
// {
//     Atom atom_o("O");
//     atom_o.setFormalCharge(-1);
//     std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
//         [atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
//         {
//             mol->replaceAtom(h_list[0], &atom_o);
//         }
//     };
//     return std::make_tuple(std::string("alkoxide ()"), 1, func);
// }

ModificationType make_alkyl(int num, int h_num = 0)
{
    char str[50];
    sprintf(str, "alkyl (num=%d, h_num=%d)", num, h_num);

    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [num = num, h_num = h_num, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int new_c;
            int last_c = h_list[0];
            for (int i = 0; i < num - 2; ++i)
            {
                new_c = mol->addAtom(&atom_c);
                add_bond(mol, last_c, new_c);
                last_c = new_c;
            }
            if (h_num == 1 && num!=0)
            {
                new_c = mol->addAtom(&atom_c);
                add_bond(mol, last_c, new_c);
            }
            else if (h_num==0 && num!=0)
            {
                mol->replaceAtom(h_list[1], &atom_c);
                add_bond(mol, last_c, h_list[1]);
            }
        }
    };
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_alkyl_nitrites()
{
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            int middle_n = mol->addAtom(&atom_n);
            add_bond(mol, h_list[0], middle_n);
            add_bond(mol, middle_n, mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("alkyl nitrites ()"), 1, func);
}

ModificationType make_alkyne(int h_num = 0)
{
    char str[50];
    sprintf(str, "alkyne (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_c);
                add_bond(mol, link_point[0], link_point[1], 3);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

// ModificationType make_alkyne()
// {
//     Atom atom_c("C");
//     std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
//         [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
//         {
//             int target_idx = get_target(mol, h_list[1]);
//             if (get_target(mol, h_list[0]) == target_idx && get_target(mol, h_list[2]) == target_idx)
//             {
//                 mol->replaceAtom(h_list[0], &atom_c);
//                 add_bond(mol, h_list[0], target_idx, 2);
//                 remove_atoms(mol, h_list.begin() + 1, h_list.end());
//             }
//         }
//     };
//     return std::make_tuple(std::string("alkyne ()"), 3, func);
// }

ModificationType make_allyl()
{
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int middle_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], middle_c);
            add_bond(mol, middle_c, mol->addAtom(&atom_c), 2);
        }
    };
    return std::make_tuple(std::string("allyl ()"), 1, func);
}

ModificationType make_amide(int h_num = 0)
{
    char str[50];
    sprintf(str, "amide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                // left
                mol->replaceAtom(link_point[0], &atom_c);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                // right
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1]);
                mol->removeAtom(link_point[2]);
                
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_amidine()
{
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_n), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_n));
        }
    };
    return std::make_tuple(std::string("amidine ()"), 1, func);
}

ModificationType make_amidrazone(bool amide_aydrazone = false)
{
    char str[50];
    sprintf(str, "amidine (amide_aydrazone=%s)", amide_aydrazone ? "true" : "false");
    
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [amide_aydrazone = amide_aydrazone, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int middle_n = mol->addAtom(&atom_n);
            if (!amide_aydrazone)
            {
                add_bond(mol, h_list[0], mol->addAtom(&atom_n), 2);
                add_bond(mol, h_list[0], middle_n);
            }
            else
            {
                add_bond(mol, h_list[0], mol->addAtom(&atom_n));
                add_bond(mol, h_list[0], middle_n, 2);
            }
            add_bond(mol, middle_n, mol->addAtom(&atom_n));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_aminal(int h_num = 0)
{
    char str[50];
    sprintf(str, "aminal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                mol->replaceAtom(link_point[1], &atom_n);
                // right_up_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // center_c = link_point[2];
                add_bond(mol, link_point[2], link_point[0]);
                add_bond(mol, link_point[2], link_point[1]);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_amine(int h_num = 0)
{
    char str[50];
    sprintf(str, "amine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_n("N");
    // atom_n.setNumRadicalElectrons(2);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINDASH);
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_amine_oxide(int h_num = 0)
{
    char str[50];
    sprintf(str, "amine oxide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_o("O");
    atom_n.setFormalCharge(1);
    atom_o.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o));
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_aminophosphine_r3(int p_h_num = 0)
{
    char str[50];
    sprintf(str, "aminophosphine (p_r_num=3, p_h_num=%d, n_h_num=0)", p_h_num);

    std::vector<int> permutation(p_h_num, -1);
    permutation.insert(permutation.end(), 3 - p_h_num, 1);
    
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINDASH);
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - p_h_num, func);
}

ModificationType make_aminophosphine_r2(int p_h_num = 0, int n_h_num = 0)
{
    char str[50];
    sprintf(str, "aminophosphine (p_r_num=2, p_h_num=%d, n_h_num=%d)", p_h_num, n_h_num);

    std::vector<int> permutation(p_h_num + n_h_num, -1);
    permutation.insert(permutation.end(), (2 - p_h_num) + (2 - n_h_num), 1);
    
    Atom atom_n("N");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    for (int bond_dir = 0; bond_dir < 3; ++bond_dir) // 0: NONE, 1: BEGINWEDGE, 2: BEGINDASH
    {
        do
        {
            func.emplace_back(
                [permutation = permutation, bond_dir = bond_dir, atom_n = atom_n, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
                {
                    std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                    int p = mol->addAtom(&atom_p);
                    // has N
                    mol->replaceAtom(link_point[0], &atom_n);
                    // n = link_point[0];
                    add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                    add_bond(mol, p, link_point[0], 1, (Bond::BondDir)(bond_dir % 3));
                    // no N
                    add_bond(mol, p, get_target(mol, link_point[2]), 1, (Bond::BondDir)((bond_dir + 1) % 3));
                    add_bond(mol, p, get_target(mol, link_point[3]), 1, (Bond::BondDir)((bond_dir + 2) % 3));
                    remove_atoms(mol, link_point.begin() + 1, link_point.end());
                }
            );
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    return std::make_tuple(std::string(str), (2 - p_h_num) + (2 - n_h_num), func);
}

ModificationType make_aminophosphine_r1(int p_h_num = 0, int n_h_num = 0)
{
    char str[50];
    sprintf(str, "aminophosphine (p_r_num=1, p_h_num=%d, n_h_num=%d)", p_h_num, n_h_num);

    std::vector<int> permutation(p_h_num + n_h_num, -1);
    permutation.insert(permutation.end(), (1 - p_h_num) + (4 - n_h_num), 1);
    
    Atom atom_n("N");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    for (int bond_dir = 0; bond_dir < 3; ++bond_dir) // 0: NONE, 1: BEGINWEDGE, 2: BEGINDASH
    {
        do
        {
            func.emplace_back(
                [permutation = permutation, bond_dir = bond_dir, atom_n = atom_n, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
                {
                    std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                    int p = mol->addAtom(&atom_p);
                    // has R
                    add_bond(mol, p, get_target(mol, link_point[2]), 1, (Bond::BondDir)(bond_dir % 3));
                    // no R
                    mol->replaceAtom(link_point[0], &atom_n);
                    // n = link_point[0];
                    add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                    add_bond(mol, p, link_point[0], 1, (Bond::BondDir)((bond_dir + 1) % 3));
                    mol->replaceAtom(link_point[1], &atom_n);
                    // n1 = link_point[1];
                    add_bond(mol, link_point[1], get_target(mol, link_point[4]));
                    add_bond(mol, p, link_point[1], 1, (Bond::BondDir)((bond_dir + 2) % 3));
                    remove_atoms(mol, link_point.begin() + 2, link_point.end());
                }
            );
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    return std::make_tuple(std::string(str), (1 - p_h_num) + (4 - n_h_num), func);
}

ModificationType make_aminophosphine_r0(int n_h_num = 0)
{
    char str[50];
    sprintf(str, "aminophosphine (p_r_num=0, p_h_num=0, n_h_num=%d)", n_h_num);

    std::vector<int> permutation(n_h_num, -1);
    permutation.insert(permutation.end(), 6 - n_h_num, 1);
    
    Atom atom_n("N");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int p = mol->addAtom(&atom_p);
                for (int bond_dir = 0; bond_dir < 3; ++bond_dir) // 0: NONE, 1: BEGINWEDGE, 2: BEGINDASH
                {
                    mol->replaceAtom(link_point[bond_dir], &atom_n);
                    // n = link_point[bond_dir];
                    add_bond(mol, link_point[bond_dir], get_target(mol, link_point[3 + bond_dir]));
                    add_bond(mol, p, link_point[bond_dir], 1, (Bond::BondDir)(bond_dir));
                }
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 6 - n_h_num, func);
}

ModificationType make_aminoxyl(int h_num = 0)
{
    char str[50];
    sprintf(str, "aminoxyl (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_o("O");
    atom_o.setNumRadicalElectrons(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_azide(bool triple_bond = true)
{
    char str[50];
    sprintf(str, "azide (triple_bond=%s)", triple_bond ? "true" : "false");

    Atom atom_n("N");
    if (triple_bond)
        atom_n.setFormalCharge(-1);
    Atom atom_n2("N");
    atom_n2.setFormalCharge(1);
    Atom atom_n3("N");
    if (!triple_bond)
        atom_n3.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [triple_bond = triple_bond, atom_n = atom_n, atom_n2 = atom_n2, atom_n3 = atom_n3] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            // left_n = h_list[0]
            int middle_n = mol->addAtom(&atom_n2);
            int right_n = mol->addAtom(&atom_n3);
            add_bond(mol, h_list[0], middle_n, triple_bond ? 1 : 2);
            add_bond(mol, middle_n, right_n, triple_bond ? 3 : 2);
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_azine(int h_num = 0)
{
    char str[50];
    sprintf(str, "azine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                int left_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], left_n, 2);
                int right_n = mol->addAtom(&atom_n);
                add_bond(mol, left_n, right_n);
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                add_bond(mol, link_point[1], right_n, 2);
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_aziridine(int h_num = 0)
{
    char str[50];
    sprintf(str, "aziridine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 5 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // up_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // left_c = link_point[1];
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                mol->replaceAtom(link_point[2], &atom_c);
                // right_c = link_point[2];
                add_bond(mol, link_point[2], get_target(mol, link_point[4]));
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], link_point[0]);
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 5 - h_num, func);
}

ModificationType make_azo(int h_num = 0)
{
    char str[50];
    sprintf(str, "azo (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                mol->replaceAtom(link_point[1], &atom_n);
                add_bond(mol, link_point[0], link_point[1], 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_azole_3_hands_1_atom(const char* name, int h_num = 0)
{
    char str[50];
    sprintf(str, "azole with 3 hands atom (atom=%s, num=1, h_num=%d)", name, h_num);

    std::vector<std::vector<bool>> compound_template = {
        {false, true, false, false}, // Imidazole
        {true, false, false, false}, // Pyrazole
    };

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    for (auto& compound : compound_template)
    {
        do
        {
            func.emplace_back(
                [permutation = permutation, bound_n = compound, atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
                {
                    std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                    mol->replaceAtom(link_point[0], &atom);
                    // middle_atom = link_point[0];
                    int added_n = mol->addAtom(&atom_n);
                    for (int i = 1; i < 4; ++i)
                        mol->replaceAtom(link_point[i], &atom_c);
                    int current_idx = link_point[0];
                    int c_idx = 1;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (bound_n[i])
                        {
                            add_bond(mol, current_idx, added_n, i % 2 == 1 ? 2 : 1);
                            current_idx = added_n;
                        }
                        else
                        {
                            add_bond(mol, current_idx, link_point[c_idx], i % 2 == 1 ? 2 : 1);
                            current_idx = link_point[c_idx++];
                        }
                    }
                    add_bond(mol, current_idx, link_point[0]);
                }
            );
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_azole_3_hands_2_atoms(const char* name, int h_num = 0)
{
    char str[50];
    sprintf(str, "azole with 3 hands atom (atom=%s, num=2, h_num=%d)", name, h_num);

    std::vector<std::vector<bool>> compound_template = {
        {true, true, false, false}, // 1,2,3-Triazole
        {true, false, true, false}, // 1,2,4-Triazole
        {true, false, false, true}, // 1,2,5-Triazole
        {false, true, true, false}, // 1,3,4-Triazole
    };

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    for (auto& compound : compound_template)
    {
        do
        {
            func.emplace_back(
                [permutation = permutation, bound_n = compound, atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
                {
                    std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                    mol->replaceAtom(link_point[0], &atom);
                    // middle_atom = link_point[0];
                    int added_n[2] = {mol->addAtom(&atom_n), mol->addAtom(&atom_n)};
                    for (int i = 1; i < 3; ++i)
                        mol->replaceAtom(link_point[i], &atom_c);
                    int current_idx = link_point[0];
                    int c_idx = 1, n_idx = 0;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (bound_n[i])
                        {
                            add_bond(mol, current_idx, added_n[n_idx], i % 2 == 1 ? 2 : 1);
                            current_idx = added_n[n_idx++];
                        }
                        else
                        {
                            add_bond(mol, current_idx, link_point[c_idx], i % 2 == 1 ? 2 : 1);
                            current_idx = link_point[c_idx++];
                        }
                    }
                    add_bond(mol, current_idx, link_point[0]);
                }
            );
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_azole_3_hands_3_atoms(const char* name, int h_num = 0)
{
    char str[50];
    sprintf(str, "azole with 3 hands atom (atom=%s, num=3, h_num=%d)", name, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom);
                 // middle_atom = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                int added_n[3] = {mol->addAtom(&atom_n), mol->addAtom(&atom_n), mol->addAtom(&atom_n)};
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], added_n[0], 2);
                add_bond(mol, added_n[0], added_n[1]);
                add_bond(mol, added_n[1], added_n[2], 2);
                add_bond(mol, added_n[2], link_point[0]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_azole_3_hands_4_atoms(const char* name)
{
    char str[50];
    sprintf(str, "azole with 3 hands atom (atom=%s, num=4, h_num=0)", name);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = 
    {
        [atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom);
            // middle_atom = link_point[0];
            int added_n[4] = {mol->addAtom(&atom_n), mol->addAtom(&atom_n), mol->addAtom(&atom_n), mol->addAtom(&atom_n)};
            add_bond(mol, h_list[0], added_n[0]);
            add_bond(mol, added_n[0], added_n[1], 2);
            add_bond(mol, added_n[1], added_n[2]);
            add_bond(mol, added_n[2], added_n[3], 2);
            add_bond(mol, added_n[3], h_list[0]);
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_azole_2_hands_1_atom(const char* name, int h_num = 0)
{
    char str[50];
    sprintf(str, "azole with 2 hands atom (atom=%s, num=1, h_num=%d)", name, h_num);

    std::vector<std::vector<bool>> compound_template = {
        {false, true, false, false}, // Oxazole
        {true, false, false, false}, // Isoxazole
    };

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    for (auto& compound : compound_template)
    {
        do
        {
            func.emplace_back(
                [permutation = permutation, bound_n = compound, atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
                {
                    std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                    int middle_atom = mol->addAtom(&atom);
                    int added_n = mol->addAtom(&atom_n);
                    for (int i = 0; i < 3; ++i)
                        mol->replaceAtom(link_point[i], &atom_c);
                    int current_idx = middle_atom;
                    int c_idx = 0;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (bound_n[i])
                        {
                            add_bond(mol, current_idx, added_n, i % 2 == 1 ? 2 : 1);
                            current_idx = added_n;
                        }
                        else
                        {
                            add_bond(mol, current_idx, link_point[c_idx], i % 2 == 1 ? 2 : 1);
                            current_idx = link_point[c_idx++];
                        }
                    }
                    add_bond(mol, current_idx, middle_atom);
                }
            );
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_azole_2_hands_2_atoms(const char* name, int h_num = 0)
{
    char str[50];
    sprintf(str, "azole with 2 hands atom (atom=%s, num=2, h_num=%d)", name, h_num);

    std::vector<std::vector<bool>> compound_template = {
        {true, true, false, false}, // 1,2,3-Oxadiazole
        {true, false, true, false}, // 1,2,4-Oxadiazole
        {true, false, false, true}, // 1,2,5-Oxadiazole
        {false, true, true, false}, // 1,3,4-Oxadiazole
    };

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    for (auto& compound : compound_template)
    {
        do
        {
            func.emplace_back(
                [permutation = permutation, bound_n = compound, atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
                {
                    std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                    int middle_atom = mol->addAtom(&atom);
                    int added_n[2] = {mol->addAtom(&atom_n), mol->addAtom(&atom_n)};
                    for (int i = 0; i < 2; ++i)
                        mol->replaceAtom(link_point[i], &atom_c);
                    int current_idx = middle_atom;
                    int c_idx = 0, n_idx = 0;
                    for (int i = 0; i < 4; ++i)
                    {
                        if (bound_n[i])
                        {
                            add_bond(mol, current_idx, added_n[n_idx], i % 2 == 1 ? 2 : 1);
                            current_idx = added_n[n_idx++];
                        }
                        else
                        {
                            add_bond(mol, current_idx, link_point[c_idx], i % 2 == 1 ? 2 : 1);
                            current_idx = link_point[c_idx++];
                        }
                    }
                    add_bond(mol, current_idx, middle_atom);
                }
            );
        } while (std::next_permutation(permutation.begin(), permutation.end()));
    }
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_azole_2_hands_3_atoms(const char* name)
{
    char str[50];
    sprintf(str, "azole with 2 hands atom (atom=%s, num=3, h_num=0)", name);

    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom(name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_n = atom_n, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            int middle_atom = mol->addAtom(&atom);
            mol->replaceAtom(h_list[0], &atom_c);
            int added_n[3] = {mol->addAtom(&atom_n), mol->addAtom(&atom_n), mol->addAtom(&atom_n)};
            add_bond(mol, middle_atom, h_list[0]);
            add_bond(mol, h_list[0], added_n[0], 2);
            add_bond(mol, added_n[0], added_n[1]);
            add_bond(mol, added_n[1], added_n[2], 2);
            add_bond(mol, added_n[2], middle_atom);
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_azoxy(int h_num = 0)
{
    char str[50];
    sprintf(str, "azoxy (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_n2("N");
    atom_n2.setFormalCharge(1);
    Atom atom_o("O");
    atom_o.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_n2 = atom_n2, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                mol->replaceAtom(link_point[1], &atom_n2);
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_basic_aluminium_r2(int h_num = 0)
{
    char str[50];
    sprintf(str, "basic aluminium (r_num=2, h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_al("Al");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_al = atom_al, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_al);
                // al = link_point[0];
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], o);
                add_bond(mol, o, mol->addAtom(&atom_h));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_basic_aluminium_r1()
{
    Atom atom_al("Al");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_al = atom_al, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_al);
            // al = h_list[0];
            int o = mol->addAtom(&atom_o);
            int o2 = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], o);
            add_bond(mol, h_list[0], o2);
            add_bond(mol, o, mol->addAtom(&atom_h));
            add_bond(mol, o2, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("basic aluminium (r_num=1, h_num=0)"), 1, func);
}

ModificationType make_benzylidene_acetal(int h_num = 0)
{
    char str[50];
    sprintf(str, "benzylidene acetal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                int right_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], right_c);
                add_bond(mol, link_point[1], right_c);
                int first_c = mol->addAtom(&atom_c);
                int next_c, last_c = first_c;
                add_bond(mol, right_c, first_c);
                // C6H5
                for (int i = 0; i < 5; ++i)
                {
                    next_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, next_c, Bond::BondType::AROMATIC);
                    last_c = next_c;
                }
                mol->addBond(last_c, first_c, Bond::BondType::AROMATIC);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_bisthiosemicarbazone(int h_num = 0)
{
    char str[50];
    sprintf(str, "bisthiosemicarbazone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c1 = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                int left_n1 = mol->addAtom(&atom_n), right_n1 = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], left_n1, 2);
                add_bond(mol, link_point[1], right_n1, 2);
                int left_n2 = mol->addAtom(&atom_n), right_n2 = mol->addAtom(&atom_n);
                add_bond(mol, left_n1, left_n2);
                add_bond(mol, right_n1, right_n2);
                int left_c2 = mol->addAtom(&atom_c), right_c2 = mol->addAtom(&atom_c);
                add_bond(mol, left_n2, left_c2);
                add_bond(mol, right_n2, right_c2);
                add_bond(mol, left_c2, mol->addAtom(&atom_s), 2);
                add_bond(mol, right_c2, mol->addAtom(&atom_s), 2);
                add_bond(mol, left_c2, mol->addAtom(&atom_n));
                add_bond(mol, right_c2, mol->addAtom(&atom_n));
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_biuret(int h_num = 0)
{
    char str[50];
    sprintf(str, "bisthiosemicarbazone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // middle_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // left_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // right_n = link_point[2];
                int left_c = mol->addAtom(&atom_c), right_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], left_c);
                add_bond(mol, link_point[0], right_c);
                add_bond(mol, left_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, right_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, left_c, mol->addAtom(&atom_n));
                add_bond(mol, right_c, mol->addAtom(&atom_n));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_boronic_acid()
{
    Atom atom_b("B");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_b = atom_b, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_b);
            // b = h_list[0];
            int o = mol->addAtom(&atom_o);
            int o2 = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], o);
            add_bond(mol, h_list[0], o2);
            add_bond(mol, o, mol->addAtom(&atom_h));
            add_bond(mol, o2, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("boronic acid ()"), 1, func);
}

ModificationType make_carbamate(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbamate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // left_o = link_point[0];
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c);
                add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                add_bond(mol, middle_c, link_point[1]);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_carbamoyl_chloride(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbamoyl chloride (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c);
                add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, middle_c, mol->addAtom(&atom_cl));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_carbazide(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbazide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n1 = link_point[0];
                int left_n2 = mol->addAtom(&atom_n), right_n2 = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], left_n2);
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, left_n2, middle_c);
                add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, middle_c, right_n2);
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n1 = link_point[1];
                add_bond(mol, right_n2, link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_carbene(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    atom_c.setNumRadicalElectrons(2);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_carbodiimide(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbodiimide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c, 2);
                add_bond(mol, link_point[1], middle_c, 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_carbonate_ester(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbonate ester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // up_o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // down_o = link_point[1];
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c);
                add_bond(mol, link_point[1], middle_c);
                add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_carbonyl(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbonyl (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_carboximidate(int h_num = 0)
{
    char str[50];
    sprintf(str, "carboximidate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n = link_point[2];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], link_point[2], 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_carboxylic_acid()
{
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // c = h_list[0];
            int right_o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], right_o);
            add_bond(mol, right_o, mol->addAtom(&atom_h));
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("carboxylic acid ()"), 1, func);
}

ModificationType make_chloroformate()
{
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            // right_o = h_list[0];
            int middle_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], middle_c);
            add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
            add_bond(mol, middle_c, mol->addAtom(&atom_cl));
        }
    };
    return std::make_tuple(std::string("chloroformate ()"), 1, func);
}

ModificationType make_cyanate(bool ocn = true)
{
    char str[50];
    sprintf(str, "cyanate (ocn=%s)", ocn ? "true" : "false");

    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [ocn = ocn, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            if (ocn)
                mol->replaceAtom(h_list[0], &atom_o);
            else
                mol->replaceAtom(h_list[0], &atom_n);
            int c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], c);
            if (ocn)
                add_bond(mol, c, mol->addAtom(&atom_n));
            else
                add_bond(mol, c, mol->addAtom(&atom_o));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_cyanate_ester()
{
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            // o = h_list[0];
            int middle_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], middle_c);
            add_bond(mol, middle_c, mol->addAtom(&atom_n), 3);
        }
    };
    return std::make_tuple(std::string("cyanate ester ()"), 1, func);
}

ModificationType make_cyanimide(int h_num = 0)
{
    char str[50];
    sprintf(str, "cyanimide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], c);
                add_bond(mol, c, mol->addAtom(&atom_n), 3);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_cyanohydrin(int h_num = 0)
{
    char str[50];
    sprintf(str, "cyanohydrin (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int middle_c = mol->addAtom(&atom_c), up_c = mol->addAtom(&atom_c);
                add_bond(mol, middle_c, get_target(mol, link_point[0]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, middle_c, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, middle_c, up_c);
                add_bond(mol, up_c, mol->addAtom(&atom_n), 3);
                int o = mol->addAtom(&atom_o);
                add_bond(mol, middle_c, o);
                add_bond(mol, o, mol->addAtom(&atom_h));
                remove_atoms(mol, link_point.begin(), link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_cyanomethyl()
{
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // right_c = h_list[0];
            int left_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], left_c);
            add_bond(mol, left_c, mol->addAtom(&atom_n), 3);
        }
    };
    return std::make_tuple(std::string("cyanomethyl ()"), 1, func);
}

ModificationType make_cyclopropyl()
{
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // bond_c = h_list[0];
            int left_c = mol->addAtom(&atom_c), right_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], left_c);
            add_bond(mol, h_list[0], right_c);
            add_bond(mol, left_c, right_c);
        }
    };
    return std::make_tuple(std::string("cyclopropyl ()"), 1, func);
}

ModificationType make_diazo(int h_num = 0, bool triple_bond = true)
{
    char str[50];
    sprintf(str, "cyanohydrin (h_num=%d, triple_bond=%s)", h_num, triple_bond ? "true" : "false");

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_n2("N");
    if (triple_bond)
    {
        atom_c.setFormalCharge(-1);
        atom_n.setFormalCharge(1);
    }
    else
    {
        atom_n.setFormalCharge(1);
        atom_n2.setFormalCharge(-1);
    }
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, triple_bond = triple_bond, atom_c = atom_c, atom_n = atom_n, atom_n2 = atom_n2] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], n, triple_bond ? 1 : 2);
                add_bond(mol, n, mol->addAtom(&atom_n2), triple_bond ? 3 : 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_diazonium(const char* halogen)
{
    char str[50];
    sprintf(str, "acyl halide (halogen=%s)", halogen);
    
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_n2("N");
    Atom atom_halogen(halogen);
    atom_halogen.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n, atom_n2 = atom_n2, atom_halogen = atom_halogen] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            // left_n = h_list[0];
            int right_n = mol->addAtom(&atom_n2);
            add_bond(mol, h_list[0], right_n, 3);
            add_bond(mol, right_n, mol->addAtom(&atom_halogen));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_dicarbonate(int h_num = 0)
{
    char str[50];
    sprintf(str, "dicarbonate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // left_o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // right_o = link_point[1];
                int left_c = mol->addAtom(&atom_c), right_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], left_c);
                add_bond(mol, left_c, mol->addAtom(&atom_o), 2);
                int middle_o = mol->addAtom(&atom_o);
                add_bond(mol, left_c, middle_o);
                add_bond(mol, middle_o, right_c);
                add_bond(mol, right_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, right_c, link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_diketopiperazine(int h_num = 0, int isomers = 0) // 0: 2,3; 1: 2,5; 2: 2,6 - isomers
{
    char str[50];
    sprintf(str, "diketopiperazine (h_num=%d, isomers=%s)", h_num, isomers == 0 ? "2,3" : (isomers == 1 ? "2,5" : "2,6"));

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, isomers = isomers, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // up_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // down_n = link_point[1];
                int right_up_c = mol->addAtom(&atom_c), right_down_c = mol->addAtom(&atom_c);
                int left_up_c = mol->addAtom(&atom_c), left_down_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], right_up_c);
                add_bond(mol, right_up_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, right_up_c, right_down_c);
                if (isomers == 0)
                    add_bond(mol, right_down_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, right_down_c, link_point[1]);
                add_bond(mol, link_point[1], left_down_c);
                if (isomers == 1)
                    add_bond(mol, left_down_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, left_down_c, left_up_c);
                if (isomers == 2)
                    add_bond(mol, left_up_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, left_up_c, link_point[0]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_diol()
{
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            // o1 = h_list[0];
            mol->replaceAtom(h_list[1], &atom_o);
            // o2 = h_list[1];
            add_bond(mol, h_list[0], mol->addAtom(&atom_h));
            add_bond(mol, h_list[1], mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("diol ()"), 2, func);
}

ModificationType make_dioxazolone()
{
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // bond_c = h_list[0];
            int n = mol->addAtom(&atom_n);
            add_bond(mol, h_list[0], n, 2);
            int left_o = mol->addAtom(&atom_o), right_o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], left_o);
            add_bond(mol, n, right_o);
            int middle_c = mol->addAtom(&atom_c);
            add_bond(mol, left_o, middle_c);
            add_bond(mol, right_o, middle_c);
            add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("dioxazolone ()"), 1, func);
}

ModificationType make_dioxirane(int h_num = 0)
{
    char str[50];
    sprintf(str, "dioxirane (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int o1 = mol->addAtom(&atom_o), o2 = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], o1);
                add_bond(mol, o1, o2);
                add_bond(mol, o2, link_point[0]);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_diphenyltriazene(int h_num = 0)
{
    char str[50];
    sprintf(str, "diphenyltriazene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                int middle_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], middle_n);
                add_bond(mol, link_point[1], middle_n, 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_disulfide(int h_num = 0)
{
    char str[50];
    sprintf(str, "disulfide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // left_s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // right_s = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_dithiocarbamate(int h_num = 0)
{
    char str[50];
    sprintf(str, "dithiocarbamate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // right_s = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c);
                add_bond(mol, middle_c, mol->addAtom(&atom_s), 2);
                add_bond(mol, middle_c, link_point[1]);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_dithiol()
{
    Atom atom_s("S");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_s = atom_s, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s1 = h_list[0];
            mol->replaceAtom(h_list[1], &atom_s);
            // s2 = h_list[1];
            add_bond(mol, h_list[0], mol->addAtom(&atom_h));
            add_bond(mol, h_list[1], mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("dithiol ()"), 2, func);
}

ModificationType make_enamine(int h_num = 0)
{
    char str[50];
    sprintf(str, "enamine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 5 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // bottom_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // up_c = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n = link_point[2];
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], get_target(mol, link_point[4]));
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 5 - h_num, func);
}

/* enediol */
ModificationType make_vicinal_diol(int h_num = 0)
{
    char str[50];
    sprintf(str, "vicinal diol (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                int left_o = mol->addAtom(&atom_o), right_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], left_o);
                add_bond(mol, left_o, mol->addAtom(&atom_h));
                add_bond(mol, link_point[1], right_o);
                add_bond(mol, right_o, mol->addAtom(&atom_h));
                add_bond(mol, link_point[0], link_point[1], 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_geminal_diol()
{
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // left_c = h_list[0];
            int right_c = mol->addAtom(&atom_c);
            int o1 = mol->addAtom(&atom_o), o2 = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], right_c, 2);
            add_bond(mol, right_c, o1);
            add_bond(mol, o1, mol->addAtom(&atom_h));
            add_bond(mol, right_c, o2);
            add_bond(mol, o2, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("geminal diol ()"), 1, func);
}
/* enediol */

ModificationType make_enediyne(int h_num = 0)
{
    char str[50];
    sprintf(str, "enediyne (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_up_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_up_c = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // left_down_c = link_point[2];
                mol->replaceAtom(link_point[3], &atom_c);
                // right_down_c = link_point[3];
                add_bond(mol, link_point[0], link_point[1], 2);
                int left_middle_c = mol->addAtom(&atom_c), right_middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], left_middle_c);
                add_bond(mol, left_middle_c, link_point[2], 3);
                add_bond(mol, link_point[1], right_middle_c);
                add_bond(mol, right_middle_c, link_point[3], 3);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_enol(int h_num = 0)
{
    char str[50];
    sprintf(str, "enol (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // middle_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], o);
                add_bond(mol, o, mol->addAtom(&atom_h));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_enol_ether(int h_num = 0)
{
    char str[50];
    sprintf(str, "enol ether (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                mol->replaceAtom(link_point[2], &atom_o);
                // o = link_point[2];
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], link_point[2]);
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_enone(int h_num = 0)
{
    char str[50];
    sprintf(str, "enone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // middle_c = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // right_c = link_point[2];
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_episulfide(int h_num = 0)
{
    char str[50];
    sprintf(str, "episulfide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int left_c = mol->addAtom(&atom_c), right_c = mol->addAtom(&atom_c);
                add_bond(mol, left_c, get_target(mol, link_point[0]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, left_c, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, right_c, get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, right_c, get_target(mol, link_point[3]), 1, Bond::BondDir::BEGINWEDGE);
                int s = mol->addAtom(&atom_s);
                add_bond(mol, left_c, s);
                add_bond(mol, right_c, s);
                add_bond(mol, left_c, right_c);
                remove_atoms(mol, link_point.begin(), link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_epoxide(int h_num = 0)
{
    char str[50];
    sprintf(str, "epoxide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int left_c = mol->addAtom(&atom_c), right_c = mol->addAtom(&atom_c);
                add_bond(mol, left_c, get_target(mol, link_point[0]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, left_c, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, right_c, get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, right_c, get_target(mol, link_point[3]), 1, Bond::BondDir::BEGINWEDGE);
                int o = mol->addAtom(&atom_o);
                add_bond(mol, left_c, o);
                add_bond(mol, right_c, o);
                add_bond(mol, left_c, right_c);
                remove_atoms(mol, link_point.begin(), link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_ester(int h_num = 0)
{
    char str[50];
    sprintf(str, "ester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // right_o = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_ether(int h_num = 0)
{
    char str[50];
    sprintf(str, "ether (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_fluorosulfonate()
{
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_f("F");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_s = atom_s, atom_o = atom_o, atom_f = atom_f] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            int s = mol->addAtom(&atom_s);
            add_bond(mol, h_list[0], s);
            add_bond(mol, s, mol->addAtom(&atom_o), 2);
            add_bond(mol, s, mol->addAtom(&atom_o), 2);
            add_bond(mol, s, mol->addAtom(&atom_f));
        }
    };
    return std::make_tuple(std::string("fluorosulfonate ()"), 1, func);
}

ModificationType make_haloalkane(const char* halide)
{
    char str[50];
    sprintf(str, "haloalkane (halide=%s)", halide);
    
    Atom atom_halide(halide);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_halide = atom_halide] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_halide);
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_halohydrin(const char* halide, int h_num = 0)
{
    char str[50];
    sprintf(str, "halohydrin (halide=%s, h_num=%d)", halide, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_h("H");
    Atom atom_halide(halide);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_h = atom_h, atom_halide = atom_halide] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_halide));
                add_bond(mol, link_point[0], link_point[1]);
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[1], o);
                add_bond(mol, o, mol->addAtom(&atom_h));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_haloketone(const char* halide, int h_num = 0)
{
    char str[50];
    sprintf(str, "haloketone (halide=%s, h_num=%d)", halide, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_halide(halide);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_halide = atom_halide] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                add_bond(mol, link_point[1], mol->addAtom(&atom_halide));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_hemithioacetal(int h_num = 0)
{
    char str[50];
    sprintf(str, "hemithioacetal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c, atom_s = atom_s, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                int c = mol->addAtom(&atom_c);
                add_bond(mol, c, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, c, mol->addAtom(&atom_h), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, c, link_point[0]);
                int o = mol->addAtom(&atom_o);
                add_bond(mol, c, o);
                add_bond(mol, o, mol->addAtom(&atom_h));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

/* hydrazides */
ModificationType make_carbohydrazides(int h_num = 0)
{
    char str[50];
    sprintf(str, "carbohydrazides (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // left_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // right_n = link_point[2];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_sulfonohydrazides(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfonohydrazides (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_s("S");
    // atom_s.setFormalCharge(4);
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_s = atom_s, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // left_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // right_n = link_point[2];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_phosphonic_dihydrazides(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphonic dihydrazides (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 7 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    // atom_p.setFormalCharge(2);
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // up_left_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // up_right_n = link_point[2];
                mol->replaceAtom(link_point[3], &atom_n);
                // down_left_n = link_point[3];
                mol->replaceAtom(link_point[4], &atom_n);
                // down_right_n = link_point[4];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], get_target(mol, link_point[5]));
                add_bond(mol, link_point[0], link_point[3]);
                add_bond(mol, link_point[3], link_point[4]);
                add_bond(mol, link_point[4], get_target(mol, link_point[6]));
                remove_atoms(mol, link_point.begin() + 5, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 7 - h_num, func);
}
/* hydrazides */

ModificationType make_hydrazone(int h_num = 0)
{
    char str[50];
    sprintf(str, "hydrazone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_h("H");
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_h = atom_h, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                int bottom_n = mol->addAtom(&atom_n), up_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], bottom_n, 2);
                add_bond(mol, bottom_n, up_n);
                add_bond(mol, up_n, mol->addAtom(&atom_h));
                add_bond(mol, up_n, mol->addAtom(&atom_h));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_hydroperoxide()
{
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            int o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], o);
            add_bond(mol, o, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("hydroperoxide ()"), 1, func);
}

ModificationType make_hydroxamic_acid(int h_num = 0)
{
    char str[50];
    sprintf(str, "hydroxamic acid (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_h("H");
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_h = atom_h, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                int right_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[1], right_o);
                add_bond(mol, right_o, mol->addAtom(&atom_h));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_hydroxylamine(int h_num = 0)
{
    char str[50];
    sprintf(str, "hydroxylamine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_h("H");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_h = atom_h, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], o);
                add_bond(mol, o, mol->addAtom(&atom_h));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_imide(int h_num = 0)
{
    char str[50];
    sprintf(str, "imide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[2]);
                add_bond(mol, link_point[1], link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_imidic_acid(int h_num = 0)
{
    char str[50];
    sprintf(str, "imidic acid (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[0], o);
                add_bond(mol, o, mol->addAtom(&atom_h));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_imidoyl_chloride(int h_num = 0)
{
    char str[50];
    sprintf(str, "imidoyl chloride (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[0], mol->addAtom(&atom_cl));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_imine(int h_num = 0)
{
    char str[50];
    sprintf(str, "imine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_iminium(int h_num = 0)
{
    char str[50];
    sprintf(str, "imiuium (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], link_point[1], 2);
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_isocyanate()
{
    Atom atom_o("O");
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            int c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], c, 2);
            add_bond(mol, c, mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("isocyanate ()"), 1, func);
}

ModificationType make_isocyanide()
{
    Atom atom_c("C");
    atom_c.setFormalCharge(-1);
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            add_bond(mol, h_list[0], mol->addAtom(&atom_c), 3);
        }
    };
    return std::make_tuple(std::string("isocyanide ()"), 1, func);
}

ModificationType make_isodiazene(int h_num = 0)
{
    char str[50];
    sprintf(str, "isodiazene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n_positive("N");
    atom_n_positive.setFormalCharge(1);
    Atom atom_n_negtive("N");
    atom_n_negtive.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n_positive = atom_n_positive, atom_n_negtive = atom_n_negtive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n_positive);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_n_negtive));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_isodiazomethane(int h_num = 0)
{
    char str[50];
    sprintf(str, "isodiazomethane (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_n_positive("N");
    atom_n_positive.setFormalCharge(1);
    Atom atom_c("C");
    atom_c.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_n_positive = atom_n_positive, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int positive_n = mol->addAtom(&atom_n_positive);
                add_bond(mol, link_point[0], positive_n);
                add_bond(mol, positive_n, mol->addAtom(&atom_c), 3);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_isothiocyanate()
{
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_n = atom_n, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            int c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], c, 2);
            add_bond(mol, c, mol->addAtom(&atom_s), 2);
        }
    };
    return std::make_tuple(std::string("isothiocyanate ()"), 1, func);
}

ModificationType make_isothiouronium()
{
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_s("S");
    Atom atom_h("H");
    Atom atom_n_poisitive("N");
    atom_n_poisitive.setFormalCharge(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_n = atom_n, atom_s = atom_s, atom_h = atom_h, atom_n_poisitive = atom_n_poisitive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            int c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], c);
            int up_n = mol->addAtom(&atom_n_poisitive), bottom_n = mol->addAtom(&atom_n);
            add_bond(mol, c, up_n, 2);
            add_bond(mol, up_n, mol->addAtom(&atom_h));
            add_bond(mol, up_n, mol->addAtom(&atom_h));
            add_bond(mol, c, bottom_n);
            add_bond(mol, bottom_n, mol->addAtom(&atom_h));
            add_bond(mol, bottom_n, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("isothiouronium ()"), 1, func);
}

ModificationType make_ketene(int h_num = 0)
{
    char str[50];
    sprintf(str, "ketene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c, 2);
                add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_ketenimine(int h_num = 0)
{
    char str[50];
    sprintf(str, "ketenimine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // bottom_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c, 2);
                add_bond(mol, middle_c, link_point[1], 2);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_ketone(int h_num = 0)
{
    char str[50];
    sprintf(str, "ketone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_ketyl(int h_num = 0)
{
    char str[50];
    sprintf(str, "ketyl (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1]);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_alpha_lactone(const char* atom_name, int h_num = 0)
{
    char str[50];
    sprintf(str, "alpha lactone (atom_name=%s, h_num=%d)", atom_name, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom(atom_name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int c2 = mol->addAtom(&atom_c), atom_idx = mol->addAtom(&atom);
                add_bond(mol, link_point[0], c2);
                add_bond(mol, c2, mol->addAtom(&atom_o), 2);
                add_bond(mol, c2, atom_idx);
                add_bond(mol, atom_idx, link_point[0]);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_beta_lactone(const char* atom_name, int h_num = 0)
{
    char str[50];
    sprintf(str, "beta lactone (atom_name=%s, h_num=%d)", atom_name, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom(atom_name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                for (int i = 0; i < 2; ++i)
                    mol->replaceAtom(link_point[i], &atom_c);
                // c1 = link_point[0];
                // c2 = link_point[1];
                for (int i = 0; i < 2; ++i)
                    add_bond(mol, link_point[i], get_target(mol, link_point[i + 2]));
                int c3 = mol->addAtom(&atom_c), atom_idx = mol->addAtom(&atom);
                add_bond(mol, c3, mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], c3);
                add_bond(mol, c3, atom_idx);
                add_bond(mol, atom_idx, link_point[0]);
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_gamma_lactone(const char* atom_name, int h_num = 0)
{
    char str[50];
    sprintf(str, "gamma lactone (atom_name=%s, h_num=%d)", atom_name, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 6 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom(atom_name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                for (int i = 0; i < 3; ++i)
                    mol->replaceAtom(link_point[i], &atom_c);
                // c1 = link_point[0];
                // c2 = link_point[1];
                // c3 = link_point[2];
                for (int i = 0; i < 3; ++i)
                    add_bond(mol, link_point[i], get_target(mol, link_point[i + 3]));
                int c4 = mol->addAtom(&atom_c), atom_idx = mol->addAtom(&atom);
                add_bond(mol, c4, mol->addAtom(&atom_o), 2);
                for (int i = 0; i < 2; ++i)
                    add_bond(mol, link_point[i], link_point[i + 1]);
                add_bond(mol, link_point[2], c4);
                add_bond(mol, c4, atom_idx);
                add_bond(mol, atom_idx, link_point[0]);
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 6 - h_num, func);
}

ModificationType make_delta_lactone(const char* atom_name, int h_num = 0)
{
    char str[50];
    sprintf(str, "delta lactone (atom_name=%s, h_num=%d)", atom_name, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 8 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom(atom_name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                for (int i = 0; i < 4; ++i)
                    mol->replaceAtom(link_point[i], &atom_c);
                // c1 = link_point[0];
                // c2 = link_point[1];
                // c3 = link_point[2];
                // c4 = link_point[3];
                for (int i = 0; i < 4; ++i)
                    add_bond(mol, link_point[i], get_target(mol, link_point[i + 4]));
                int c5 = mol->addAtom(&atom_c), atom_idx = mol->addAtom(&atom);
                add_bond(mol, c5, mol->addAtom(&atom_o), 2);
                for (int i = 0; i < 3; ++i)
                    add_bond(mol, link_point[i], link_point[i + 1]);
                add_bond(mol, link_point[3], c5);
                add_bond(mol, c5, atom_idx);
                add_bond(mol, atom_idx, link_point[0]);
                remove_atoms(mol, link_point.begin() + 4, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 8 - h_num, func);
}

ModificationType make_epislon_lactone(const char* atom_name, int h_num = 0)
{
    char str[50];
    sprintf(str, "epsilon lactone (atom_name=%s, h_num=%d)", atom_name, h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 10 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom(atom_name);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                for (int i = 0; i < 5; ++i)
                    mol->replaceAtom(link_point[i], &atom_c);
                // c1 = link_point[0];
                // c2 = link_point[1];
                // c3 = link_point[2];
                // c4 = link_point[3];
                // c5 = link_point[4];
                for (int i = 0; i < 5; ++i)
                    add_bond(mol, link_point[i], get_target(mol, link_point[i + 5]));
                int c6 = mol->addAtom(&atom_c), atom_idx = mol->addAtom(&atom);
                add_bond(mol, c6, mol->addAtom(&atom_o), 2);
                for (int i = 0; i < 4; ++i)
                    add_bond(mol, link_point[i], link_point[i + 1]);
                add_bond(mol, link_point[4], c6);
                add_bond(mol, c6, atom_idx);
                add_bond(mol, atom_idx, link_point[0]);
                remove_atoms(mol, link_point.begin() + 5, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 10 - h_num, func);
}

ModificationType make_methanedithiol(int h_num = 0)
{
    char str[50];
    sprintf(str, "methanedithiol (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                mol->replaceAtom(link_point[1], &atom_s);
                // right_up_s = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // center_c = link_point[2];
                add_bond(mol, link_point[2], link_point[0]);
                add_bond(mol, link_point[2], link_point[1]);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_methine(int h_num = 0)
{
    char str[50];
    sprintf(str, "methine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_h("H");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_h = atom_h, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int target_idx = get_target(mol, link_point[1]);
                if (target_idx == get_target(mol, link_point[2]))
                {
                    mol->replaceAtom(link_point[0], &atom_c);
                    // c = link_point[0];
                    add_bond(mol, link_point[0], get_target(mol, link_point[1]), 2);
                    add_bond(mol, link_point[0], mol->addAtom(&atom_h));
                    remove_atoms(mol, link_point.begin() + 1, link_point.end());
                }
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_methylene_bridge()
{
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], get_target(mol, h_list[1]));
            mol->removeAtom(h_list[1]);
        }
    };
    return std::make_tuple(std::string("methylene bridge ()"), 2, func);
}

ModificationType make_methylenedioxy(int h_num = 0)
{
    char str[50];
    sprintf(str, "methylenedioxy (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // up_o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // bottom_o = link_point[1];
                int c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], c);
                add_bond(mol, link_point[1], c);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

// ModificationType make_methylidene()
// {
//     Atom atom_c("C");
//     std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
//         [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
//         {
//             int c = mol->addAtom(&atom_c);
//             add_bond(mol, c, get_target(mol, h_list[0]), 2);
//             mol->removeAtom(h_list[0]);
//         }
//     };
//     return std::make_tuple(std::string("methylidene ()"), 1, func);
// }

ModificationType make_methylidene()
{
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            int target_idx = get_target(mol, h_list[1]);
            if (get_target(mol, h_list[0]) == target_idx)
            {
                mol->replaceAtom(h_list[0], &atom_c);
                add_bond(mol, h_list[0], target_idx);
                mol->removeAtom(h_list[1]);
            }
        }
    };
    return std::make_tuple(std::string("methylidene ()"), 2, func);
}

ModificationType make_nitrate_ester()
{
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_o("O");
    Atom atom_o_negtive("O");
    atom_o_negtive.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n, atom_o = atom_o, atom_o_negtive = atom_o_negtive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            int n = mol->addAtom(&atom_n);
            add_bond(mol, h_list[0], n);
            add_bond(mol, n, mol->addAtom(&atom_o), 2);
            add_bond(mol, n, mol->addAtom(&atom_o_negtive));
        }
    };
    return std::make_tuple(std::string("nitrate ester ()"), 1, func);
}

ModificationType make_nitrene()
{
    Atom atom_n("N");
    // atom_n.setNumRadicalElectrons(4);
    // atom_n.setFormalCharge(2);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
        }
    };
    return std::make_tuple(std::string("nitrene ()"), 1, func);
}

ModificationType make_nitrile_ylide(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitrile ylide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_c("C");
    Atom atom_c_negtive("C");
    atom_c_negtive.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_c = atom_c, atom_c_negtive = atom_c_negtive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c_negtive);
                // right_c = link_point[1];
                int n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], n, 3);
                add_bond(mol, n, link_point[1]);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_nitrilimine(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitrilimine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n_positive("N");
    atom_n_positive.setFormalCharge(1);
    Atom atom_n_negtive("N");
    atom_n_negtive.setFormalCharge(-1);
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n_positive = atom_n_positive, atom_n_negtive = atom_n_negtive, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n_negtive);
                // right_n = link_point[1];
                int left_n = mol->addAtom(&atom_n_positive);
                add_bond(mol, link_point[0], left_n, 3);
                add_bond(mol, left_n, link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_nitro()
{
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_o("O");
    Atom atom_o_negtive("O");
    atom_o_negtive.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n, atom_o = atom_o, atom_o_negtive = atom_o_negtive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            // n = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o_negtive));
        }
    };
    return std::make_tuple(std::string("nitro ()"), 1, func);
}

ModificationType make_nitroalkene(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitroalkene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_o_negtive("O");
    atom_o_negtive.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_c = atom_c, atom_o = atom_o, atom_o_negtive = atom_o_negtive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                int n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], n);
                add_bond(mol, n, mol->addAtom(&atom_o), 2);
                add_bond(mol, n, mol->addAtom(&atom_o_negtive));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_nitroamine(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitroamine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_n_positive("N");
    atom_n_positive.setFormalCharge(1);
    Atom atom_o("O");
    Atom atom_o_negtive("O");
    atom_o_negtive.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_n_positive = atom_n_positive, atom_o = atom_o, atom_o_negtive = atom_o_negtive] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n = link_point[0];
                int positive_n = mol->addAtom(&atom_n_positive);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], positive_n);
                add_bond(mol, positive_n, mol->addAtom(&atom_o), 2);
                add_bond(mol, positive_n, mol->addAtom(&atom_o_negtive));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_nitrolic_acid()
{
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_n_positive("N");
    atom_n_positive.setFormalCharge(1);
    Atom atom_o("O");
    Atom atom_o_negtive("O");
    atom_o_negtive.setFormalCharge(-1);
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_n = atom_n, atom_n_positive = atom_n_positive, atom_o = atom_o, atom_o_negtive = atom_o_negtive, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // c = h_list[0];
            int n = mol->addAtom(&atom_n), positive_n = mol->addAtom(&atom_n_positive), o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], n, 2);
            add_bond(mol, n, o);
            add_bond(mol, o, mol->addAtom(&atom_h));
            add_bond(mol, h_list[0], positive_n);
            add_bond(mol, positive_n, mol->addAtom(&atom_o), 2);
            add_bond(mol, positive_n, mol->addAtom(&atom_o_negtive));
        }
    };
    return std::make_tuple(std::string("nitrolic acid ()"), 1, func);
}

ModificationType make_nitronate(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitronate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_o("O");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // up_o = link_point[1];
                mol->replaceAtom(link_point[2], &atom_o);
                // bottom_o = link_point[2];
                int n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], n, 2);
                add_bond(mol, n, link_point[1]);
                add_bond(mol, n, link_point[2]);
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_nitrone(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitrone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_n("N");
    atom_n.setFormalCharge(1);
    Atom atom_o("O");
    atom_o.setFormalCharge(-1);
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_o = atom_o, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_nitrosamine(int h_num = 0)
{
    char str[50];
    sprintf(str, "nitrosamine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // up_n = link_point[0];
                int bottom_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], bottom_n);
                add_bond(mol, bottom_n, mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_nitroso()
{
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_n);
            // n = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("nitroso ()"), 1, func);
}

ModificationType make_s_nitrosothiol()
{
    Atom atom_n("N");
    Atom atom_o("O");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_n = atom_n, atom_o = atom_o, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            int n = mol->addAtom(&atom_n);
            add_bond(mol, h_list[0], n);
            add_bond(mol, n, mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("S-Nitrosothiol ()"), 1, func);
}

ModificationType make_organic_acid_anhydride(int h_num = 0)
{
    char str[50];
    sprintf(str, "organic acid anhydride (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                int middle_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], middle_o);
                add_bond(mol, link_point[1], middle_o);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_organic_peroxide(int h_num = 0)
{
    char str[50];
    sprintf(str, "organic peroxide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // left_o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // right_o = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_orthoester(int h_num = 0)
{
    char str[50];
    sprintf(str, "orthoester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // up_o = link_point[1];
                mol->replaceAtom(link_point[2], &atom_o);
                // right_o = link_point[2];
                mol->replaceAtom(link_point[3], &atom_o);
                // bottom_o = link_point[3];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], link_point[2]);
                add_bond(mol, link_point[0], link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_oxaziridine(int h_num = 0)
{
    char str[50];
    sprintf(str, "oxaziridine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // c = link_point[1];
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], o);
                add_bond(mol, link_point[1], o);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphaalkene(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphaalkene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // c = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphaalkyne()
{
    Atom atom_c("C");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_p), 3);
        }
    };
    return std::make_tuple(std::string("phosphaalkyne ()"), 1, func);
}

ModificationType make_phosphate(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // wedge_o = link_point[1];
                mol->replaceAtom(link_point[2], &atom_o);
                // dash_o = link_point[2];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p, 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, link_point[2], p, 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphinate(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphinate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o = link_point[0];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, p, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, p, get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphine(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p= atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINDASH);
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphine_imide(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphine imide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_p("P");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_phosphine_oxide(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphine oxide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_p("P");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphinite(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphinite (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_p("P");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o = link_point[0];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, p, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, p, get_target(mol, link_point[2]), 1, Bond::BondDir::BEGINWEDGE);
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphinous_acids(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphinous acids (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_p("P");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int p = mol->addAtom(&atom_p);
                add_bond(mol, p, get_target(mol, link_point[0]), 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, p, get_target(mol, link_point[1]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
                add_bond(mol, p, mol->addAtom(&atom_h));
                remove_atoms(mol, link_point.begin(), link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_phosphite_ester(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphite ester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_p("P");
    atom_p.setNumRadicalElectrons(2);
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                mol->replaceAtom(link_point[1], &atom_o);
                mol->replaceAtom(link_point[2], &atom_o);              
                int p = mol->addAtom(&atom_p);
                add_bond(mol, p, link_point[0]);
                add_bond(mol, p, link_point[1]);
                add_bond(mol, p, link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphonate(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphonate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_p("P");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                mol->replaceAtom(link_point[1], &atom_o);
                mol->replaceAtom(link_point[2], &atom_p);
                add_bond(mol, link_point[2], link_point[0]);
                add_bond(mol, link_point[2], link_point[1]);
                add_bond(mol, link_point[2], mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphonite(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphonite (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_p("P");
    atom_p.setNumRadicalElectrons(2);
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                mol->replaceAtom(link_point[1], &atom_o);
                mol->replaceAtom(link_point[2], &atom_p);
                add_bond(mol, link_point[2], link_point[0]);
                add_bond(mol, link_point[2], link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_phosphonium(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphonium (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_p("P");
    atom_p.setFormalCharge(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

/* Phosphoramidate */
ModificationType make_phosphorodiamidate(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphorodiamidate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 5 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n1 = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n2 = link_point[2];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p);
                add_bond(mol, link_point[2], p);
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                add_bond(mol, link_point[2], get_target(mol, link_point[4]));
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 5 - h_num, func);
}

ModificationType make_phosphoramidate(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphoramidate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n = link_point[2];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p);
                add_bond(mol, link_point[2], p);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}
/* Phosphoramidate */

ModificationType make_phosphoramides(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphoramides (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 6 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n2 = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n3 = link_point[2];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p);
                add_bond(mol, link_point[2], p);
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                add_bond(mol, link_point[1], get_target(mol, link_point[4]));
                add_bond(mol, link_point[2], get_target(mol, link_point[5]));
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 6 - h_num, func);
}

ModificationType make_phosphoramidite(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphoramidite (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    atom_p.setNumRadicalElectrons(2);
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // n = link_point[2];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p);
                add_bond(mol, link_point[2], p);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_phosphorane(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphorane (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 5 - h_num, 1);
    
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_p);
                // p = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], get_target(mol, link_point[3]), 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, link_point[0], get_target(mol, link_point[4]), 1, Bond::BondDir::BEGINDASH);
                remove_atoms(mol, link_point.begin() + 1, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 5 - h_num, func);
}

ModificationType make_phosphorochloridate(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphorochloridate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p);
                add_bond(mol, p, mol->addAtom(&atom_cl));
                add_bond(mol, p, mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

/* Phosphorochloridite */
ModificationType make_phosphochloridite()
{
    Atom atom_o("O");
    Atom atom_p("P");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_p = atom_p, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            // o = h_list[0];
            int p = mol->addAtom(&atom_p);
            add_bond(mol, h_list[0], p);
            add_bond(mol, p, mol->addAtom(&atom_cl));
            add_bond(mol, p, mol->addAtom(&atom_cl));
        }
    };
    return std::make_tuple(std::string("phosphochloridite ()"), 1, func);
}

ModificationType make_phosphodichloridite(int h_num = 0)
{
    char str[50];
    sprintf(str, "phosphodichloridite (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_p("P");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_p = atom_p, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, link_point[0], p);
                add_bond(mol, link_point[1], p);
                add_bond(mol, p, mol->addAtom(&atom_cl));
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}
/* Phosphorochloridite */

ModificationType make_phosphoryl()
{
    Atom atom_o("O");
    Atom atom_o_negtive("O");
    atom_o_negtive.setFormalCharge(-1);
    Atom atom_p("P");
    atom_p.setFormalCharge(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_o_negtive = atom_o_negtive, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_p);
            // p = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_o_negtive));
            add_bond(mol, h_list[0], mol->addAtom(&atom_o_negtive));
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
        }
    };
    return std::make_tuple(std::string("phosphoryl ()"), 1, func);
}

ModificationType make_propenyl()
{
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // c = h_list[0];
            int middle_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], middle_c, 2);
            add_bond(mol, middle_c, mol->addAtom(&atom_c));
        }
    };
    return std::make_tuple(std::string("propenyl ()"), 1, func);
}

/* Quinone methide */
ModificationType make_para_quinone_methide()
{
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            int target = get_target(mol, h_list[0]);
            if (target == get_target(mol, h_list[1]))
            {
                int ori_c = mol->addAtom(&atom_c);
                add_bond(mol, ori_c, target, 2);
                std::vector<int> c;
                for (int i = 0; i < 5; ++i)
                    c.emplace_back(mol->addAtom(&atom_c));
                add_bond(mol, ori_c, c[0]);
                add_bond(mol, c[0], c[1], 2);
                add_bond(mol, c[1], c[2]);
                add_bond(mol, c[2], mol->addAtom(&atom_o), 2);
                add_bond(mol, c[2], c[3]);
                add_bond(mol, c[3], c[4], 2);
                add_bond(mol, c[4], ori_c);
                remove_atoms(mol, h_list.begin(), h_list.end());
            }
        }
    };
    return std::make_tuple(std::string("para quinone methide ()"), 2, func);
}

ModificationType make_ortho_quinone_methide()
{
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            int target = get_target(mol, h_list[0]);
            if (target == get_target(mol, h_list[1]))
            {
                int ori_c = mol->addAtom(&atom_c);
                add_bond(mol, ori_c, target, 2);
                std::vector<int> c;
                for (int i = 0; i < 5; ++i)
                    c.emplace_back(mol->addAtom(&atom_c));
                add_bond(mol, ori_c, c[0]);
                add_bond(mol, c[0], c[1], 2);
                add_bond(mol, c[1], c[2]);
                add_bond(mol, c[2], c[3], 2);
                add_bond(mol, c[3], c[4]);
                add_bond(mol, c[4], mol->addAtom(&atom_o), 2);
                add_bond(mol, c[4], ori_c);
                remove_atoms(mol, h_list.begin(), h_list.end());
            }
        }
    };
    return std::make_tuple(std::string("ortho quinone methide ()"), 2, func);
}
/* Quinone methide */

ModificationType make_reductone(int h_num = 0)
{
    char str[50];
    sprintf(str, "reductone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                int middle_c = mol->addAtom(&atom_c);
                int left_o = mol->addAtom(&atom_o), middle_o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], left_o);
                add_bond(mol, left_o, mol->addAtom(&atom_h));
                add_bond(mol, link_point[0], middle_c, 2);
                add_bond(mol, middle_c, middle_o);
                add_bond(mol, middle_o, mol->addAtom(&atom_h));
                add_bond(mol, middle_c, link_point[1]);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_schiff_base(int h_num = 0)
{
    char str[50];
    sprintf(str, "schiff base (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1], 2);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_selenenic_acid()
{
    Atom atom_se("Se");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_se = atom_se, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_se);
            // se = h_list[0];
            int o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], o);
            add_bond(mol, o, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("selenenic acid ()"), 1, func);
}

ModificationType make_selenol()
{
    Atom atom_se("Se");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_se = atom_se, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_se);
            // se = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("selenol ()"), 1, func);
}

ModificationType make_selenonic_acid()
{
    Atom atom_se("Se");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_se = atom_se, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_se);
            // se = h_list[0];
            int o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], o);
            add_bond(mol, o, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("selenonic acid ()"), 1, func);
}

ModificationType make_selone(int h_num = 0)
{
    char str[50];
    sprintf(str, "selone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_se("Se");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_se = atom_se] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_se), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_semicarbazide(int h_num = 0)
{
    char str[50];
    sprintf(str, "semicarbazide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 5 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // middle_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // right_n = link_point[2];
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                int c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], c);
                add_bond(mol, c, mol->addAtom(&atom_o), 2);
                add_bond(mol, c, link_point[1]);
                add_bond(mol, link_point[1], link_point[2]);
                add_bond(mol, link_point[2], get_target(mol, link_point[4]));
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 5 - h_num, func);
}

ModificationType make_semicarbazone(int h_num = 0)
{
    char str[50];
    sprintf(str, "semicarbazone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 5 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // middle_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // right_n = link_point[2];
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                int left_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], left_n, 2);
                add_bond(mol, left_n, link_point[1]);
                int right_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[1], right_c);
                add_bond(mol, right_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, right_c, link_point[2]);
                add_bond(mol, link_point[2], get_target(mol, link_point[4]));
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 5 - h_num, func);
}

ModificationType make_silyl_enol_ether(int h_num = 0)
{
    char str[50];
    sprintf(str, "silyl enol ether (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 6 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_si("Si");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_si = atom_si] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_si);
                // si = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // left_c = link_point[2];
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                add_bond(mol, link_point[0], get_target(mol, link_point[4]));
                int o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], o);
                add_bond(mol, o, link_point[1]);
                add_bond(mol, link_point[1], link_point[2], 2);
                add_bond(mol, link_point[2], get_target(mol, link_point[5]));
                remove_atoms(mol, link_point.begin() + 3, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 6 - h_num, func);
}

ModificationType make_silyl_ether(int h_num = 0)
{
    char str[50];
    sprintf(str, "silyl ether (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_si("Si");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_o = atom_o, atom_si = atom_si] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_si);
                // si = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], get_target(mol, link_point[3]));
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_sulfamoyl_fluoride(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfamoyl fluoride (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_f("F");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_f = atom_f, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int s = mol->addAtom(&atom_s);
                add_bond(mol, link_point[0], s);
                add_bond(mol, s, mol->addAtom(&atom_o), 2);
                add_bond(mol, s, mol->addAtom(&atom_o), 2);
                add_bond(mol, s, mol->addAtom(&atom_f));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfenamide(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfenamide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // s = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1]);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_sulfenic_acid()
{
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_s = atom_s, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            int o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], o);
            add_bond(mol, o, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("Sulfenic acid ()"), 1, func);
}

ModificationType make_sulfenyl_chloride()
{
    Atom atom_s("S");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_s = atom_s, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_cl));
        }
    };
    return std::make_tuple(std::string("sulfenyl chloride ()"), 1, func);
}

ModificationType make_sulfide(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfilimine(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfilimine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1], 2);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_sulfinamide(bool double_bond = true, int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfinamide (double_bond=%s, h_num=%d)", double_bond ? "true" : "false", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_s("S");
    if (!double_bond)
        atom_s.setFormalCharge(1);
    Atom atom_n("N");
    Atom atom_o("O");
    if (!double_bond)
        atom_o.setFormalCharge(-1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, double_bond = double_bond, atom_s = atom_s, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), double_bond ? 2 : 1);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_sulfinic_acid()
{
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_s = atom_s, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            int right_o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], right_o);
            add_bond(mol, right_o, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("sulfinic acid ()"), 1, func);
}

ModificationType make_sulfite_ester(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfite ester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                int s = mol->addAtom(&atom_s);
                add_bond(mol, s, link_point[0], 1, Bond::BondDir::BEGINWEDGE);
                add_bond(mol, s, link_point[1], 1, Bond::BondDir::BEGINDASH);
                add_bond(mol, s, mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfonamide(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfonamide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_sulfonanilide(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfonanilide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
                int begin_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[1], begin_c);
                int current_c, last_c = begin_c;
                for (int i = 0; i < 5; ++i)
                {
                    current_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, current_c, Bond::BondType::AROMATIC);
                    last_c = current_c;
                }
                mol->addBond(last_c, begin_c, Bond::BondType::AROMATIC);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfonate(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfonate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // up_o = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfone(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfonic_acid()
{
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_s = atom_s, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            int right_o = mol->addAtom(&atom_o);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], right_o);
            add_bond(mol, right_o, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("sulfonic acid ()"), 1, func);
}

ModificationType make_sulfonyl_halide(const char* halide)
{
    char str[50];
    sprintf(str, "sulfonyl halide (halide=%s)", halide);
    
    Atom atom_o("O");
    Atom atom_s("S");
    Atom atom_halide(halide);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_o = atom_o, atom_s = atom_s, atom_halide = atom_halide] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_halide));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_sulfoxide(int h_num = 0)
{
    char str[50];
    sprintf(str, "sulfoxide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                // s = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_telluroketone(int h_num = 0)
{
    char str[50];
    sprintf(str, "telluroketone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_te("Te");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_te = atom_te] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_te), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_tellurol()
{
    Atom atom_te("Te");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_te = atom_te, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_te);
            // te = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("tellurol ()"), 1, func);
}

ModificationType make_thiadiazoles(int h_num = 0, int isomers = 0) // 0: 1,2,3; 1: 1,2,4; 2: 1,2,5; 3: 1,3,4 - isomers
{
    char str[50];
    sprintf(str, "thiadiazoles (h_num=%d, isomers=%s)", h_num, isomers == 0 ? "1,2,3" : (isomers == 1 ? "1,2,4" : (isomers == 2 ? "1,2,5" : "1,3,4")));

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, isomers = isomers, atom_c = atom_c, atom_n = atom_n, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                int s = mol->addAtom(&atom_s);
                int current_atom, last_atom = s, link_point_idx = 0;
                for (int i = 2; i < 6; ++i)
                {
                    if ((i == 2 && (isomers == 0 || isomers == 1 || isomers == 2)) ||
                        (i == 3 && (isomers == 0 || isomers == 3)) ||
                        (i == 4 && (isomers == 1 || isomers == 3)) ||
                        (i == 5 && isomers == 2))
                        current_atom = mol->addAtom(&atom_n);
                    else
                    {
                        mol->replaceAtom(link_point[link_point_idx], &atom_c);
                        current_atom = link_point[link_point_idx++];
                    }
                    add_bond(mol, last_atom, current_atom, (i == 3 || i == 5) ? 2 : 1);
                    last_atom = current_atom;
                }
                add_bond(mol, last_atom, s);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_thial()
{
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_s = atom_s, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // c = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_s), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string("thial ()"), 1, func);
}

/* Thioacetal */
ModificationType make_thioacetal(int h_num = 0)
{
    char str[50];
    sprintf(str, "thioacetal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o = link_point[1];
                mol->replaceAtom(link_point[2], &atom_s);
                // s = link_point[2];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_dithioacetal(int h_num = 0)
{
    char str[50];
    sprintf(str, "dithioacetal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // s1 = link_point[1];
                mol->replaceAtom(link_point[2], &atom_s);
                // s2 = link_point[2];
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[0], link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}
/* Thioacetal */

ModificationType make_thioacyl_chloride()
{
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_cl("Cl");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_s = atom_s, atom_cl = atom_cl] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // c = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_s), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_cl));
        }
    };
    return std::make_tuple(std::string("thioacyl chloride ()"), 1, func);
}

ModificationType make_thioamide(int h_num = 0)
{
    char str[50];
    sprintf(str, "thioamide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_s), 2);
                add_bond(mol, link_point[0], link_point[1]);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_thiocarbamate(bool o_isomers = true, int h_num = 0)
{
    char str[50];
    sprintf(str, "%s-thiocarbamate (h_num=%d)", o_isomers ? "O" : "S", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, o_isomers = o_isomers, atom_c = atom_c, atom_s = atom_s, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], o_isomers ? &atom_o : &atom_s);
                // o or s = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // n = link_point[1];
                int c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], c);
                add_bond(mol, c, mol->addAtom(o_isomers ? &atom_s : &atom_o), 2);
                add_bond(mol, c, link_point[1]);
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_thiocarboxylic_acid(bool thione_form = true)
{
    char str[50];
    sprintf(str, "thiocarboxylic acid (thione_form=%s)", thione_form ? "true" : "false");

    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [thione_form = thione_form, atom_c = atom_c, atom_s = atom_s, atom_o = atom_o, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // c = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(thione_form ? &atom_s : &atom_o), 2);
            int atom = mol->addAtom(thione_form ? &atom_o : &atom_s);
            add_bond(mol, h_list[0], atom);
            add_bond(mol, atom, mol->addAtom(&atom_h));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_thiocyanate()
{
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_s = atom_s, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            int c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], c);
            add_bond(mol, c, mol->addAtom(&atom_n), 3);
        }
    };
    return std::make_tuple(std::string("thiocyanate ()"), 1, func);
}

ModificationType make_thioester(int h_num = 0)
{
    char str[50];
    sprintf(str, "thioester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // s = link_point[1];
                add_bond(mol, link_point[0], mol->addAtom(&atom_o), 2);
                add_bond(mol, link_point[0], link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_thioketal(int h_num = 0)
{
    char str[50];
    sprintf(str, "thioketal (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                mol->replaceAtom(link_point[1], &atom_s);
                // right_up_s = link_point[1];
                mol->replaceAtom(link_point[2], &atom_c);
                // center_c = link_point[2];
                add_bond(mol, link_point[2], link_point[0]);
                add_bond(mol, link_point[2], link_point[1]);
                add_bond(mol, link_point[2], get_target(mol, link_point[3]));
                mol->removeAtom(link_point[3]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_thioketene(int h_num = 0)
{
    char str[50];
    sprintf(str, "thioketene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], c, 2);
                add_bond(mol, c, mol->addAtom(&atom_s), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_thioketone(int h_num = 0)
{
    char str[50];
    sprintf(str, "thioketone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                add_bond(mol, link_point[0], mol->addAtom(&atom_s), 2);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_thiol(int h_num = 0)
{
    char str[50];
    sprintf(str, "thiol (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_s);
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_thiophosphate(int h_num = 0)
{
    char str[50];
    sprintf(str, "thiophosphate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_p("P");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_s = atom_s, atom_o = atom_o, atom_p = atom_p] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o1 = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o2 = link_point[1];
                mol->replaceAtom(link_point[2], &atom_o);
                // o3 = link_point[2];
                int p = mol->addAtom(&atom_p);
                add_bond(mol, p, mol->addAtom(&atom_s), 2);
                add_bond(mol, p, link_point[0]);
                add_bond(mol, p, link_point[1]);
                add_bond(mol, p, link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_thiourea(bool thione_form = true, int h_num = 0)
{
    char str[50];
    sprintf(str, "thiourea (thione_form=%s, h_num=%d)", thione_form ? "true" : "false", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_s("S");
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, thione_form = thione_form, atom_s = atom_s, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                if (!thione_form)
                    mol->replaceAtom(link_point[2], &atom_s);
                    // s = link_point[2];
                int c = mol->addAtom(&atom_c);
                if (thione_form)
                    add_bond(mol, c, mol->addAtom(&atom_s), 2);
                else
                    add_bond(mol, c, link_point[2]);
                add_bond(mol, link_point[0], get_target(mol, link_point[thione_form ? 2 : 3]));
                if (thione_form)
                    add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                add_bond(mol, c, link_point[0]);
                add_bond(mol, c, link_point[1], thione_form ? 1 : 2);
                remove_atoms(mol, link_point.begin() + (thione_form ? 2 : 3), link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

/* Tosylhydrazone */
ModificationType make_tosyl()
{
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_s);
            // s = h_list[0];
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            add_bond(mol, h_list[0], mol->addAtom(&atom_o), 2);
            int begin_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], begin_c);
            int current_c, last_c = begin_c;
            for (int i = 0; i < 5; ++i)
            {
                current_c = mol->addAtom(&atom_c);
                mol->addBond(last_c, current_c, Bond::BondType::AROMATIC);
                if (i == 2)
                    add_bond(mol, current_c, mol->addAtom(&atom_c));
                last_c = current_c;
            }
            mol->addBond(last_c, begin_c, Bond::BondType::AROMATIC);
        }
    };
    return std::make_tuple(std::string("tosyl ()"), 1, func);
}

ModificationType make_tosylate()
{
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_s = atom_s, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            // o = h_list[0];
            int s = mol->addAtom(&atom_s);
            add_bond(mol, s, h_list[0]);
            add_bond(mol, s, mol->addAtom(&atom_o), 2);
            add_bond(mol, s, mol->addAtom(&atom_o), 2);
            int begin_c = mol->addAtom(&atom_c);
            add_bond(mol, s, begin_c);
            int current_c, last_c = begin_c;
            for (int i = 0; i < 5; ++i)
            {
                current_c = mol->addAtom(&atom_c);
                mol->addBond(last_c, current_c, Bond::BondType::AROMATIC);
                if (i == 2)
                    add_bond(mol, current_c, mol->addAtom(&atom_c));
                last_c = current_c;
            }
            mol->addBond(last_c, begin_c, Bond::BondType::AROMATIC);
        }
    };
    return std::make_tuple(std::string("tosylate ()"), 1, func);
}

ModificationType make_tosylhydrazone(int h_num = 0)
{
    char str[50];
    sprintf(str, "tosylhydrazone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    Atom atom_o("O");
    Atom atom_n("N");
    Atom atom_h("H");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s, atom_o = atom_o, atom_n = atom_n, atom_h = atom_h] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                add_bond(mol, link_point[0], get_target(mol, link_point[1]));
                int left_n = mol->addAtom(&atom_n), right_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], right_n, 2);
                add_bond(mol, right_n, left_n);
                add_bond(mol, left_n, mol->addAtom(&atom_h));
                int s = mol->addAtom(&atom_s);
                add_bond(mol, s, left_n);
                add_bond(mol, s, mol->addAtom(&atom_o), 2);
                add_bond(mol, s, mol->addAtom(&atom_o), 2);
                int begin_c = mol->addAtom(&atom_c);
                add_bond(mol, s, begin_c);
                int current_c, last_c = begin_c;
                for (int i = 0; i < 5; ++i)
                {
                    current_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, current_c, Bond::BondType::AROMATIC);
                    if (i == 2)
                        add_bond(mol, current_c, mol->addAtom(&atom_c));
                    last_c = current_c;
                }
                mol->addBond(last_c, begin_c, Bond::BondType::AROMATIC);
                mol->removeAtom(link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}
/* Tosylhydrazone */

ModificationType make_triazenes(int h_num = 0)
{
    char str[50];
    sprintf(str, "triazenes (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);
    
    Atom atom_n("N");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_n = atom_n] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                int middle_n = mol->addAtom(&atom_n);
                add_bond(mol, link_point[0], middle_n);
                add_bond(mol, middle_n, link_point[1], 2);
                mol->removeAtom(link_point[2]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 3 - h_num, func);
}

ModificationType make_triuret(int h_num = 0)
{
    char str[50];
    sprintf(str, "triuret (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 6 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_up_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // right_up_n = link_point[1];
                mol->replaceAtom(link_point[2], &atom_n);
                // left_bottom_n = link_point[2];
                mol->replaceAtom(link_point[3], &atom_n);
                // right_bottom_n = link_point[3];
                add_bond(mol, link_point[0], get_target(mol, link_point[4]));
                add_bond(mol, link_point[1], get_target(mol, link_point[5]));
                int left_c = mol->addAtom(&atom_c), right_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], left_c);
                add_bond(mol, left_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, left_c, link_point[2]);
                add_bond(mol, link_point[1], right_c);
                add_bond(mol, right_c, mol->addAtom(&atom_o), 2);
                add_bond(mol, right_c, link_point[3]);
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[2], middle_c);
                add_bond(mol, link_point[3], middle_c);
                add_bond(mol, middle_c, mol->addAtom(&atom_o), 2);
                remove_atoms(mol, link_point.begin() + 4, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 6 - h_num, func);
}

ModificationType make_urea(bool thione_form = true, int h_num = 0)
{
    char str[50];
    sprintf(str, "urea (thione_form=%s, h_num=%d)", thione_form ? "true" : "false", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_o("O");
    Atom atom_n("N");
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, thione_form = thione_form, atom_o = atom_o, atom_n = atom_n, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);
                // left_n = link_point[0];
                mol->replaceAtom(link_point[1], &atom_n);
                // right_n = link_point[1];
                if (!thione_form)
                    mol->replaceAtom(link_point[2], &atom_o);
                    // s = link_point[2];
                int c = mol->addAtom(&atom_c);
                if (thione_form)
                    add_bond(mol, c, mol->addAtom(&atom_o), 2);
                else
                    add_bond(mol, c, link_point[2]);
                add_bond(mol, link_point[0], get_target(mol, link_point[thione_form ? 2 : 3]));
                if (thione_form)
                    add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                add_bond(mol, c, link_point[0]);
                add_bond(mol, c, link_point[1], thione_form ? 1 : 2);
                remove_atoms(mol, link_point.begin() + (thione_form ? 2 : 3), link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

ModificationType make_vanillyl()
{
    Atom atom_c("C");
    Atom atom_h("H");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_h = atom_h, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            // right_c = h_list[0];
            int begin_c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], begin_c);
            int current_c, last_c = begin_c;
            for (int i = 0; i < 5; ++i)
            {
                current_c = mol->addAtom(&atom_c);
                mol->addBond(last_c, current_c, Bond::BondType::AROMATIC);
                if (i == 1)
                {
                    int o = mol->addAtom(&atom_o);
                    add_bond(mol, current_c, o);
                    add_bond(mol, o, mol->addAtom(&atom_c));
                }
                if (i == 2)
                {
                    int o = mol->addAtom(&atom_o);
                    add_bond(mol, current_c, o);
                    add_bond(mol, o, mol->addAtom(&atom_h));
                }
                last_c = current_c;
            }
            mol->addBond(last_c, begin_c, Bond::BondType::AROMATIC);
        }
    };
    return std::make_tuple(std::string("vanillyl ()"), 1, func);
}

ModificationType make_vinyl()
{
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            add_bond(mol, h_list[0], mol->addAtom(&atom_c), 2);
        }
    };
    return std::make_tuple(std::string("vinyl ()"), 1, func);
}

ModificationType make_vinylene(int h_num = 0)
{
    char str[50];
    sprintf(str, "vinylene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_vinylidene(int h_num = 0)
{
    char str[50];
    sprintf(str, "vinylidene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 2);
                add_bond(mol, link_point[0], get_target(mol, link_point[2]));
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}

/* Xanthate */
ModificationType make_xanthate(const char* atom_name)
{
    char str[50];
    sprintf(str, "xanthate (atom_name=%s)", atom_name);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_s("S");
    Atom atom_s_negtive("S");
    // atom_s_negtive.setFormalCharge(-1);
    Atom atom(atom_name);
    // atom.setFormalCharge(1);
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c, atom_o = atom_o, atom_s = atom_s, atom_s_negtive = atom_s_negtive, atom = atom] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_o);
            int c = mol->addAtom(&atom_c);
            add_bond(mol, h_list[0], c);
            add_bond(mol, c, mol->addAtom(&atom_s), 2);
            int right_s = mol->addAtom(&atom_s_negtive);
            add_bond(mol, c, right_s);
            add_bond(mol, right_s, mol->addAtom(&atom));
        }
    };
    return std::make_tuple(std::string(str), 1, func);
}

ModificationType make_xanthate_ester(int h_num = 0)
{
    char str[50];
    sprintf(str, "xanthate ester (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                // o = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // s = link_point[1];
                int c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], c);
                add_bond(mol, c, mol->addAtom(&atom_s), 2);
                add_bond(mol, c, link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}
/* Xanthate */

/* Ynolate */
ModificationType make_ynol(int h_num = 0)
{
    char str[50];
    sprintf(str, "ynol (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // bottom_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_o);
                // o = link_point[1];
                int up_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], up_c, 3);
                add_bond(mol, up_c, link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_ynolate(int h_num = 0)
{
    char str[50];
    sprintf(str, "ynolate (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 4 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    Atom atom_si("Si");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o, atom_si = atom_si] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // bottom_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_si);
                // si = link_point[1];
                add_bond(mol, link_point[1], get_target(mol, link_point[2]));
                add_bond(mol, link_point[1], get_target(mol, link_point[3]));
                int up_c = mol->addAtom(&atom_c), o = mol->addAtom(&atom_o);
                add_bond(mol, link_point[0], up_c, 3);
                add_bond(mol, up_c, o);
                add_bond(mol, o, link_point[1]);
                remove_atoms(mol, link_point.begin() + 2, link_point.end());
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 4 - h_num, func);
}
/* Ynolate */

ModificationType make_ynone(int h_num = 0)
{
    char str[50];
    sprintf(str, "ynone (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                int middle_c = mol->addAtom(&atom_c);
                add_bond(mol, link_point[0], middle_c, 3);
                add_bond(mol, middle_c, link_point[1]);
                add_bond(mol, link_point[1], mol->addAtom(&atom_o), 2);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_benzyl(int h_num = 0)
{
    char str[50];
    sprintf(str, "benzyl (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_c);
                // right_c = link_point[1];
                add_bond(mol, link_point[0], link_point[1], 1);
                int next_c, last_c = link_point[1];
                for (int i = 0; i < 5; ++i)
                {
                    next_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, next_c, Bond::BondType::AROMATIC);
                    last_c = next_c;
                }
                mol->addBond(last_c, link_point[1], Bond::BondType::AROMATIC);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_benzene(int h_num = 0)
{
    char str[50];
    sprintf(str, "benzene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // left_c = link_point[0];
                int next_c, last_c = link_point[0];
                for (int i = 0; i < 5; ++i)
                {
                    next_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, next_c, Bond::BondType::AROMATIC);
                    last_c = next_c;
                }
                mol->addBond(last_c, link_point[0], Bond::BondType::AROMATIC);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_sulfide_ether(int h_num = 0)
{
    char str[50];
    sprintf(str, "make_sulfide_ether (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);
    
    Atom atom_c("C");
    Atom atom_s("S");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_s = atom_s] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                // c = link_point[0];
                mol->replaceAtom(link_point[1], &atom_s);
                // s = link_point[1];
                add_bond(mol, link_point[0], link_point[1]);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_methylmorpholine(int h_num = 0)
{
    char str[50];
    sprintf(str, "methylmorpholine (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);

    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);

                int up_c = mol->addAtom(&atom_c), up_c2 = mol->addAtom(&atom_c);
                mol->addBond(link_point[0], up_c, Bond::BondType::SINGLE);
                mol->addBond(up_c, up_c2, Bond::BondType::SINGLE);

                int down_c = mol->addAtom(&atom_c), down_c2 = mol->addAtom(&atom_c);
                mol->addBond(link_point[0], down_c, Bond::BondType::SINGLE);
                mol->addBond(down_c, down_c2, Bond::BondType::SINGLE);

                int o = mol->addAtom(&atom_o);
                mol->addBond(up_c2, o, Bond::BondType::SINGLE);
                mol->addBond(down_c2, o, Bond::BondType::SINGLE);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_dioxolane(int h_num = 0)
{
    char str[50];
    sprintf(str, "dioxolane (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);

    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_o);
                mol->replaceAtom(link_point[1], &atom_o);

                int up_c = mol->addAtom(&atom_c), down_c = mol->addAtom(&atom_c);
                mol->addBond(link_point[0], up_c, Bond::BondType::SINGLE);
                mol->addBond(link_point[1], down_c, Bond::BondType::SINGLE);
                mol->addBond(up_c, down_c, Bond::BondType::SINGLE);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_test4(int pos = 0,int h_num = 0)
{
    char str[50];
    sprintf(str, "test4 (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);

    Atom atom_c("C");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, pos = pos, atom_c = atom_c, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);

                if (pos == 0) {
                    mol->replaceAtom(link_point[0], &atom_o);
                    mol->replaceAtom(link_point[1], &atom_c);
                }
                else {
                    mol->replaceAtom(link_point[1], &atom_o);
                    mol->replaceAtom(link_point[0], &atom_c);
                }

                int c = mol->addAtom(&atom_c);
                mol->addBond(link_point[0], c, Bond::BondType::SINGLE);
                mol->addBond(link_point[1], c, Bond::BondType::SINGLE);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}


ModificationType make_Imide(int h_num = 0)
{
    char str[50];
    sprintf(str, "Imide (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 2 - h_num, 1);

    Atom atom_c("C");
    Atom atom_n("N");
    Atom atom_o("O");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c, atom_n = atom_n, atom_o = atom_o] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_n);

                int up_c = mol->addAtom(&atom_c), up_c2 = mol->addAtom(&atom_c);
                int up_o = mol->addAtom(&atom_o);
                mol->addBond(link_point[0], up_c, Bond::BondType::SINGLE);
                mol->addBond(up_c, up_o, Bond::BondType::DOUBLE);
                mol->addBond(up_c, up_c2, Bond::BondType::SINGLE);

                int down_c = mol->addAtom(&atom_c), down_c2 = mol->addAtom(&atom_c);
                int down_o = mol->addAtom(&atom_o);
                mol->addBond(link_point[0], down_c, Bond::BondType::SINGLE);
                mol->addBond(down_c, down_o, Bond::BondType::DOUBLE);
                mol->addBond(down_c, down_c2, Bond::BondType::SINGLE);

                int last_c = up_c2;
                int next_c;
                // C6H5
                for (int i = 0; i < 4; ++i)
                {
                    next_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, next_c, Bond::BondType::AROMATIC);
                    last_c = next_c;
                }
                mol->addBond(last_c, down_c2, Bond::BondType::AROMATIC);
                mol->addBond(down_c2, up_c2, Bond::BondType::AROMATIC);
            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_naphthalene(int h_num = 0)
{
    char str[50];
    sprintf(str, "Naphthalene (h_num=%d)", h_num);

    std::vector<int> permutation(h_num, -1);
    permutation.insert(permutation.end(), 3 - h_num, 1);

    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func;
    do
    {
        func.emplace_back(
            [permutation = permutation, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
            {
                std::vector<uint8_t> link_point = calculate_link_point(mol, permutation, h_list);
                mol->replaceAtom(link_point[0], &atom_c);
                mol->replaceAtom(link_point[1], &atom_c);

                int last_c = link_point[0];
                int next_c;
                // C6H5
                for (int i = 0; i < 2; ++i)
                {
                    next_c = mol->addAtom(&atom_c);
                    mol->addBond(last_c, next_c, Bond::BondType::AROMATIC);
                    last_c = next_c;
                }
                mol->addBond(last_c, link_point[1], Bond::BondType::AROMATIC);
                // mol->addBond(link_point[1], link_point[0], Bond::BondType::AROMATIC);

            }
        );
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return std::make_tuple(std::string(str), 2 - h_num, func);
}

ModificationType make_ben1(int h_num = 0)
{
    char str[50];
    sprintf(str, "ben1 (h_num=%d)", h_num);
    
    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            auto neighbors = mol->getAtomNeighbors(mol->getAtomWithIdx(h_list[0]));
            auto first_atom = (*mol)[*boost::begin(neighbors)];
            neighbors = mol->getAtomNeighbors(mol->getAtomWithIdx(first_atom->getIdx()));

            if (first_atom->getSymbol() != "C")
                return;

            int num_h = 0;
            std::vector<int> h_idx;
            for (const auto &nbri : boost::make_iterator_range(neighbors))
                if ((*mol)[nbri]->getSymbol() == "H")
                {
                    h_idx.push_back((*mol)[nbri]->getIdx());
                    // std::cout << (*mol)[nbri]->getIdx() << " ";
                    num_h++;
                }
                
            // std::cout << "\n";

            if (num_h == 3)
            {
                int first_c = first_atom->getIdx();
                mol->replaceAtom(h_idx[0], &atom_c);
                mol->replaceAtom(h_idx[1], &atom_c);
                mol->removeAtom(h_idx[2]);

                auto* bond = mol->getBondBetweenAtoms(first_c, h_idx[0]);
                Bond* new_bond = new Bond(Bond::BondType::DOUBLE);
                mol->replaceBond(bond->getIdx(), new_bond);

                int prev_c = h_idx[0];
                int next_c;
                // C6H5
                for (int i = 0; i < 3; ++i)
                {
                    next_c = mol->addAtom(&atom_c);
                    if (i%2 == 1)
                        mol->addBond(prev_c, next_c, Bond::BondType::DOUBLE);
                    else
                        mol->addBond(prev_c, next_c, Bond::BondType::SINGLE);
                    prev_c = next_c;
                }

                mol->addBond(prev_c, h_idx[1], Bond::BondType::DOUBLE);

                bond = mol->getBondBetweenAtoms(h_idx[1], first_c);
                delete new_bond;
                new_bond = new Bond(Bond::BondType::SINGLE);
                mol->replaceBond(bond->getIdx(), new_bond);
                delete new_bond;
            }
            // remove_atoms(mol, h_list.begin(), h_list.end());
        }
    };
    return std::make_tuple(std::string(str), 2, func);
}

ModificationType make_dialkyl(int num, int h_num = 0)
{
    char str[50];
    sprintf(str, "dialkyl (num=%d, h_num=%d)", num, h_num);

    Atom atom_c("C");
    std::vector<std::function<void(RWMOL_SPTR, std::vector<uint8_t>&)>> func = {
        [num = num, h_num = h_num, atom_c = atom_c] (RWMOL_SPTR mol, std::vector<uint8_t>& h_list) mutable
        {
            mol->replaceAtom(h_list[0], &atom_c);
            int new_c;
            int last_c = h_list[0];
            for (int i = 0; i < num - 2; ++i)
            {
                new_c = mol->addAtom(&atom_c);
                add_bond(mol, last_c, new_c);
                last_c = new_c;
            }
            if (h_num == 1 && num!=0)
            {
                new_c = mol->addAtom(&atom_c);
                add_bond(mol, last_c, new_c);
            }
            else if (h_num==0 && num!=0)
            {
                mol->replaceAtom(h_list[1], &atom_c);
                add_bond(mol, last_c, h_list[1]);
            }
        }
    };
    return std::make_tuple(std::string(str), 2 - h_num, func);
}
