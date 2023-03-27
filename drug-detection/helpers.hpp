#include <tuple>
#include <set>
#include <map>
#include <algorithm>
#include <shared_mutex>
#include <stdexcept>
#include "graph.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
using namespace RDKit; 

struct RingSearch
{
  public:
    bool find(int key)
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return data.find(key) != data.end();
    }

    RWMOL_SPTR operator[](int key)
    {
        std::shared_lock<std::shared_mutex> lock(mutex);
        return data[key];
    }

    void emplace(int key, RWMOL_SPTR val)
    {
        std::unique_lock<std::shared_mutex> lock(mutex);
        data.emplace(key, val);
    }

  private:
    mutable std::shared_mutex mutex;
    std::map<int, RWMOL_SPTR> data;
};
RingSearch ring_search;

double calculate_weight(RWMOL_SPTR mol)
{
    double mass = 0;
    for (auto& atom : mol->atoms())
        mass += atom->getMass();
    return mass;
}

bool is_planar(RWMOL_SPTR mol)
{
    graphP thegraph = gp_New();
    if (gp_InitGraph(thegraph, mol->getNumAtoms()) != 0)
    {
        std::cerr << "[ERROR] gp_InitGraph() failed\n";
        // throw std::invalid_argument("Molecule could not vlidate by is_planar()");
        return true;
    }
    int node1, node2;
    for (auto& bond : mol->bonds())
    {
        node1 = bond->getBeginAtomIdx();
        node2 = bond->getEndAtomIdx();
        gp_AddEdge(thegraph, node1, 0, node2, 0);
    }
    int result = gp_Embed(thegraph, 1);
    gp_Free(&thegraph);
    if (result == -3)
        return false;
    return true;
}

std::vector<Atom*> match_atoms(RWMOL_SPTR parent, std::vector<int> match)
{
    std::vector<Atom*> atoms;
    atoms.reserve(match.size());
    for (auto& idx : match)
        atoms.emplace_back(parent->getAtomWithIdx(idx));
    return atoms;
}

std::vector<Bond*> match_bonds(RWMOL_SPTR parent, std::vector<int> match, RWMOL_SPTR pattern)
{
    std::vector<Bond*> bonds;
    bonds.reserve(pattern->getNumBonds());
    for (auto& bond : pattern->bonds())
        bonds.emplace_back(parent->getBondBetweenAtoms(match[bond->getBeginAtomIdx()], match[bond->getEndAtomIdx()]));
    return bonds;
}

std::vector<std::vector<int>> match_mol_id(RWMOL_SPTR mol, RWMOL_SPTR patt)
{
    std::vector<MatchVectType> matches = SubstructMatch(*mol, *patt);
    std::vector<std::vector<int>> matchids(matches.size());
    for (unsigned int i = 0; i < matches.size(); ++i)
    {
        matchids[i].reserve(matches[i].size());
        for (auto& item : matches[i])
            matchids[i].emplace_back(item.second);
    }
    return matchids;
}

std::tuple<std::vector<int>, int, int> SSSR(RWMOL_SPTR mol, bool force = false)
{
    std::vector<int> ringcounts;
    if (!force && mol->getPropIfPresent("ringcounts", ringcounts))
    {
        std::pair<int, int> buf;
        mol->getProp("sharedcounts", buf);
        return std::make_tuple(ringcounts, buf.first, buf.second);
    }
    
    std::set<int> assigned_atoms;
    std::set<int> assigned_bonds;
    std::set<int> shared_atoms;
    std::set<int> shared_bonds;
    RingInfo* ring_info = mol->getRingInfo();
    unsigned int nringatom = 0;
    for (auto& atom : mol->atoms())
        if (ring_info->numAtomRings(atom->getIdx()) != 0)
            ++nringatom;
    unsigned int nringbond = 0;
    for (auto& bond : mol->bonds())
        if (ring_info->numBondRings(bond->getIdx()) != 0)
            ++nringbond;
    
    // Loop over all possible ring sizes
    std::vector<int> nrings(nringatom > 8 ? nringatom : 8, 0);
    std::string buf;
    buf.reserve(100);
    std::set<int> bondids;
    for (unsigned int i = 3; i <= nringatom; ++i)
    {
        if (assigned_atoms.size() == nringatom && assigned_bonds.size() == nringbond)
            break;
        if (!ring_search.find(i))
        {
            buf = "*~1";
            for (unsigned int j = 0; j < i - 1; ++j)
                buf += "~*";
            ring_search.emplace(i, RWMOL_SPTR(SmartsToMol(buf + '1')));
        }
        // Find all instances of ring size i
        std::vector<std::vector<int>> matches = match_mol_id(mol, ring_search[i]);
        bool atom_is_subset = true;
        bool bond_is_subset = true;
        for (auto& match : matches)
        {
            bondids.clear();
            for (auto& bond : match_bonds(mol, match, ring_search[i]))
                bondids.emplace(bond->getIdx());
            // Count this ring only if some of its atoms or bonds
            // have not already been assigned to a smaller ring
            for (auto& id : match)
                if (assigned_atoms.find(id) == assigned_atoms.end())
                {
                    atom_is_subset = false;
                    break;
                }
            if (atom_is_subset)
            {
                for (auto& id : bondids)
                    if (assigned_bonds.find(id) == assigned_bonds.end())
                    {
                        bond_is_subset = false;
                        break;
                    }
            }
            if (!(atom_is_subset && bond_is_subset))
            {
                ++nrings[i - 3];
                std::vector<int> intersection;
                std::sort(match.begin(), match.end());
                std::set_intersection(assigned_atoms.begin(), assigned_atoms.end()
                                     , match.begin(), match.end()
                                     , std::back_inserter(intersection));
                shared_atoms.insert(intersection.begin(), intersection.end());
                intersection.clear();
                std::set_intersection(assigned_bonds.begin(), assigned_bonds.end()
                                     , bondids.begin(), bondids.end()
                                     , std::back_inserter(intersection));
                shared_bonds.insert(intersection.begin(), intersection.end());
                assigned_atoms.insert(match.begin(), match.end());
                assigned_bonds.insert(bondids.begin(), bondids.end());
            }
        }
    }
    mol->setProp("ringcounts", nrings);
    mol->setProp("sharedcounts", std::make_pair((int)shared_atoms.size(), (int)shared_bonds.size()));
    return std::make_tuple(nrings, shared_bonds.size(), shared_bonds.size());
}

std::tuple<std::vector<std::set<int>>, std::vector<std::set<int>>> SSSR_getrings(RWMOL_SPTR mol, bool force = false)
{
    int num_ring;
    std::vector<std::set<int>> ringatoms;
    std::vector<std::set<int>> ringbonds;
    if (!force && mol->getPropIfPresent("numSSSRrings", num_ring))
    {
        mol->getProp("atomSSSR", ringatoms);
        mol->getProp("bondSSSR", ringbonds);
        return std::make_tuple(ringatoms, ringbonds);
    }
    
    std::set<int> assigned_atoms;
    std::set<int> assigned_bonds;
    std::set<int> shared_atoms;
    std::set<int> shared_bonds;
    RingInfo* ring_info = mol->getRingInfo();
    unsigned int nringatom = 0;
    for (auto& atom : mol->atoms())
        if (ring_info->numAtomRings(atom->getIdx()) != 0)
            ++nringatom;
    unsigned int nringbond = 0;
    for (auto& bond : mol->bonds())
        if (ring_info->numBondRings(bond->getIdx()) != 0)
            ++nringbond;
    
    // Loop over all possible ring sizes
    std::vector<int> nrings(nringatom > 8 ? nringatom : 8, 0);
    std::string buf;
    buf.reserve(100);
    std::set<int> bondids;
    for (unsigned int i = 3; i <= nringatom; ++i)
    {
        if (assigned_atoms.size() == nringatom && assigned_bonds.size() == nringbond)
            break;
        if (!ring_search.find(i))
        {
            buf = "*~1";
            for (unsigned int j = 0; j < i - 1; ++j)
                buf += "~*";
            ring_search.emplace(i, RWMOL_SPTR(SmartsToMol(buf + '1')));
        }
        // Find all instances of ring size i
        std::vector<std::vector<int>> matches = match_mol_id(mol, ring_search[i]);
        bool atom_is_subset = true;
        bool bond_is_subset = true;
        for (auto& match : matches)
        {
            if (assigned_atoms.size() == nringatom && assigned_bonds.size() == nringbond)
                break;
            bondids.clear();
            for (auto& bond : match_bonds(mol, match, ring_search[i]))
                bondids.emplace(bond->getIdx());
            // Count this ring only if some of its atoms or bonds
            // have not already been assigned to a smaller ring
            for (auto& id : match)
                if (assigned_atoms.find(id) == assigned_atoms.end())
                {
                    atom_is_subset = false;
                    break;
                }
            if (atom_is_subset)
            {
                for (auto& id : bondids)
                    if (assigned_bonds.find(id) == assigned_bonds.end())
                    {
                        bond_is_subset = false;
                        break;
                    }
            }
            if (!(atom_is_subset && bond_is_subset))
            {
                ringatoms.emplace_back(match.begin(), match.end());
                ringbonds.emplace_back(bondids);
                ++nrings[i - 3];
                std::vector<int> intersection;
                std::set_intersection(assigned_atoms.begin(), assigned_atoms.end()
                                     , ringatoms.back().begin(), ringatoms.back().end()
                                     , std::back_inserter(intersection));
                shared_atoms.insert(intersection.begin(), intersection.end());
                intersection.clear();
                std::set_intersection(assigned_bonds.begin(), assigned_bonds.end()
                                     , bondids.begin(), bondids.end()
                                     , std::back_inserter(intersection));
                shared_bonds.insert(intersection.begin(), intersection.end());
                assigned_atoms.insert(match.begin(), match.end());
                assigned_bonds.insert(bondids.begin(), bondids.end());
            }
        }
    }
    mol->setProp("ringcounts", nrings);
    mol->setProp("sharedcounts", std::make_pair((int)shared_atoms.size(), (int)shared_bonds.size()));
    mol->setProp("numSSSRrings", (int)ringatoms.size());
    mol->setProp("atomSSSR", ringatoms);
    mol->setProp("bondSSSR", ringbonds);
    return std::make_tuple(ringatoms, ringbonds);
}

std::tuple<std::set<int>, std::set<int>> find_bridge_atoms(RWMOL_SPTR mol, bool force = false)
{
    std::vector<std::set<int>> rings = std::get<0>(SSSR_getrings(mol, force));
    std::set<int> bridge_atoms;
    std::set<int> multiple_bridge;
    std::set<int> isct;
    for (unsigned int i = 0; i < rings.size(); ++i)
    {
        if (rings[i].size() > 8)
            continue;
        for (unsigned int j = i + 1; j < rings.size(); ++j)
        {
            if (rings[j].size() > 8)
                continue;
            std::set_intersection(rings[i].begin(), rings[i].end()
                                , rings[j].begin(), rings[j].end()
                                , std::inserter(isct, isct.begin()));
            if (isct.size() > 2)
            {
                for (auto& id : isct)
                {
                    if (bridge_atoms.find(id) != bridge_atoms.end())
                        multiple_bridge.insert(id);
                    else
                        bridge_atoms.insert(id);
                }
            }
        }
    }
    return std::make_tuple(bridge_atoms, multiple_bridge);
}

std::set<Bond*> find_bridge_bonds(RWMOL_SPTR mol, bool force = false)
{
    std::vector<std::set<int>> rings = std::get<1>(SSSR_getrings(mol, force));
    std::set<Bond*> bridge_bonds;
    std::set<int> isct;
    Bond* bd;
    for (unsigned int i = 0; i < rings.size(); ++i)
    {
        if (rings[i].size() > 8)
            continue;
        for (unsigned int j = i + 1; j < rings.size(); ++j)
        {
            if (rings[j].size() > 8)
                continue;
            std::set_intersection(rings[i].begin(), rings[i].end()
                                , rings[j].begin(), rings[j].end()
                                , std::inserter(isct, isct.begin()));
            if (isct.size() > 1)
            {
                for (auto& id : isct)
                {
                    bd = mol->getBondWithIdx(id);
                    if (bridge_bonds.find(bd) == bridge_bonds.end())
                        bridge_bonds.insert(bd);
                }
            }
        }
    }
    return bridge_bonds;
}

std::set<Atom*> find_bridge_heads(RWMOL_SPTR mol, bool force = false)
{
    std::vector<std::set<int>> rings = std::get<0>(SSSR_getrings(mol, force));
    std::vector<std::set<int>> bridges;
    std::set<int> isct;
    for (unsigned int i = 0; i < rings.size(); ++i)
    {
        if (rings[i].size() > 8)
            continue;
        for (unsigned int j = i + 1; j < rings.size(); ++j)
        {
            if (rings[j].size() > 8)
                continue;
            isct.clear();
            std::set_intersection(rings[i].begin(), rings[i].end()
                                 , rings[j].begin(), rings[j].end()
                                 , std::inserter(isct, isct.begin()));
            if (isct.size() > 2)
                bridges.emplace_back(isct);
        }
    }
    std::set<Atom*> bridge_heads;
    std::set<int> nbid;
    int idd;
    for (auto& bridge : bridges)
    {
        for (auto& atid : bridge)
        {
            nbid.clear();
            for (const auto &nbri : boost::make_iterator_range(mol->getAtomNeighbors(mol->getAtomWithIdx(atid))))
            {
                idd = (int)(*mol)[nbri]->getIdx();
                if (bridge.find(idd) != bridge.end())
                    nbid.emplace(idd);
            }
            if (nbid.size() == 1)
                bridge_heads.emplace(mol->getAtomWithIdx(atid));
        }
    }
    return bridge_heads;
}

void print_molecule(RWMOL_SPTR mol)
{
    for (auto& atom : mol->atoms())
    {
        std::cout << atom->getIdx() << " " << atom->getSymbol() << ": ";
        for (const auto &nbri : boost::make_iterator_range(mol->getAtomNeighbors(atom)))
            std::cout << (*mol)[nbri]->getIdx() << " " << (*mol)[nbri]->getSymbol() << ", ";
        std::cout << "\n";
    }
}