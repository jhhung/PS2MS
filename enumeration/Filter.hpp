#include <string>
#include <functional>
#include <vector>
#include <set>
#include <numeric>
#include <map>
#include "helpers.hpp"

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/MolOps.h>
using namespace RDKit; 

enum FILTER_NAME
{
    BREDT_RULE
  , BETA_HATERO_CARBONYL
  , TOO_MANY_AMINE
  , MESSY_RING
  , HET_HET_HET
  , DOUBLE_BOND_IN_3_RING
  , HET_HET_NO
  , DOUBLE_BOND_IN_4_RING
  , HEMIAMINAL_2
  , HETERO_SR_5
  , HETERO_SR_2
  , HEMIAMINAL_3
  , HETERO_SR_1
  , HETERO_SR_3
  , HET_HET_NN
  , HEMIAMINAL_1
  , HEMIACETAL
  , INTRAMOL
  , POLYCYCLIC_1
  , UNSATURATION_BRIDGE
  , NON_AROMATIC_HALOGEN
  , HETERO_TRIPLE_BOND
  , FC_2
  , FC_1
  , ENOL_ETHER
  , SSSR_TOO_BIG
  , ORTHO_1
  , MIXED_7
  , TOPO1_34FUSE
  , TRIPLE_BOND_IN_RING
  , HET_HET_OO
  , TOPO1_33FUSE
  , HETERO_SR_4
  , C_DOUBLE_N
  , ENAMINE
  , TOPO1_44FUSE
  , TRIPLE_BOND_TOO_CLOSE
  , ALLENE
  , AMINAL
  , HET_HET_AROMATIC_2
  , ENOL
  , ORTHO_2
  , HET_HET_AROMATIC
  , HETERO_SR_7
  , TOPO1_34SPIRO
  , TOPO1_44SPIRO
  , HETERO_SR_6
  , ANHYDRIDE
  , NON_PLANAR_BOYES
  , TOPO1_33SPIRO
  , TOO_MANY_TERMINAL_DOUBLE_BOND_C
  , TOO_BIG
  , NON_PLANAR_EULER
  , TOO_MANY_RING
  , ATOM_COUNT
  , ACIDTAUT
  , DECARBOXY_1
  , DECARBOXY_2
  , GEMINAL_1
  , GEMINAL_2
  , HET_HET_NS
  , HET_HET_OS_1
  , HET_HET_OS_2
  , HET_HET_OS_3
  , HET_HET_XO
  , HET_HET_XN
  , HET_HET_SS
  , HETERO_SULFUR_WITHOUT_SULFONE
  , MIXED_1
  , MIXED_2
  , MIXED_3
  , MIXED_4
  , MIXED_5
  , MIXED_6
  , ORTHOESTER
  , TOPO1_44BRIDGE
  , NON_AROMATIC_NITRO
  , HETERO_AROMATIC_HALOGEN
  , POLYCYCLIC_2
  , BETA_KETO_CARBOXYL
  , TERMINAL_SULFUR_1
  , TERMINAL_SULFUR_2
  , HET_HET_OS
  , TOO_MANY_STEREOCENTER
  , MULTIPLY_BRIDGED_RING
  , NO_N_OR_O
  , TOO_MANY_HALOGEN
  , TOO_MANY_ALDEHYDE
  , TOO_MANY_ENE
  , TOO_SMALL
  , DRUG_LIKENESS
  , HET_SULFONE_HET
  , PERHALO_KETONE
  , POLYCYCLIC_SULFUR
  , HALOPYRIMIDINE
  , LONGALKYL
  , TERMINAL_DOUBLE_BOND
};

class NewFilter
{
  public:
    NewFilter() {}

    NewFilter(const std::string& myname, const std::function<std::string(RWMOL_SPTR)>& filter_routine)
    : name      (myname)
    , routine   (filter_routine)
    {      
    }

    void set_filter_routine(const std::function<std::string(RWMOL_SPTR)>& filter_routine)
    {
        routine = filter_routine;
    }

    virtual std::string do_filt(RWMOL_SPTR mol) const
    {
        return routine(mol);
    }

    std::string get_name() const
    {
        return name;
    }

    std::string operator()(RWMOL_SPTR mol) const
    {
        return do_filt(mol);
    }

  protected:
    NewFilter(const std::string& myname)
    : name  (myname)
    {}

    std::string name;
    std::function<std::string(RWMOL_SPTR)> routine;
};

class NewPatternFilter : public NewFilter
{
  public:
    NewPatternFilter(const std::string& myname, RWMOL_SPTR patt, RWMOL_SPTR_VECT exc = RWMOL_SPTR_VECT())
    : NewFilter     (myname)
    , pattern       (patt)
    , exceptions    (exc)
    {
    }

    void set_filter_pattern(RWMOL_SPTR pattern)
    {
        this->pattern = pattern;
    }

    void set_exceptions(RWMOL_SPTR_VECT exc)
    {
        exceptions = exc;
    }

    std::vector<std::vector<int>> filter_with_exceptions(RWMOL_SPTR mol) const
    {
        std::vector<std::vector<int>> matchids = match_mol_id(mol, pattern);
        if (exceptions.empty())
            return matchids;
        
        std::vector<int> remove_idx;
        remove_idx.reserve(matchids.size());
        std::vector<std::vector<int>> exmatches;
        bool is_subset;
        // Remove matches that are substructures of exceptions
        for (auto& exception : exceptions)
        {
            exmatches = match_mol_id(mol, exception);
            // for each exception substructure found:
            for (auto& em : exmatches)
            {
                std::set exceptids(em.begin(), em.end());
                remove_idx.clear();
                for (int i = (int)matchids.size() - 1; i >= 0; --i)
                {
                    is_subset = true;
                    for (auto& id : matchids[i])
                    {
                        if (exceptids.find(id) == exceptids.end())
                        {
                            is_subset = false;
                            break;
                        }
                    }
                    if (is_subset)
                        remove_idx.emplace_back(i);
                }
                for (auto& idx : remove_idx)
                    matchids.erase(matchids.begin() + idx);
                if (matchids.size() == 0)
                    break;
            }
            if (matchids.size() == 0)
                break;
        }
        return matchids;
    }

    std::string do_filt(RWMOL_SPTR mol) const 
    {
        if (!exceptions.empty())
        {
            if (!filter_with_exceptions(mol).empty())
                return name;
        }
        else
        {
            MatchVectType result;
            if (SubstructMatch(*mol, *pattern, result))
                return name;
        }
        return std::string();
    }

  protected:
    RWMOL_SPTR pattern;
    RWMOL_SPTR_VECT exceptions;
};

class Filter
{
  public:
    Filter()
    : Filter({}, {})
    {
    }

    virtual std::string geom_filter(RWMOL_SPTR mol) = 0;

    std::string operator()(RWMOL_SPTR mol)
    {
        // First, try geometry filter
        std::string failure;// = geom_filter(mol);
        if (failure.empty())
        {
            // Run through all filters
            // MolOps::removeHs(*mol);
            for (auto& [name, filter] : filters)
            {
                failure = (*filter)(mol);
                if (!failure.empty())
                    break;
            }
        }
        // mol->setProp("filtered", true);
        // mol->setProp("failedfilter", failure);
        return failure;
    }

  protected:
    Filter(RWMOL_SPTR_VECT no, RWMOL_SPTR_VECT opt)
    : max_rings             (8)
    , max_ring_size         (8)
    , ring_size_exceptions  (1)
    , constraint_init       (false)
    , set_neutral_ph        (false)
    , max_weight            (0.0)
    , max_atom              (-1)
    , gen_structs           (false)
    , max_try               (10)
    , no_sulfone            (no)
    , opt_sulfone           (opt)
    {
        init_filters();
    }

    void init_filters()
    {
        filters.emplace(TOO_BIG, new NewFilter("Too big",
            [max_atom = max_atom, max_weight = max_weight] (RWMOL_SPTR mol) -> std::string
            {
                if (max_atom > 0 && (int)mol->getNumHeavyAtoms() > max_atom)
                    return "Too big";
                if (max_weight > 0 && calculate_weight(mol) > max_weight)
                    return "Too big";
                return std::string();
            }
        ));
        filters.emplace(NON_PLANAR_EULER, new NewFilter("Non-planar graph (Euler critereon)",
            [] (RWMOL_SPTR mol) -> std::string
            {
                unsigned int nnodes = mol->getNumHeavyAtoms();
                if (nnodes >= 3 && mol->getNumBonds() > 3 * nnodes - 6)
                    return "Non-planar graph (Euler critereon)";
                return std::string();
            }
        )); 
        filters.emplace(NON_PLANAR_BOYES, new NewFilter("Non-planar graph (Boyes)",
            [] (RWMOL_SPTR mol) -> std::string
            {
                if (!is_planar(mol))
                    return "Non-planar graph (Boyes)";
                return std::string();
            }
        )); 
        filters.emplace(TOO_MANY_RING, new NewFilter("Too many rings",
            [max_rings = max_rings] (RWMOL_SPTR mol) -> std::string
            {
                if (mol->getNumAtoms() > max_rings)
                {
                    std::vector<int> nrings = std::get<0>(SSSR(mol));
                    if (std::accumulate(nrings.begin(), nrings.end(), 0) > (int)max_rings)
                        return "Too many rings";
                }
                return std::string();
            }
        )); 
        filters.emplace(SSSR_TOO_BIG, new NewFilter("SSSR ring bigger than max allowed size",
            [max_ring_size = max_ring_size, ring_size_exceptions = ring_size_exceptions] (RWMOL_SPTR mol) -> std::string
            {
                if (mol->getNumAtoms() > max_ring_size)
                {
                    std::vector<int> nrings = std::get<0>(SSSR(mol));
                    if (std::accumulate(nrings.begin() + (max_ring_size - 2), nrings.end(), 0) > (int)ring_size_exceptions)
                        return "SSSR ring bigger than max allowed size";
                }
                return std::string();
            }
        )); 
    }

    using NEW_FILTER_SPTR = boost::shared_ptr<NewFilter>;
    const std::uint32_t max_rings;
    const std::uint32_t max_ring_size;
    const std::int32_t ring_size_exceptions;
    const bool constraint_init;
    const bool set_neutral_ph;
    const double max_weight;
    const std::int32_t max_atom;
    const bool gen_structs;
    const std::int32_t max_try;
    RWMOL_SPTR_VECT no_sulfone;
    RWMOL_SPTR_VECT opt_sulfone;
    std::map<FILTER_NAME, NEW_FILTER_SPTR> filters;

    static RWMOL_SPTR thiourea;
};

RWMOL_SPTR Filter::thiourea = RWMOL_SPTR(SmartsToMol("NC(=S)N"));
