#include "GDBFilter.hpp"
#include <GraphMol/Depictor/RDDepictor.h>

class SMUFilter : public GDBFilter
{
  public:
    SMUFilter() : GDBFilter()
    {
        Filter::opt_sulfone.emplace_back(SmartsToMol("[#6]S[#6]"));
        init_filters();
    }

    std::string geom_filter(RWMOL_SPTR mol)
    {
        std::string failname = GDBFilter::geom_filter(mol);
        if (!failname.empty())
            return failname;

        std::vector<RDGeom::Point3D> coords;
        coords.reserve(mol->getNumAtoms());
        Conformer& conf = mol->getConformer(RDDepict::compute2DCoords(*mol));
        for (unsigned int i = 0;i < coords.size(); ++i)
            coords[i] = conf.getAtomPos(i);
        // Check sp2 carbons for planarity
        RDGeom::Point3D normalvec;
        for (auto& match : match_mol_id(mol, RWMOL_SPTR(SmartsToMol("*C(*)=*"))))
        {
            RDGeom::Point3D& center = coords[match[1]];
            normalvec = center.directionVector(coords[match[0]]).crossProduct(
                center.directionVector(coords[match[2]])
            );
            normalvec.normalize();
            if (std::abs(normalvec.dotProduct(center.directionVector(coords[match[3]]))) > 0.15)
                return "Non-planar sp2 system";
        }
        // check triple bonds for linearity
        for (auto& match : match_mol_id(mol, RWMOL_SPTR(SmartsToMol("**#*"))))
        {
            RDGeom::Point3D& center = coords[match[1]];
            if (center.directionVector(coords[match[0]]).dotProduct(
                    center.directionVector(coords[match[2]])
                ) > -0.999657325)
                return "Non-linear sp system";
        }
        // If here, geometry passes
        return failname;
    }
  protected:
    void init_filters()
    {
        Filter::filters[BREDT_RULE] = Filter::NEW_FILTER_SPTR(new NewFilter("Bredt's rule",
            [] (RWMOL_SPTR mol) -> std::string
            {
                std::set<Atom*> bridge_heads = find_bridge_heads(mol);
                for (auto& atom : bridge_heads)
                    for (const auto &nbri : boost::make_iterator_range(mol->getAtomBonds(atom)))
                        if ((*mol)[nbri]->getBondType() == Bond::BondType::AROMATIC
                             || (*mol)[nbri]->getBondTypeAsDouble() > 1)
                            return "Bredt violation";
                return std::string();
            }
        )); 
        Filter::filters.emplace(TERMINAL_SULFUR_1, new NewPatternFilter(
            "TerminalSulfur1", RWMOL_SPTR(SmartsToMol("[S&H1]"))
        )); 
        Filter::filters.emplace(TERMINAL_SULFUR_2, new NewPatternFilter(
            "TerminalSulfur2", RWMOL_SPTR(SmartsToMol("[!O]=S")), 
            { Filter::thiourea }
        )); 
        static_cast<NewPatternFilter*>(Filter::filters[HET_HET_NN].get())->set_exceptions({
            RWMOL_SPTR(SmartsToMol("[#6]N([#6])N=C"))
          , RWMOL_SPTR(SmartsToMol("C=N[N&H1][#6]"))
          , RWMOL_SPTR(SmartsToMol("C=N[N&H2]"))
          , RWMOL_SPTR(SmartsToMol("C=NN=C"))
          , RWMOL_SPTR(SmartsToMol("cN=Nc"))
          , RWMOL_SPTR(SmartsToMol("NNC=[N,O]"))
          , RWMOL_SPTR(SmartsToMol("[N&H2]Nc"))
        });
        Filter::filters.erase(HET_HET_OS_1);
        Filter::filters.erase(HET_HET_OS_2);
        Filter::filters.erase(HET_HET_OS_3);
        Filter::filters.erase(ORTHOESTER);
        Filter::filters.emplace(HET_HET_OS, new NewPatternFilter(
            "HetHet_OS", RWMOL_SPTR(SmartsToMol("[#8][#16]"))
        )); 
        Filter::filters[HET_HET_NS] = Filter::NEW_FILTER_SPTR(new NewPatternFilter(
            "HetHet_NS", RWMOL_SPTR(SmartsToMol("NS")), 
            { RWMOL_SPTR(SmartsToMol("NS(=O)(=O)c"))
            , RWMOL_SPTR(SmartsToMol("cNS(=O)(=O)"))}
        )); 
        Filter::filters.emplace(TOO_MANY_STEREOCENTER, new NewFilter(">3 stereocenters",
            [] (RWMOL_SPTR mol) -> std::string
            {
                MolOps::assignStereochemistry(*mol, false, true, false);
                std::vector<char> chiral_list;
                std::string code;
                unsigned int num = 0;
                for (auto& atom : mol->atoms() )
                    if (atom->getPropIfPresent("_CIPCode", code) && code == "C")
                        ++num;
                if (num > 3)
                    return ">3 stereocenters";
                return std::string();
            }
        )); 
        Filter::filters.emplace(TRIPLE_BOND_TOO_CLOSE, new NewFilter("Triple bonds too close",
            [] (RWMOL_SPTR mol) -> std::string
            {
                std::vector<Bond*> tbonds;
                tbonds.reserve(mol->getNumBonds());
                for (auto& bond : mol->bonds())
                    if (bond->getBondType() == Bond::BondType::TRIPLE)
                        tbonds.emplace_back(bond);
                Bond* tb;
                std::vector<unsigned int(Bond::*)() const> funcs
                    = { &Bond::getBeginAtomIdx, &Bond::getEndAtomIdx };
                unsigned int dist;
                while (tbonds.size() > 1)
                {
                    tb = tbonds.back();
                    tbonds.pop_back();
                    for (auto& bond : tbonds)
                    {
                        for (auto& tb_func : funcs)
                            for (auto& bond_func : funcs)
                            {
                                dist = MolOps::getShortestPath(*mol, (tb->*tb_func)(), (bond->*bond_func)()).size();
                                if (dist < 8)
                                    return "Triple bonds too close";
                                if (dist >= 9)
                                    break;
                            }
                    }
                }
                return std::string();
            }
        )); 
        Filter::filters.emplace(MULTIPLY_BRIDGED_RING, new NewFilter("Multiply-bridged ring",
            [] (RWMOL_SPTR mol) -> std::string
            {
                std::set<int> adatoms;
                for (auto& ad : match_mol_id(mol, adamantane))
                    adatoms.insert(ad.begin(), ad.end());
                std::vector<std::set<int>> rings = std::get<0>(SSSR_getrings(mol));
                unsigned int nring = rings.size();
                std::vector<bool> bridged(nring, false);
                std::set<int> isct;
                bool is_subset;
                for (unsigned int i = 0; i < nring; ++i)
                {
                    if (rings[i].size() > 8)
                        continue;
                    for (unsigned int j = i + 1; j < nring; ++j)
                    {
                        if (rings[j].size() > 8)
                            continue;
                        isct.clear();
                        std::set_intersection(rings[i].begin(), rings[i].end()
                                            , rings[j].begin(), rings[j].end()
                                            , std::inserter(isct, isct.begin()));
                        if (isct.size() > 2)
                        {
                            is_subset = true;
                            for (auto& item : isct)
                                if (adatoms.find(item) == adatoms.end())
                                {
                                    is_subset = false;
                                    break;
                                }
                            if (!is_subset)
                            {
                                if (bridged[i] || bridged[j])
                                    return "Multiply-bridged ring";
                                else
                                {
                                    bridged[i] = true;
                                    bridged[j] = true;
                                }
                            }
                        }
                    }
                }
                return std::string();
            }
        )); 
        Filter::filters.emplace(MESSY_RING, new NewFilter("Messy rings",
            [] (RWMOL_SPTR mol) -> std::string
            {
                std::set<int> multiple_bridge = std::get<1>(find_bridge_atoms(mol));
                unsigned int n_share_bridge = multiple_bridge.size();
                if (n_share_bridge > 1)
                    n_share_bridge -= 2 * SubstructMatch(*mol, *adamantane).size();
                if (n_share_bridge > 1)
                    return "Messy rings";
                return std::string();
            }
        ));
        Filter::filters.emplace(UNSATURATION_BRIDGE, new NewFilter("Unsaturation in bridge",
            [] (RWMOL_SPTR mol) -> std::string
            {
                for (auto& bond : find_bridge_bonds(mol))
                    if (bond->getBondTypeAsDouble() > 1 || bond->getBondType() == Bond::BondType::AROMATIC)
                        return "Unsaturation in bridge";
                return std::string();
            }
        ));
    }

    static RWMOL_SPTR adamantane;
};

RWMOL_SPTR SMUFilter::adamantane = RWMOL_SPTR(SmartsToMol("C3C4CC5CC(C4)CC3C5"));