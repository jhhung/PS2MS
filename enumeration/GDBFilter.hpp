#include "Filter.hpp"
#include <cmath>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <boost/range/iterator_range.hpp>

class GDBFilter : public Filter
{
  public:
    GDBFilter()
    : Filter    ({Filter::thiourea, RWMOL_SPTR(SmartsToMol("s"))}, {})
    {
        init_filters();
    }

    virtual std::string geom_filter(RWMOL_SPTR mol)
    {
        // Check that molecule has geometry
        MolOps::addHs(*mol); // could move to do_modification()
        if (DGeomHelpers::EmbedMolecule(*mol, 0, -1, true, false, 2.0, true, 1, 0, 0.001, false, true, true, true) == -1)
            return std::string(); // no geometry

        // Check each conformer
        std::vector<ROMol::ConformerIterator>badconfs;
        badconfs.reserve(mol->getNumConformers());
        for (auto it = mol->beginConformers(); it != mol->endConformers(); ++it)
        {
            for (auto& atom : mol->atoms())
            {
                if (atom->getAtomicNum() == 6 and atom->getDegree() == 4)
                {
                    RDGeom::Point3D& center = (*it)->getAtomPos(atom->getIdx());
                    std::vector<RDGeom::Point3D> pos;
                    pos.reserve(4);
                    int j = 0;
                    for (const auto &nbri : boost::make_iterator_range(mol->getAtomNeighbors(atom)))
                    {
                        pos.emplace_back((*it)->getAtomPos((*mol)[nbri]->getIdx()));
                        pos.back() -= center;
                        pos.back().normalize();
                        if (j > 0)
                            pos[j] -= pos[0];
                        ++j;
                    }
                    double valume = std::abs(pos[1].dotProduct(pos[2].crossProduct(pos[3])) / 6.0);
                    if (valume < 0.345)
                    {
                        badconfs.emplace_back(it);
                        break;
                    }
                }
            }
        }
        // If all conformers fail, filter the molecule
        if (badconfs.size() == mol->getNumConformers())
            return "SAV";
        else if (mol->getNumConformers() != 0)
        {
            for (auto& conf : badconfs)
                mol->removeConformer((*conf)->getId());
        }
        // If we're here, then it passed the filter
        MolOps::removeHs(*mol);
        return std::string();
    }
  protected:
    void init_filters()
    {
        NewPatternFilter bredt_violation("bredt violation", RWMOL_SPTR(SmartsToMol("[R]@;=,:[R&x3](@[R])@[R]")), {RWMOL_SPTR(SmartsToMol("[R]@[R&x3](@[R])@[x3,x4]"))});
        Filter::filters.emplace(BREDT_RULE, new NewFilter("Bredt's rule",
            [bredt_violation = bredt_violation] (RWMOL_SPTR mol) -> std::string
            {
                std::vector<std::vector<int>> matches = bredt_violation.filter_with_exceptions(mol);
                if (!matches.empty())
                {
                    bool macrocycle;
                    RingInfo* ring_info = mol->getRingInfo();
                    for (auto& match : matches)
                    {
                        macrocycle = false;
                        for (auto& atom_idx : match)
                        {
                            if (ring_info->minAtomRingSize(atom_idx) >= 8)
                            {
                                macrocycle = true;
                                break;
                            }
                        }
                        if (!macrocycle)
                            return "Bredt's rule";
                    }
                }
                return std::string();
            }
        )); 
        RWMOL_SPTR nitro(SmartsToMol("[N+]([O-])=O"));
        RWMOL_SPTR nitrile(SmartsToMol("C#N"));
        RWMOL_SPTR sulfone(SmartsToMol("O=[S&H0]=O"));
        Filter::filters.emplace(ATOM_COUNT, new NewFilter("atomcounts",
            [nitro = nitro, nitrile = nitrile, sulfone = sulfone] (RWMOL_SPTR mol) -> std::string
            {
                unsigned int nnitros = SubstructMatch(*mol, *nitro).size();
                unsigned int nnitriles = SubstructMatch(*mol, *nitrile).size();
                unsigned int nsulfones = SubstructMatch(*mol, *sulfone).size();
                unsigned int c_num = 0, n_num = 0, o_num = 0, s_num = 0, halogen_num = 0;
                for (auto& atom : mol->atoms())
                {
                    switch (atom->getAtomicNum())
                    {
                        case 6:
                            ++c_num;
                            break;
                        case 7:
                            ++n_num;
                            break;
                        case 8:
                            ++o_num;
                            break;
                        case 16:
                            ++s_num;
                            break;
                        case 9:
                        case 17:
                        case 35:
                        case 53:
                        case 85:
                        case 117:
                            ++halogen_num;
                            break;
                    }
                }
                
                unsigned int carbon = c_num + nsulfones;
                unsigned int nitrogen = n_num - nnitriles - nnitros;
                unsigned int oxygen = o_num - 2 * nsulfones + nnitros;
                if (carbon == 0)
                    return "No carbon atoms";
                if (nitrogen / carbon > 0.6)
                    return "N/C ratio too high";
                if (oxygen / carbon > 0.666)
                    return "O/C ratio too high";
                if ((nitrogen + oxygen) / carbon > 1)
                    return "(N+O)/C ratio too high";
                if (s_num / carbon > 0.333)
                    return "S/C ratio too high";
                if (halogen_num / carbon > 0.5)
                    return "Halogen/C ratio too high";
                return std::string();
            }
        )); 
        Filter::filters.emplace(TRIPLE_BOND_IN_RING, new NewFilter("triple bond in ring",
            [] (RWMOL_SPTR mol) -> std::string
            {
                RingInfo* ring_info = mol->getRingInfo();
                for (auto& bond : mol->bonds())
                {
                    if (ring_info->numBondRings(bond->getIdx()) != 0
                        && ring_info->minBondRingSize(bond->getIdx()) < 9
                        && bond->getBondType() == Bond::BondType::TRIPLE)
                        return "triple bond in ring";
                }
                return std::string();
            }
        )); 
        Filter::filters.emplace(ALLENE, new NewPatternFilter(
            "allene", RWMOL_SPTR(SmartsToMol("[!O]=*=*"))
        )); 
        Filter::filters.emplace(ACIDTAUT, new NewPatternFilter(
            "acidtaut", RWMOL_SPTR(SmartsToMol("C=CO(O)"))
        )); 
        Filter::filters.emplace(AMINAL, new NewPatternFilter(
            "aminal", RWMOL_SPTR(SmartsToMol("[N&X3][C&X4][N,O]"))
        )); 
        Filter::filters.emplace(C_DOUBLE_N, new NewPatternFilter(
            "C=N", RWMOL_SPTR(SmartsToMol("[C&X3]=[N&X2]")), 
            { RWMOL_SPTR(SmartsToMol("[#7]C=N"))
            , RWMOL_SPTR(SmartsToMol("[C&X3]=N[N,O,n,o]"))
            , RWMOL_SPTR(SmartsToMol("cn"))}
        )); 
        Filter::filters.emplace(DECARBOXY_1, new NewPatternFilter(
            "Decarboxy1", RWMOL_SPTR(SmartsToMol("CC(=O)[C&X4]C([O&H1])=O"))
        )); 
        Filter::filters.emplace(DECARBOXY_2, new NewPatternFilter(
            "Decarboxy2", RWMOL_SPTR(SmartsToMol("O=[C&H1][C&X4]C([O&H1])=O"))
        )); 
        Filter::filters.emplace(ENAMINE, new NewPatternFilter(
            "Enamine", RWMOL_SPTR(SmartsToMol("[C&X3]=C[N&X3]")), 
            { RWMOL_SPTR(SmartsToMol("[C&X3]=CNC=[O,N]"))
            , RWMOL_SPTR(SmartsToMol("[O,N]=CC=C[N&X3]"))
            , RWMOL_SPTR(SmartsToMol("[c&X3]c[n,N]"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1C([C,O])=C([C,O])NC([C,O])=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=C([C,O])NC([C,O])=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=[CH1]NC([C,O])=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=C([C,O])NC([C,O])=[CH1]1"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=[CH1]NC([C,O])=[CH1]1"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=[CH1]NC[CH1]=[CH1]1"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=[CH1]NC[CH1]=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1[CH1]=C([C,O])N[CH1]=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1C([C,O])=[CH1]NC([C,O])=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1C([C,O])=[CH1]N[CH1]=C1([C,O])"))
            , RWMOL_SPTR(SmartsToMol("[C,N,O]1C=CNC=C1"))
            , RWMOL_SPTR(SmartsToMol("N#CC=CN"))
            , RWMOL_SPTR(SmartsToMol("C=CNS(=O)=O"))
            , RWMOL_SPTR(SmartsToMol("NC=CS(=O)=O"))
            , RWMOL_SPTR(SmartsToMol("C=CNC=S"))}
        )); 
        Filter::filters.emplace(ENOL, new NewPatternFilter(
            "Enol", RWMOL_SPTR(SmartsToMol("*=C[O&H1]")), 
            { RWMOL_SPTR(SmartsToMol("O=C[O&H1]"))
            , RWMOL_SPTR(SmartsToMol("[a]=[c,C][O&H1]"))
            , RWMOL_SPTR(SmartsToMol("[a]C=C[O&H1]"))}
        )); 
        Filter::filters.emplace(FC_1, new NewPatternFilter(
            "FC1", RWMOL_SPTR(SmartsToMol("F[C&X4][N,O]")), 
            { RWMOL_SPTR(SmartsToMol("Fc([a])=[N,O]")) }
        )); 
        Filter::filters.emplace(FC_2, new NewPatternFilter(
            "FC2", RWMOL_SPTR(SmartsToMol("FC=[N,O]")), 
            { RWMOL_SPTR(SmartsToMol("Fc([a])=[N,O]")) }
        )); 
        Filter::filters.emplace(GEMINAL_1, new NewPatternFilter(
            "Geminal1", RWMOL_SPTR(SmartsToMol("[N&H2,O&H1][C&X4][N&H2,O&H1]"))
        )); 
        Filter::filters.emplace(GEMINAL_2, new NewPatternFilter(
            "Geminal2", RWMOL_SPTR(SmartsToMol("[C&X4][N&H1][C&X4][N&H1][C&X4]"))
        )); 
        Filter::filters.emplace(HEMIACETAL, new NewPatternFilter(
            "Hemiacetal", RWMOL_SPTR(SmartsToMol("[O&H1][C&X4][N,O]"))
        )); 
        Filter::filters.emplace(HEMIAMINAL_1, new NewPatternFilter(
            "Hemiaminal1", RWMOL_SPTR(SmartsToMol("[N&H2][C&X4]O"))
        )); 
        Filter::filters.emplace(HEMIAMINAL_2, new NewPatternFilter(
            "Hemiaminal2", RWMOL_SPTR(SmartsToMol("C[N&H1][C&X4]O"))
        )); 
        Filter::filters.emplace(HEMIAMINAL_3, new NewPatternFilter(
            "Hemiaminal3", RWMOL_SPTR(SmartsToMol("CN(C)[C&X4]O"))
        )); 
        Filter::filters.emplace(HET_HET_OO, new NewPatternFilter(
            "HetHet_OO", RWMOL_SPTR(SmartsToMol("[#8][#8]"))
        ));
        Filter::filters.emplace(HET_HET_NS, new NewPatternFilter(
            "HetHet_NS", RWMOL_SPTR(SmartsToMol("NSO")), 
            { RWMOL_SPTR(SmartsToMol("C=NS[C,S]"))
            , RWMOL_SPTR(SmartsToMol("C=N[S&H1]"))
            , RWMOL_SPTR(SmartsToMol("[N,O]=CNSC"))
            , RWMOL_SPTR(SmartsToMol("[N,O,S]=CN[S&H1]"))
            , RWMOL_SPTR(SmartsToMol("[N+]([O-])=O"))}
        )); 
        Filter::filters.emplace(HET_HET_OS_1, new NewPatternFilter(
            "HetHet_OS1", RWMOL_SPTR(SmartsToMol("C=COSN")), 
            { RWMOL_SPTR(SmartsToMol("S1(=O)O[C,S](=O)**1"))
            , RWMOL_SPTR(SmartsToMol("S1(=O)O[C,S](=O)***1"))}
        )); 
        Filter::filters.emplace(HET_HET_OS_2, new NewPatternFilter(
            "HetHet_OS2", RWMOL_SPTR(SmartsToMol("NS[O&H1]")), 
            { RWMOL_SPTR(SmartsToMol("S1(=O)O[C,S](=O)**1"))
            , RWMOL_SPTR(SmartsToMol("S1(=O)O[C,S](=O)***1"))}
        )); 
        Filter::filters.emplace(HET_HET_OS_3, new NewPatternFilter(
            "HetHet_OS3", RWMOL_SPTR(SmartsToMol("[O&H1]SO[C,O,N]")), 
            { RWMOL_SPTR(SmartsToMol("S1(=O)O[C,S](=O)**1"))
            , RWMOL_SPTR(SmartsToMol("S1(=O)O[C,S](=O)***1"))}
        )); 
        Filter::filters.emplace(HETERO_TRIPLE_BOND, new NewPatternFilter(
            "Hetero triple bond", RWMOL_SPTR(SmartsToMol("[N,O]C#*"))
        ));
        Filter::filters.emplace(HET_HET_NN, new NewPatternFilter(
            "HetHet_NN", RWMOL_SPTR(SmartsToMol("N~N")), 
            { RWMOL_SPTR(SmartsToMol("[#6]N([#6])N=[#6]"))
            , RWMOL_SPTR(SmartsToMol("[#6]=N[N&H1][#6]"))
            , RWMOL_SPTR(SmartsToMol("[c,C]=N[N&H2]"))
            , RWMOL_SPTR(SmartsToMol("[#6]=NN=[#6]"))
            , RWMOL_SPTR(SmartsToMol("NNC=[N,O,S]"))
            , RWMOL_SPTR(SmartsToMol("NNS=O"))}
        )); 
        Filter::filters.emplace(HET_HET_XO, new NewPatternFilter(
            "HetHet_XO", RWMOL_SPTR(SmartsToMol("[F,Cl,Br,I]O"))
        ));
        Filter::filters.emplace(HET_HET_XN, new NewPatternFilter(
            "HetHet_XN", RWMOL_SPTR(SmartsToMol("[F,Cl,Br,I]N"))
        ));
        Filter::filters.emplace(HET_HET_NO, new NewPatternFilter(
            "HetHet_NO", RWMOL_SPTR(SmartsToMol("N-,=O")), 
            { RWMOL_SPTR(SmartsToMol("[#6]=NO[C,S]"))
            , RWMOL_SPTR(SmartsToMol("[#6]=N[O&H1]"))
            , RWMOL_SPTR(SmartsToMol("[N,O,S]=[C,S]NO[C,S]"))
            , RWMOL_SPTR(SmartsToMol("[N,O,S]=CN[O&H]"))
            , RWMOL_SPTR(SmartsToMol("[N+]([O-])=O"))}
        )); 
        Filter::filters.emplace(HET_HET_SS, new NewPatternFilter(
            "HetHet_SS", RWMOL_SPTR(SmartsToMol("[#16][#16]"))
        ));
        Filter::filters.emplace(HET_HET_AROMATIC, new NewPatternFilter(
            "HetHet_Aromatic", RWMOL_SPTR(SmartsToMol("[n,o][n,o][n,o][n,o][n,o]"))
        ));
        Filter::filters.emplace(HET_HET_HET, new NewPatternFilter(
            "HetHetHet", RWMOL_SPTR(SmartsToMol("N[O,N]N"))
        ));
        Filter::filters.emplace(HETERO_SULFUR_WITHOUT_SULFONE, new NewPatternFilter(
            "Hetero-sulfur without sulfone", RWMOL_SPTR(SmartsToMol("*[S&X2][!#6]"))
        ));
        Filter::filters.emplace(HETERO_SR_1, new NewPatternFilter(
            "Hetero-SR1", RWMOL_SPTR(SmartsToMol("[O,N]1C[O,N]1"))
        ));
        Filter::filters.emplace(HETERO_SR_2, new NewPatternFilter(
            "Hetero-SR2", RWMOL_SPTR(SmartsToMol("C1[N,O][N,O]C1"))
        ));
        Filter::filters.emplace(HETERO_SR_3, new NewPatternFilter(
            "Hetero-SR3", RWMOL_SPTR(SmartsToMol("C1[N,O]C[N,O]1"))
        ));
        Filter::filters.emplace(HETERO_SR_4, new NewPatternFilter(
            "Hetero-SR4", RWMOL_SPTR(SmartsToMol("[C,N]=*1***1"))
        ));
        Filter::filters.emplace(HETERO_SR_5, new NewPatternFilter(
            "Hetero-SR5", RWMOL_SPTR(SmartsToMol("*=*1**1"))
        ));
        Filter::filters.emplace(HETERO_SR_6, new NewPatternFilter(
            "Hetero-SR6", RWMOL_SPTR(SmartsToMol("O=C1**C1=O"))
        ));
        Filter::filters.emplace(HETERO_SR_7, new NewPatternFilter(
            "Hetero-SR7", RWMOL_SPTR(SmartsToMol("O=[C,S]1*[C,S](=O)*1"))
        ));
        RWMOL_SPTR intramol_patt11(SmartsToMol("[C&X4][N&H2]"));
        RWMOL_SPTR intramol_patt12(SmartsToMol("[N&X3][N&H2]"));
        RWMOL_SPTR intramol_patt21(SmartsToMol("O=[S,C](C)C"));
        RWMOL_SPTR intramol_patt22(SmartsToMol("O=[S&H1,C&H1]C"));
        Filter::filters.emplace(INTRAMOL, new NewFilter("Intramol",
            [intramol_patt11 = intramol_patt11, intramol_patt12 = intramol_patt12, 
                intramol_patt21 = intramol_patt21, intramol_patt22 = intramol_patt22] (RWMOL_SPTR mol) -> std::string
            {
                MatchVectType result;
                if ( (SubstructMatch(*mol, *intramol_patt11, result)
                     || SubstructMatch(*mol, *intramol_patt12, result))
                    &&
                     (SubstructMatch(*mol, *intramol_patt21, result)
                     || SubstructMatch(*mol, *intramol_patt22, result)))
                    return "Intramol";
                return std::string();
            }
        )); 
        Filter::filters.emplace(MIXED_1, new NewPatternFilter(
            "Mixed1", RWMOL_SPTR(SmartsToMol("C=CON"))
        ));
        Filter::filters.emplace(MIXED_2, new NewPatternFilter(
            "Mixed2", RWMOL_SPTR(SmartsToMol("C=COC(N)=[O,N]"))
        ));
        Filter::filters.emplace(MIXED_3, new NewPatternFilter(
            "Mixed3", RWMOL_SPTR(SmartsToMol("C=COC(O)=N"))
        ));
        Filter::filters.emplace(MIXED_4, new NewPatternFilter(
            "Mixed4", RWMOL_SPTR(SmartsToMol("O=C(N)[O&H1]"))
        ));
        Filter::filters.emplace(MIXED_5, new NewPatternFilter(
            "Mixed5", RWMOL_SPTR(SmartsToMol("O=C(N)[O&H1][C,O,N]"))
        ));
        Filter::filters.emplace(MIXED_6, new NewPatternFilter(
            "Mixed6", RWMOL_SPTR(SmartsToMol("[C,S](=O)O[C,S]=O")), 
            { RWMOL_SPTR(SmartsToMol("O=[C,S]1**[C,S](O1)=O"))
            , RWMOL_SPTR(SmartsToMol("O=[C,S]1***[C,S](O1)=O"))}
        )); 
        Filter::filters.emplace(MIXED_7, new NewPatternFilter(
            "Mixed7", RWMOL_SPTR(SmartsToMol("OCOCO"))
        ));
        Filter::filters.emplace(ORTHO_1, new NewPatternFilter(
            "Ortho1", RWMOL_SPTR(SmartsToMol("[N,O]C(C)(C)N"))
        ));
        Filter::filters.emplace(ORTHO_2, new NewPatternFilter(
            "Ortho2", RWMOL_SPTR(SmartsToMol("[N,O][C&H1](C)N"))
        ));
        Filter::filters.emplace(ORTHOESTER, new NewPatternFilter(
            "Orthoester", RWMOL_SPTR(SmartsToMol("C([N,O])([N,O])([N,O])"))
        ));
        Filter::filters.emplace(TOPO1_33SPIRO, new NewPatternFilter(
            "Topo1-33spiro", RWMOL_SPTR(SmartsToMol("*~1~*~*~1~2~*~*2"))
        ));
        Filter::filters.emplace(TOPO1_33FUSE, new NewPatternFilter(
            "Topo1-33fuse", RWMOL_SPTR(SmartsToMol("*~1~*~2~*1~*2"))
        ));
        Filter::filters.emplace(TOPO1_34SPIRO, new NewPatternFilter(
            "Topo1-34spiro", RWMOL_SPTR(SmartsToMol("*~1~*~*~*~1~2~C~*2"))
        ));
        Filter::filters.emplace(TOPO1_34FUSE, new NewPatternFilter(
            "Topo1-34fuse", RWMOL_SPTR(SmartsToMol("*~1~*~*~2~*1~*2"))
        ));
        Filter::filters.emplace(TOPO1_44SPIRO, new NewPatternFilter(
            "Topo1-44spiro", RWMOL_SPTR(SmartsToMol("*~1~*~*~*~1~2~*~*~*2"))
        ));
        Filter::filters.emplace(TOPO1_44FUSE, new NewPatternFilter(
            "Topo1-44fuse", RWMOL_SPTR(SmartsToMol("*~1~*~*~2*~1~*~*2"))
        ));
        Filter::filters.emplace(TOPO1_44BRIDGE, new NewPatternFilter(
            "Topo1-44bridge", RWMOL_SPTR(SmartsToMol("*~1~*~2~*~*~1~*~2"))
        ));
        Filter::filters.emplace(NON_AROMATIC_NITRO, new NewPatternFilter(
            "Non-aromatic nitro", RWMOL_SPTR(SmartsToMol("[!a][N+](=O)[O-]"))
        ));
        Filter::filters.emplace(NON_AROMATIC_HALOGEN, new NewPatternFilter(
            "Non-aromatic halogen", RWMOL_SPTR(SmartsToMol("[!c][F,Cl,Br,I]"))
        ));
        Filter::filters.emplace(HETERO_AROMATIC_HALOGEN, new NewPatternFilter(
            "Hetero-aromatic halogen", RWMOL_SPTR(SmartsToMol("[!c]c[F,Cl,Br,I]"))
        ));
        Filter::filters.emplace(DOUBLE_BOND_IN_4_RING, new NewPatternFilter(
            "Double bond in 4 ring", RWMOL_SPTR(SmartsToMol("*~1~*=,#,:*~*1")), 
            { RWMOL_SPTR(SmartsToMol("[A]1[A][a]~[a]1")) }
        )); 
        Filter::filters.emplace(DOUBLE_BOND_IN_3_RING, new NewPatternFilter(
            "Double bond in 3 ring", RWMOL_SPTR(SmartsToMol("*~1~*=,#,:*1"))
        ));
        Filter::filters.emplace(POLYCYCLIC_1, new NewPatternFilter(
            "Polycyclic1", RWMOL_SPTR(SmartsToMol("*@N(@*)@[O,N]"))
        ));
        Filter::filters.emplace(POLYCYCLIC_2, new NewPatternFilter(
            "Polycyclic2", RWMOL_SPTR(SmartsToMol("*@N(@*)@C@=*"))
        ));
        Filter::filters.emplace(BETA_KETO_CARBOXYL, new NewPatternFilter(
            "Beta keto carboxyl", RWMOL_SPTR(SmartsToMol("[O&H1]C(=O)[C&X4]C(=O)[!N]"))
        ));
    }
};