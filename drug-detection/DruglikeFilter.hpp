#include "SMUFilter.hpp"
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/Lipinski.h>

class DruglikeFilter : public SMUFilter
{
  public:
    DruglikeFilter() : SMUFilter()
    {
        Filter::no_sulfone.emplace_back(SmartsToMol("C1SC(=S)NC1=O"));
        Filter::no_sulfone.emplace_back(SmartsToMol("O=C(*)S*"));
        Filter::no_sulfone.emplace_back(SmartsToMol("O=C(*)OS*"));
        Filter::no_sulfone.emplace_back(SmartsToMol("*S[Cl,F,I,Br]"));
        init_filters();
    }
  protected:
    void init_filters()
    {
        Filter::filters.emplace(NO_N_OR_O, new NewFilter("No nitrogen or oxygen",
            [] (RWMOL_SPTR mol) -> std::string
            {
                bool found = false;
                int num;
                for (auto& atom : mol->atoms())
                {
                    num = atom->getAtomicNum();
                    if (num == 7 || num == 8)
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    return "No nitrogen or oxygen";
                return std::string();
            }
        )); 
        Filter::filters.emplace(TOO_MANY_HALOGEN, new NewFilter("More than 5 halogens",
            [] (RWMOL_SPTR mol) -> std::string
            {
                int total = 0;
                for (auto& atom : mol->atoms())
                {
                    switch (atom->getAtomicNum())
                    {
                        case 9:
                        case 17:
                        case 35:
                        case 53:
                        case 85:
                        case 117:
                            total += 1;
                            break;
                    }
                    if (total > 5)
                        break;
                }
                if (total > 5)
                    return "More than 5 halogens";
                return std::string();
            }
        )); 
        RWMOL_SPTR aldehyde(SmartsToMol("O=[CH]"));
        Filter::filters.emplace(TOO_MANY_ALDEHYDE, new NewFilter("More than 2 aldehydes",
            [aldehyde = aldehyde] (RWMOL_SPTR mol) -> std::string
            {
                if (SubstructMatch(*mol, *aldehyde).size() > 2)
                    return "More than 2 aldehydes";
                return std::string();
            }
        )); 
        RWMOL_SPTR terminal_2_bond_c(SmartsToMol("O=[CH]"));
        Filter::filters.emplace(TOO_MANY_TERMINAL_DOUBLE_BOND_C, new NewFilter("More than one terminal double-bonded carbon",
            [terminal_2_bond_c = terminal_2_bond_c] (RWMOL_SPTR mol) -> std::string
            {
                if (SubstructMatch(*mol, *terminal_2_bond_c).size() > 1)
                    return "More than one terminal double-bonded carbon";
                return std::string();
            }
        )); 
        RWMOL_SPTR enes(SmartsToMol("[A;!O;!S]=[A;!O;!S]"));
        Filter::filters.emplace(TOO_MANY_ENE, new NewFilter("Too many enes",
            [enes = enes] (RWMOL_SPTR mol) -> std::string
            {
                if (SubstructMatch(*mol, *enes).size() > 3)
                    return "Too many enes";
                return std::string();
            }
        )); 
        Filter::filters.emplace(TOO_SMALL, new NewFilter("MW<100 Dalton",
            [] (RWMOL_SPTR mol) -> std::string
            {
                if (calculate_weight(mol) < 100)
                    return "MW<100 Dalton";
                return std::string();
            }
        )); 
        RWMOL_SPTR amine1(SmartsToMol("[N+0&H2]"));
        RWMOL_SPTR amine2(SmartsToMol("[A][N+0&H1][A]"));
        RWMOL_SPTR amine3(SmartsToMol("[A][N+0&H0]([A])[A]"));
        Filter::filters.emplace(TOO_MANY_AMINE, new NewFilter(">3 basic amines",
            [amine1 = amine1, amine2 = amine2, amine3 = amine3] (RWMOL_SPTR mol) -> std::string
            {
                if (SubstructMatch(*mol, *amine1).size()
                    + SubstructMatch(*mol, *amine2).size()
                    + SubstructMatch(*mol, *amine3).size() > 3)
                    return ">3 basic amines";
                return std::string();
            }
        )); 
        Filter::filters.emplace(DRUG_LIKENESS, new NewFilter("Drug likeness paramters",
            [] (RWMOL_SPTR mol) -> std::string
            {
                if (Descriptors::calcClogP(*mol) > 7)
                    return "CLogP > 7";
                if (Descriptors::calcNumHBA(*mol) > 10)
                    return "HBA > 10";
                if (Descriptors::calcNumHBD(*mol) > 5)
                    return "HBA > 5";
                if (Descriptors::calcNumRotatableBonds(*mol) > 11)
                    return "Rotatable bonds > 11";
                return std::string();
            }
        )); 
        Filter::filters.emplace(HET_SULFONE_HET, new NewPatternFilter(
            "het_sulfone_het", RWMOL_SPTR(SmartsToMol("[!#6]S(=O)(=O)[!#6]")), 
            { RWMOL_SPTR(SmartsToMol("[#7]S(=O)(=O)[#7]")) }
        )); 
        Filter::filters.emplace(ENOL_ETHER, new NewPatternFilter(
            "Enol ether", RWMOL_SPTR(SmartsToMol("C=C[OH0]"))
        )); 
        Filter::filters[HETERO_TRIPLE_BOND] = Filter::NEW_FILTER_SPTR(new NewPatternFilter(
            "Hetero triple bond", RWMOL_SPTR(SmartsToMol("[!#6]C#*"))
        )); 
        Filter::filters.emplace(PERHALO_KETONE, new NewPatternFilter(
            "perhalo-ketone", RWMOL_SPTR(SmartsToMol("[F,Cl,Br,I][c,C]([F,Cl,Br,I])([F,Cl,Br,I])[c,C]=O[c,C]"))
        )); 
        static_cast<NewPatternFilter*>(Filter::filters[TERMINAL_SULFUR_2].get())->set_exceptions({
            Filter::thiourea
          , RWMOL_SPTR(SmartsToMol("C1SC(=S)NC1=O"))
        });
        Filter::filters.emplace(POLYCYCLIC_SULFUR, new NewPatternFilter(
            "Polycyclic_sulfur", RWMOL_SPTR(SmartsToMol("*@N(@S)@*"))
        ));
        Filter::filters.emplace(HALOPYRIMIDINE, new NewPatternFilter(
            "Halopyrimidine", RWMOL_SPTR(SmartsToMol("c1cnc([F,Cl,Br,I])nc1"))
        ));
        Filter::filters[HETERO_AROMATIC_HALOGEN] = Filter::NEW_FILTER_SPTR(new NewPatternFilter(
            "Hetero-aromatic halogen", RWMOL_SPTR(SmartsToMol("[n,o]c[F,Cl,Br,I]"))
        ));
        static_cast<NewPatternFilter*>(Filter::filters[NON_AROMATIC_HALOGEN].get())->set_exceptions({
            RWMOL_SPTR(SmartsToMol("cC(F)(F)F"))
          , RWMOL_SPTR(SmartsToMol("cC([Cl])([Cl])[Cl]"))
          , RWMOL_SPTR(SmartsToMol("cC([Cl])([Br])[Br]"))
        });
        Filter::filters.emplace(ANHYDRIDE, new NewPatternFilter(
            "Anhydride", RWMOL_SPTR(SmartsToMol("O=[C&H0]O[C&H0]=O"))
        ));
        Filter::filters.emplace(BETA_HATERO_CARBONYL, new NewPatternFilter(
            "Beta-heterosubstituted carbonyl", RWMOL_SPTR(SmartsToMol("*C(=O)[C&H2]C([!C])C")),
            { RWMOL_SPTR(SmartsToMol("*C(=O)[C&H2]C([!C]*)C")) }
        ));
        Filter::filters.emplace(LONGALKYL, new NewPatternFilter(
            "longalkyl", RWMOL_SPTR(SmartsToMol("[C&H2][C&H2][C&H2][C&H2][C&H2][C&H2][C&H3]"))
        ));
        Filter::filters.emplace(TERMINAL_DOUBLE_BOND, new NewPatternFilter(
            "Terminal double bond", RWMOL_SPTR(SmartsToMol("[!S;!O;!C]=[!S;!O]*")),
            { RWMOL_SPTR(SmartsToMol("**=**")) }
        ));
        Filter::filters[HET_HET_AROMATIC] = Filter::NEW_FILTER_SPTR(new NewPatternFilter(
            "HetHet_Aromatic", RWMOL_SPTR(SmartsToMol("n[o,s]")),
            { RWMOL_SPTR(SmartsToMol("o1nc[c,n]c1")) }
        ));
        Filter::filters.emplace(HET_HET_AROMATIC_2, new NewPatternFilter(
            "HetHet_Aromatic2", RWMOL_SPTR(SmartsToMol("nnn")),
            { RWMOL_SPTR(SmartsToMol("n1nn*c1")) }
        ));
    }
};
