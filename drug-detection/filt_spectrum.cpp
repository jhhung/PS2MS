#include <iostream>
#include <fstream>
#include <string.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include "DruglikeFilter.hpp"
using namespace RDKit;

using FILTER_SPTR = boost::shared_ptr<Filter>;

void help_and_exit()
{
    std::cerr << "Usage: standalone_filter <smiles_filename> <output_filename> <enable_druglike_filter>\n\n";
    std::cerr << "Args:\n";
    std::cerr << "\tsmiles_filename: Molecules to filter. Each line containes only one SMILES (target molecule).\n";
    std::cerr << "\toutput_filename: Molecules that passed filter will write into this file.\n";
    std::cerr << "\tenable_druglike_filter (Optional): Whether to enable druglike filter. 0: disable, 1:enable (default)\n\n";
    exit(-1);
}

int main(int argc, char **argv)
{
    if (argc < 3 || argc > 4)
        help_and_exit();
    
    std::ifstream infile(argv[1]);
    std::ofstream outfile(argv[2]);
    FILTER_SPTR filter;
    if (!infile)
    {
        std::cerr << "ERROR: " << argv[1] << " doesn't exist.\n\n";
        help_and_exit();
    }
    else if (argc == 4)
    {
        if (strcmp(argv[3], "0") == 0)
            filter.reset(new SMUFilter());
        else if (strcmp(argv[3], "1") == 0)
            filter.reset(new DruglikeFilter());
        else
            help_and_exit();
    }
    else
        filter.reset(new DruglikeFilter());
    
    std::string row, smiles, reason;
    RWMOL_SPTR mol;
    infile >> row; // skip header
    while (infile >> row)
    {
        smiles = row.substr(0, row.find(","));
        try
        {
            mol.reset(SmilesToMol(smiles));
            if (mol)
            {
                reason = (*filter)(mol);
                if (reason.empty())
                    outfile << row << "\n";
                else
                    std::cout << smiles << ": " << reason << "\n";
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
    }
}