# PS<sup>2</sup>MS

Repository for NYCU JHHLab NPS detection project. 

The repository for the paper ["PS<sup>2</sup>MS"](): A Deep Learning-Based Prediction System for Detecting Novel New Psychoactive Substances Using Mass Spectrometry.

## Informations

PS<sup>2</sup>MS is designed specifically to address the limitations of identifying the emergence of unidentified novel illicit drugs. 
PS<sup>2</sup>MS builds a synthetic NPS database by enumerating possible derivatives based on the core structure of a preselected illicit drug.
The system leverages two deep learning tools, NEIMS and DeepEI, to generate mass spectra and chemical fingerprints, respectively.
Finally, PS<sup>2</sup>MS calculates the integrated similarity scores(SMSF) between the unknown analyte and the derivatives from synthetic database and yields a list of potential NPS identities for the analyte.


## Requirement

- GNU [g++-10](https://gcc.gnu.org/gcc-10/) or higher
- [CMake 3.16.0](https://cmake.org/download/) or higher to build the enumeration step and the detection step
- [python 3.6.9](https://www.python.org/downloads/) or higher to the run deep learing tools
- [rdkit](https://www.rdkit.org/docs/Install.html)
  - build the c++ code from the source and install python package from conda


## Enumeration
- PS<sup>2</sup>MS will generate a synthetic database by substituting hydrogen atoms in the core structure of a given drug with functional groups.
- The synthetic database will conduct matching calculations with suspicious analytes using mass spectrometry and chemical fingerprints.
- See how to build and run the enumeration [here](enumeration/README.md)


## NEIMS
- PS<sup>2</sup>MS employs NEIMS to predict the mass spectrum of the compounds in synthetic database
- The system uses NEIMS to predict the spectrum of compounds of synthetic database.
- See how to train the NEIMS [here](neims/README.md)


## DeepEI
- PS<sup>2</sup>MS employs DeepEI to predict the fingerprint of the unknown analyte.
- The system also use the function of DeepEI to calculate the fingerprint of compounds of synthetic database.
- The system merges the fingerprint into the msp file.
- See how to predict fingerprint [here](DeepEI/README.md)


## Drug detection
- The final step of PS<sup>2</sup>MS is to compare the analyte and the synthetic database.
- The system will compare the spectrum and chemical fingerprint between compounds and generate a list of the hundred most similar compounds which are ranked by similarity score.
- See how to build and run drug detection [here](drug-detection/README.md)

