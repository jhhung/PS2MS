# NPS detection

Repository for NYCU JHHLab NPS detection project

## DeepEI
- Original paper: ["Predicting a Molecular Fingerprint from an Electron Ionization Mass Spectrum with Deep Neural Networks"](https://pubs.acs.org/doi/10.1021/acs.analchem.0c01450)
- Source code: [DeepEI Github](https://github.com/hcji/DeepEI)
### Usage:
- Install conda env 
```bash
conda env create --name deepei -f DeepEI/environment.yml
```
- Run prediction
```bash
cd DeepEI
conda activate deepei

python predict.py --help
```

## NEIMS
- Original paper: ["Rapid Prediction of Electron–Ionization Mass Spectrometry Using Neural Networks"](https://pubs.acs.org/doi/10.1021/acscentsci.9b00085)
- Source code: [NEIMS Github](https://github.com/brain-research/deep-molecular-massspec/issues)
### Usage:
- Install conda env

## Drug-detection
- By Samuel

### Usage:
```bash
cd drug-detection
./build/build_database <split_count> <split_idx> <output_dir>

# for example: ./build/build_database 1 0 /tmp
```