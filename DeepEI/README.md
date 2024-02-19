# DeepEI
- PS<sup>2</sup>MS employs DeepEI to predict the fingerprint of the unknown analyte. 
- Original paper: ["Predicting a Molecular Fingerprint from an Electron Ionization Mass Spectrum with Deep Neural Networks"](https://pubs.acs.org/doi/10.1021/acs.analchem.0c01450)
- Source code: [DeepEI Github](https://github.com/hcji/DeepEI)
- Original markdown: [DeepEI ReadMe.md](https://github.com/hcji/DeepEI/blob/master/ReadMe.md)

## Usage:

### Install conda env 
```bash
conda env create --name deepei -f DeepEI/environment.yml
```

### Run prediction
```bash
cd DeepEI
conda activate deepei

python predict.py --help
```

### Combine the mass spectrum and fingerprint
```bash
python3 merge_fp_into_msp.py \
  ${MS.msp} \
  ${FP.pkl} \
  ${result.msp}
```


