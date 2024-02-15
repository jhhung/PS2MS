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
- PS<sup>2</sup>MS build the synthetic database by enumerating possible derivatives.

### Usage:

#### Install conda env 
```bash
conda env create --name drug_detection -f drug-detection/environment.yml
```

#### Build project
```
cd drug-detection
conda activate drug_detection

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
```

#### Run permutation
```bash
cd drug-detection
conda activate drug_detection

./build/build_database <split_count> <split_idx> <output_dir>

# for example: ./build/build_database 200 0 /tmp
```


## NEIMS
- PS<sup>2</sup>MS employs NEIMS to predict the mass spectrum of the compounds in synthetic database
- Original paper: ["Rapid Prediction of Electron–Ionization Mass Spectrometry Using Neural Networks"](https://pubs.acs.org/doi/10.1021/acscentsci.9b00085)
- Source code: [NEIMS Github](https://github.com/brain-research/deep-molecular-massspec/issues)
- Quickstart of retraining model from NEIMS: [Model retrain quickstart](https://github.com/brain-research/deep-molecular-massspec/blob/main/Model_Retrain_Quickstart.md)

### Usage:

#### Prepare dataset
```bash
python3 make_train_test_split.py --main_sdf_name=path/to/mainlib_merge.SDF --replicates_sdf_name=path/to/test.SDF --output_master_dir=path/to/output/dir/spectra_tf_records
```
Arguments: 
* Required
  * `--main_sdf_name`: The file comprises mass spectra of all compounds, encompassing both the training and test sets, presented in SDF format.
  * `--replicates_sdf_name`: The file comprises mass spectra of compounds of test set, presented in SDF format.
  * `--output_master_dir`: The output directory of preprocessing data.

#### Training
```bash
python3 molecule_estimator.py \
  --dataset_config_file=path/to/output/dir/spectra_tf_records/query_replicates_val_predicted_replicates_val.json \
  --train_steps=10000 \
  --model_dir=path/of/output/model/models/output \
  --hparams=make_spectra_plots=True,mask=False,mass_power=0 --alsologtostderr
```
Arguments:
* Required
  * `--dataset_config_file`: The path of preprocessed dataset
  * `--model_dir`: The path of model

#### Predicting
```bash
python3 make_spectra_prediction.py \
  --input_file=path/to/test_SMILES.txt \
  --output_file=path/to/test.SDF \
  --weights_dir=path/of/output/model/models/output
```
Arguments:
* Required
  * `--input_file`: The SMILES of test set, presented in txt format.
  * `--output_file`: The predict mass spectrum of test set, presented in SDF format.
  * `--weights_dir`: The path of model.

#### Transfer file
```bash
python3 sdf_to_msp.py ${input_file} ${is_predict_spectrum}
```
Arguments:
* Required
  * Convert files from SDF format to msp format.
  * `${is_predict_spectrum}`: a boolean value. If the input file contains predict mass spectrum, this value should beset to True.

## DeepEI
- PS<sup>2</sup>MS employs DeepEI to predict the fingerprint of the unknown analyte. 
- Original paper: ["Predicting a Molecular Fingerprint from an Electron Ionization Mass Spectrum with Deep Neural Networks"](https://pubs.acs.org/doi/10.1021/acs.analchem.0c01450)
- Source code: [DeepEI Github](https://github.com/hcji/DeepEI)

### Usage:

#### Install conda env 
```bash
conda env create --name deepei -f DeepEI/environment.yml
```

#### Run prediction
```bash
cd DeepEI
conda activate deepei

python predict.py --help
```

#### Combine the mass spectrum and fingerprint
```bash
python3 merge_fp_into_msp.py \
  ${MS.msp} \
  ${FP.pkl} \
  ${result.msp}
```


## Drug detection
- The final step of PS<sup>2</sup>MS to compare the analyte and the synthetic database
- We utilize the data type from the [cfm-id-code](https://bitbucket.org/wishartlab/cfm-id-code/src/master/cfm/) project

### Usage

#### Build project
Follow the step in INSTALL.txt

#### Run detection
* Condition 1: The analytes are known.
  * `database_file`: The file contains the mass spectrum and the fingerprint of enumerated compounds and is in msp format.
  * `test_file`: The file contains the mass spectrum and the fingerprint of the analytes and is in msp format.
  * `result.txt`: The file records the ranking performance of every analyte. For every analyte, the first line is the rank and the SMSF score of the answer in enumerate database. The following lines are the SMILES and SMSF score of molecular which has higher rank than the answer in enumerate database.
  * `scores.txt`: The file records the similarity scores of every analyte. For an analyte, it records the cosine similarity of mass spectrum, the Jaccard similarity of mass spectrum, the similarity of fingerprint and the SMSF score of with the corresponding one in enumerate database.
  * `restrict_mw`: A boolean value. If true, the system will solely compute the similarity for compounds whose molecular weight falls within the range of the analyte's molecular weight plus or minus one..
  * `top_n_of_JS_of_MS`: A positive integer. This value determines the number of the highest peaks used to compute the Jaccard similarity of mass spectra between two compounds. If set to zero, the system will not calculate the Jaccard similarity of mass sepctrum.
```bash
./build/drug-comparation/cfm_comparation \
  ${database_file} \
  ${test_file} \
  ${result.txt} \
  ${scores.txt} \
  ${restrict_mw} ${top_n_of_JS_of_MS}
```

* Condition 2: The analytes are unknown.
  * `database_file`: The file contains the mass spectrum and the fingerprint of enumerated compounds and is in msp format.
  * `test_file`: The file contains the mass spectrum and the fingerprint of the analytes and is in msp format.
  * `result.txt`: The file records the result of every analyte. For every analyte, the top 100 compounds is recorded with the SMILES and SMSF score.
  * `restrict_mw`: A boolean value. If true, the system will solely compute the similarity for compounds whose molecular weight falls within the range of the analyte's molecular weight plus or minus one..
  * `top_n_of_JS_of_MS`: A positive integer. This value determines the number of the highest peaks used to compute the Jaccard similarity of mass spectra between two compounds. If set to zero, the system will not calculate the Jaccard similarity of mass sepctrum.
```bash
./build/drug-comparation/cfm_comparation_without_answer \
  ${database_file} \
  ${test_file} \
  ${result.txt} \
  ${restrict_mw} ${top_n_of_JS_of_MS}
```

