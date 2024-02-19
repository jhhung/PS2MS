# NEIMS
- PS<sup>2</sup>MS employs NEIMS to predict the mass spectrum of the compounds in synthetic database
- Original paper: ["Rapid Prediction of Electron–Ionization Mass Spectrometry Using Neural Networks"](https://pubs.acs.org/doi/10.1021/acscentsci.9b00085)
- Source code: [NEIMS Github](https://github.com/brain-research/deep-molecular-massspec/issues)
- Quickstart of retraining model from NEIMS: [Model retrain quickstart](https://github.com/brain-research/deep-molecular-massspec/blob/main/Model_Retrain_Quickstart.md)

## Usage:

### Prepare dataset
```bash
python3 make_train_test_split.py --main_sdf_name=path/to/mainlib_merge.SDF --replicates_sdf_name=path/to/test.SDF --output_master_dir=path/to/output/dir/spectra_tf_records
```
Arguments: 
* Required
  * `--main_sdf_name`: The file comprises mass spectra of all compounds, encompassing both the training and test sets, presented in SDF format.
  * `--replicates_sdf_name`: The file comprises mass spectra of compounds of test set, presented in SDF format.
  * `--output_master_dir`: The output directory of preprocessing data.

### Training
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

### Predicting
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

### Transfer file
```bash
python3 sdf_to_msp.py ${input_file} ${is_predict_spectrum}
```
Arguments:
* Required
  * Convert files from SDF format to msp format.
  * `${is_predict_spectrum}`: a boolean value. If the input file contains predict mass spectrum, this value should beset to True.


