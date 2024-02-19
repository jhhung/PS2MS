# Drug detection
- The final step of PS<sup>2</sup>MS is to compare the analyte and the synthetic database
- We utilize the data type from the [cfm-id-code](https://bitbucket.org/wishartlab/cfm-id-code/src/master/cfm/) project

## Usage

### Build project
Follow the step in INSTALL.txt

### Run detection

There are two executable files. One is for verification, which provides rankings and related scores if the object to be tested is known. The other is for testing, used when the identity of the test object is unknown, to obtain results.

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

