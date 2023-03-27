r"""Makes spectra prediction using model and writes predictions to SDF.

Make predictions using our trained model. Example of how to run:

# Save weights to a models directory
$ MODEL_WEIGHTS_DIR=/tmp/neims_model
$ cd $MODEL_WEIGHTS_DIR
$ wget https://storage.googleapis.com/deep-molecular-massspec/massspec_weights/massspec_weights.zip  # pylint: disable=line-too-long
$ unzip massspec_weights.zip

$ python make_spectra_prediction.py \
--input_file=examples/pentachlorobenzene.sdf \
--output_file=/tmp/neims_model/annotated.sdf \
--weights_dir=$MODEL_WEIGHTS_DIR/massspec_weights
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

from absl import app
from absl import flags
from absl import logging

import spectra_predictor
from rdkit import Chem
from rdkit.Chem import Descriptors

FLAGS = flags.FLAGS
flags.DEFINE_string('input_file', 'input.sdf',
                    'Name of input file for predictions.')
flags.DEFINE_string('weights_dir',
                    '/usr/local/massspec_weights',
                    'Name of directory that stores model weights.')
flags.DEFINE_string('output_file', 'annotated.sdf',
                    'Name of output file for predictions.')


def main(_):
  logging.info('Loading weights from %s', FLAGS.weights_dir)
  predictor = spectra_predictor.NeimsSpectraPredictor(
      model_checkpoint_dir=FLAGS.weights_dir, hparams_str=open(os.path.join(FLAGS.weights_dir, "command_line_arguments.txt")).read()[39:])

  logging.info('Loading molecules from %s', FLAGS.input_file)
  # mols_from_file = spectra_predictor.get_mol_list_from_sdf(
  #     FLAGS.input_file, predictor._hparams)
  mols_from_file = spectra_predictor.get_mol_list_from_smiles(FLAGS.input_file)
  # mols_from_file = spectra_predictor.get_mol_list_from_csv(FLAGS.input_file)
  if not mols_from_file[0].HasProp("SMILES"):
    for mol in mols_from_file:
        mol.SetProp("SMILES", Chem.MolToSmiles(mol))
  if not mols_from_file[0].HasProp("EXACT MASS"):
    for mol in mols_from_file:
        mol.SetProp("EXACT MASS", str(Descriptors.ExactMolWt(mol)))
  fingerprints, mol_weights, masks = predictor.get_inputs_for_model_from_mol_list(
      mols_from_file)

  logging.info('Making predictions ...')
  spectra_predictions = predictor.make_spectra_prediction(
      fingerprints, mol_weights, masks)

  logging.info('Updating molecules in place with predictions.')
  spectra_predictor.update_mols_with_spectra(mols_from_file,
                                             spectra_predictions)

  logging.info('Writing predictions to %s', FLAGS.output_file)
  with open(FLAGS.output_file, 'w') as f:
    spectra_predictor.write_rdkit_mols_to_sdf(mols_from_file, f)

if __name__ == '__main__':
  app.run(main)
