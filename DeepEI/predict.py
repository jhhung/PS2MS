import os
import pandas as pd
import numpy as np
from io import StringIO as sio
from scipy.sparse import load_npz, csr_matrix
import tensorflow as tf
from tensorflow.keras.models import model_from_json, load_model
from tqdm import tqdm
import pathlib
import argparse

def array2bits(arr):
    return ''.join(['1' if pos else '0' for pos in arr.tolist()])

def bits2array(bits):
    return np.array(list(map(int, bits)), dtype=np.int64)

def parse_sdf(sdf_file):
    from rdkit import Chem
    spr = Chem.SDMolSupplier(sdf_file)
    smiles_list = []
    spectrum = []
    for mol in spr:
        smiles = mol.GetProp('SMILES')
        smiles_list.append(smiles)
        peaks = mol.GetProp('MASS SPECTRAL PEAKS')
        df_peaks = pd.read_csv(sio(peaks), header=None, sep=' ')
        # print(df_peaks.dtypes) # debug
        tmp = np.zeros((2000,))
        tmp[df_peaks[0]] = df_peaks[1]
        spectrum.append(tmp)
    spectrum = np.vstack(spectrum)
    print(spectrum.shape)
    assert spectrum.shape[0] == len(smiles_list)

    return spectrum, smiles_list

def parse_msp(msp_file):
    smiles_list = []
    spec_list = []
    with open(msp_file, 'r') as f:
        while line := f.readline():
            if "ID: " in line:
                smiles_list.append(line.strip().split()[1])
            if "Num Peaks:" in line:
                rows = []
                while row := f.readline():
                    if row.strip() != "":
                        rows.append(row)
                    else:
                        break
                df = pd.read_csv(sio("".join(rows)), header=None, sep=' ')
                print(df.shape)
                vec = np.zeros((2000,))
                vec[df[0]] = df[1]
                spec_list.append(vec)
    
    spec_list = np.vstack(spec_list)
    assert len(smiles_list) == spec_list.shape[0]
    
    return spec_list, smiles_list

def predict_fingerprint(spec, weight_dir: pathlib.Path):
    files = os.listdir(weight_dir / 'weight')
    rfp = np.array([int(f.split('.')[0]) for f in files if '.h5' in f])
    rfp = np.sort(rfp)
    assert len(rfp) == 646
    
    # rfp = set(rfp).intersection(set(fpkeep))
    # rfp = np.sort(list(rfp)).astype(int)
    
    files = [str(f) + '.h5' for f in rfp]
    modjs = open(weight_dir / 'model.json', 'r').read()
    model = model_from_json(modjs)
    pred_fp = np.zeros((spec.shape[0], len(files)))
    for i, f in enumerate(tqdm(files)):
        model.load_weights(weight_dir / 'weight' / f)
        pred = np.round(model.predict(spec))[:,0]
        pred_fp[:,i] = pred  
    return pred_fp

def predict(file_name, checkpoint, save_path: pathlib.Path, file_type='npz', smiles_file=None):
    if file_type == 'sdf':
        spectrum, smiles_list = parse_sdf(sdf_file=file_name)
    elif file_type == 'msp':
        spectrum, smiles_list = parse_msp(msp_file=file_name)
    elif file_type == 'npz':
        spectrum: csr_matrix = load_npz(file_name)
        spectrum = spectrum.toarray()
        smiles_list = open(smiles_file, 'r').read().strip().split('\n')
        assert spectrum.shape[0] == len(smiles_list), f"data size not matched, {spectrum.shape[0]} != {len(smiles_list)}"
    else:
        raise NotImplementedError
    
    pred_fp = predict_fingerprint(spectrum, pathlib.Path(checkpoint))
    
    pred_fp_list = []
    for idx in range(pred_fp.shape[0]):
        pred_fp_list.append(array2bits(pred_fp[idx]))
    df_output = pd.DataFrame({
        'smiles': smiles_list,
        'fingerprint': pred_fp_list
    })
    save_path.parent.mkdir(parents=True, exist_ok=True)
    df_output.to_csv(save_path, index=False)

def config_gpu_mem(mem_size=1024):
    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        try:
            tf.config.set_logical_device_configuration(
                gpus[0],
                [tf.config.LogicalDeviceConfiguration(memory_limit=mem_size)])
            logical_gpus = tf.config.list_logical_devices('GPU')
            print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
        except RuntimeError as e:
            print(e)

def get_args():
    parser = argparse.ArgumentParser("DeepEI prediction")
    parser.add_argument("-i", "--input_file", type=pathlib.Path, required=True, help="input file path")
    parser.add_argument("-s", "--smiles_file", type=pathlib.Path, help="if file_type is npz, smiles file is required", default=None)
    parser.add_argument("-f", "--file_type", type=str, choices=['sdf', 'msp', 'npz', 'auto'], default='auto')
    parser.add_argument("-c", "--checkpoint", type=pathlib.Path, required=True, help="dir to the checkpoint")
    parser.add_argument("-o", "--output_file", type=pathlib.Path, required=True, help="output file path")
    parser.add_argument("-m", "--mem_size", type=int, default=1024, help="gpu memory size")
    args = parser.parse_args()
    assert args.input_file.exists(), "input file not exists"

    if args.file_type == 'auto':
        if args.input_file.suffix.upper() == '.SDF':
            args.file_type = 'sdf'
        elif args.input_file.suffix == '.msp':
            args.file_type = 'msp'
        elif args.input_file.suffix == '.npz':
            args.file_type = 'npz'
            assert args.smiles_file.exists(), "smiles file is required when file_type is npz"
        else:
            raise NotImplementedError("cant automatically recognize the file type")
    return args

if __name__ == '__main__':
    args = get_args()
    config_gpu_mem(mem_size=args.mem_size)
    predict(args.input_file, args.checkpoint, args.output_file, args.file_type, args.smiles_file)
