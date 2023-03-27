import os
from rdkit import Chem
import pandas as pd
import numpy as np
from io import StringIO as sio
from scipy.sparse import load_npz, csr_matrix
from model.model import Simple
import torch

device = 'cuda' if torch.cuda.is_available() else 'cpu'

def array2bits(arr):
    return ''.join(['1' if pos else '0' for pos in arr.tolist()])

def parse_sdf(sdf_file):
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
    # TODO:
    
    return None, None

def load_model(checkpoint, ModelType=Simple, **kwarg):
    model = ModelType(
        **kwarg
    )
    state_dict = torch.load(checkpoint, map_location='cpu')['state_dict']
    model.load_state_dict(state_dict)
    model = model.to(device)
    model.eval()
    return model

def predict(file_name, checkpoint, save_path, file_type='sdf'):
    if file_type == 'sdf':
        spectrum, smiles_list = parse_sdf(sdf_file=file_name)
    elif file_type == 'msp':
        spectrum, smiles_list = parse_msp(msp_file=file_name)
    elif file_type == 'npz':
        spectrum: csr_matrix = load_npz(file_name + "_spec.npz")
        spectrum = spectrum.toarray()
        smiles_list = open(file_name + "_smiles.txt", 'r').read().strip().split('\n')
    else:
        raise NotImplementedError
    
    model = load_model(
        checkpoint=checkpoint,
        ModelType=Simple,
        input_dim=2000,
        output_dim=646,
        layers=6,
        res_width=512,
        dropout=0.2
    )

    with torch.no_grad():
        pred_fp = torch.sigmoid(model(torch.tensor(spectrum, device=device, dtype=torch.float32))) > 0.5
        pred_fp = pred_fp.detach().cpu().numpy()
    
    pred_fp_list = []
    for idx in range(pred_fp.shape[0]):
        pred_fp_list.append(array2bits(pred_fp[idx]))
    df_output = pd.DataFrame({
        'smiles': smiles_list,
        'fingerprint': pred_fp_list
    })
    df_output.to_csv(save_path, index=False)

def pred_real_sample(sample_dir, checkpoint, save_dir):
    model = load_model(checkpoint=checkpoint)
    for f in os.listdir(sample_dir):
        print(f)
        with open(os.path.join(sample_dir, f, 'normalized_spec.msp'), 'r') as msp:
            for _ in range(4):
                msp.readline()
            df = pd.read_csv(sio(msp.read()), sep=' ', header=None)
        x = np.zeros((2000,), dtype=np.float32)
        x[df[0]] = df[1]
        x /= np.max(x) + 10**-6

        with torch.no_grad():
            pred_fp = torch.sigmoid(model(torch.tensor(x, device=device, dtype=torch.float32).view(1, -1))) > 0.5
            pred_fp = pred_fp.detach().cpu().numpy()

        pred_fp_str = array2bits(pred_fp[0])
        with open(os.path.join(sample_dir, f, 'normalized_spec.msp'), 'r') as msp:
            origin_lines = msp.readlines()
        with open(os.path.join(save_dir, f + '.msp'), 'w') as output_file:
            output_file.writelines(origin_lines[:3] + ['FP: ' + pred_fp_str + '\n'] + origin_lines[3:])

if __name__ == '__main__':
    pass