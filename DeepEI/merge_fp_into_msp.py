import pandas as pd
from scipy.sparse import load_npz
from read_msp_lib import *
import sys

msp_file = open(sys.argv[1])
fp_table = pd.read_pickle(sys.argv[2]) # smiles, fingerprint, mw
outfile = open(sys.argv[3], 'w')

i = 0
spectrum = read_a_spectrum(msp_file)
while len(spectrum) != 0:
    fp = fp_table[fp_table['smiles'] == spectrum['NAME']]['fingerprint'].values
    if len(fp) != 0:
        spectrum["FP"] = "".join([str(item) for item in fp[0]])
    write_a_spectrum(outfile, spectrum)
    spectrum = read_a_spectrum(msp_file)
    if i % 10000 == 0:
        print("Num: {0}".format(i), end='\r')
    i += 1
print()

