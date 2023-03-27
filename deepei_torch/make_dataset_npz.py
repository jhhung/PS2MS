from scipy.sparse import load_npz, save_npz, csr_matrix, vstack
import os
import pandas as pd
import pickle
import numpy as np
import logging
import argparse
from parse_spectrum import parse_sdf

def get_args() -> argparse.Namespace:
	main_parser = argparse.ArgumentParser('make_dataset_npz')
	sub_manager = main_parser.add_subparsers(help="subcommands", dest='subcommands')

	make_parser = sub_manager.add_parser('make')
	make_parser.add_argument('-f', '--fingerprint', required=True)
	make_parser.add_argument('-s', '--spectrum', required=True)
	make_parser.add_argument('-p', '--prefix', required=True)

	concat_parser = sub_manager.add_parser('concat')
	concat_parser.add_argument("-f", "--first", help="prefix of first dataset", required=True)
	concat_parser.add_argument("-s", "--second", help="prefix of second dataset", required=True)
	concat_parser.add_argument('-p', '--prefix', required=True)

	args = main_parser.parse_args()
	print(args)
	return args

def pd2csr_matrix(pkl) -> tuple[csr_matrix, pd.Series]:
	logger = logging.getLogger('pd2csr_matrix')
	logger.setLevel(logging.DEBUG)

	logger.info("load pkl...")
	df = pd.read_pickle(pkl)
	fp_series = df['fingerprint']

	logger.info("convert to np array...")
	fp_numpy = np.array(list(map(np.array, fp_series)))

	logger.info("convert to csr matrix...")
	fp_sparse = csr_matrix(fp_numpy)

	logger.info(f"matrix shape: {fp_sparse.shape}")
	assert fp_sparse.shape == fp_numpy.shape

	return fp_sparse, df['smiles']



if __name__ == "__main__":
	logging.basicConfig(format='[%(asctime)s] [%(name)s] <%(levelname)s>: %(message)s')
	logger = logging.getLogger('main')
	logger.setLevel(logging.DEBUG)

	args = get_args()

	if args.subcommands == "make":
		if not os.path.exists(os.path.dirname(args.prefix)):
			raise FileExistsError(f"Dir '{os.path.dirname(args.prefix)}' not exist, can't create output")

		logger.info("Making Fingerprint dataset...")
		fp_sparse, fp_smiles = pd2csr_matrix(args.fingerprint)

		logger.info("Making Spectrum dataset...")
		logger.info("parse sdf...")
		spec_numpy, spec_smiles = parse_sdf(args.spectrum)

		# check matched (for the shared part) like:
		##s1 s2
		#------
		# AA AA
		# BB CC
		#    DD
		# this case DD will not be taken into consideration
		for idx, (fp_s, spec_s) in enumerate(zip(fp_smiles, spec_smiles)):
			if fp_s != spec_s:
				logger.error(f"fp_smiles: {fp_s}, spec_smiles: {spec_s}, idx: {idx} !!")
				exit(1)
		

		logger.info("make csr matrix...")
		spec_sparse = csr_matrix(spec_numpy)
		assert spec_sparse.shape == spec_numpy.shape
		max_common = min(fp_sparse.shape[0], spec_sparse.shape[0])
		
		# take the shared part
		if max_common != fp_sparse.shape[0] or max_common != spec_sparse.shape[0]:
			fp_sparse = fp_sparse[:max_common, :]
			spec_sparse = spec_sparse[:max_common, :]
		logger.info(f"matrix shape: fp: {fp_sparse.shape}, spec: {spec_sparse.shape}")

		logger.info(f"Saving to {args.prefix} ...")
		save_npz(f"{args.prefix}_fingerprint.npz", fp_sparse)
		save_npz(f"{args.prefix}_spec.npz", spec_sparse)
		logger.info("Done.")

	elif args.subcommands == "concat":
		if not os.path.exists(args.first + "_fingerprint.npz") or not os.path.exists(args.second + "_fingerprint.npz") or \
			not os.path.exists(args.first + "_spec.npz") or not os.path.exists(args.second + "_spec.npz"):
			raise FileNotFoundError(f"Dataset not found!")
		if not os.path.exists(os.path.dirname(args.prefix)):
			raise FileExistsError(f"Dir '{os.path.dirname(args.prefix)}' not exist, can't create output")
		
		first_fp = load_npz(args.first + "_fingerprint.npz")
		first_spec = load_npz(args.first + "_spec.npz")
		second_fp = load_npz(args.second + "_fingerprint.npz")
		second_spec = load_npz(args.second + "_spec.npz")

		for npz_matrix in [first_fp, first_spec, second_fp, second_spec]:
			logger.info(f"Shape: {npz_matrix.shape}")
		
		
		logger.info("Stacking fp...")
		concated_fp = vstack([first_fp, second_fp])
		logger.info("Stacking spec...")
		concated_spec = vstack([first_spec, second_spec])

		save_npz(args.prefix + "_fingerprint.npz", concated_fp)
		save_npz(args.prefix + "_spec.npz", concated_spec)
		logger.info("Done.")
	else:
		raise NotImplementedError

