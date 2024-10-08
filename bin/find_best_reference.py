#!/usr/bin/env python3
import os 
import pandas as pd 
from glob import glob 
from Bio import Entrez
import argparse


def collect_refseeker_files(input_path):
	print("Collecting referenceseeker files ...")

	refseeker_files = glob(os.path.join(input_path, '*', '*_referenceseeker.tsv'))

	if len(refseeker_files) == 0:
		print("ERROR: No referenceseeker files found in input directory.")
		raise ValueError

	return refseeker_files

def find_best_reference(filepaths: list) -> str:

	print("Finding best reference ...")

	dfs = []

	for filepath in filepaths:
		refseeker_df = pd.read_csv(filepath, sep='\t')
		
		refseeker_df.insert(0, 'sample_id', os.path.basename(filepath).split("_")[0])
		refseeker_df['z_mash']    = - (refseeker_df['Mash Distance'] - refseeker_df['Mash Distance'].mean()) / refseeker_df['Mash Distance'].std()
		refseeker_df['z_ani']     = (refseeker_df['ANI'] - refseeker_df['ANI'].mean()) / refseeker_df['ANI'].std()
		refseeker_df['z_con_dna'] = (refseeker_df['Con. DNA'] - refseeker_df['Con. DNA'].mean()) / refseeker_df['Con. DNA'].std()
		refseeker_df['z_total']   = refseeker_df['z_mash'] + refseeker_df['z_ani'] + refseeker_df['z_con_dna']
		refseeker_df = refseeker_df.sort_values('z_total',ascending=False)
		dfs.append(refseeker_df)

	concat_df = pd.concat(dfs)

	best_references = concat_df.groupby('#ID')['z_total'].sum().sort_values(ascending=False)

	return best_references.index[0]

def write_reference_file(assembly_id, outfile_name):

	print("Searching for top reference assembly ...")

	success = False
	attempt = 1
	while not success and attempt <= 3:
		try: 

			with Entrez.esearch(db='nucleotide', term=assembly_id, retmax=1, idtype='acc') as handle:
				results = Entrez.read(handle)
			success = True

		except RuntimeError as e :
			print(f"Encountered RuntimeEerror with Entrez esearch on attempt #{attempt}. Retrying ...")
			attempt += 1
	
	if not success:
		raise ValueError

	if "IdList" not in results or len(results['IdList']) == 0:
		print("ERROR: No GenBank entries found") 
		raise ValueError

	if len(results['IdList']) > 1:
		print("ERROR: Multiple GenBank entries found") 
		raise ValueError

	best_reference_accno = results['IdList'][0]


	print("Fetching and writing FASTA reference file ...")

	success = False
	attempt = 1
	while not success and attempt <= 3:
		try:
			with Entrez.efetch(db='nucleotide', id=best_reference_accno, rettype='fasta', retmode='text') as infile, open(outfile_name + '.fa', 'w') as outfile:
				for line in infile.readlines():
					outfile.write(line)
			
			success = True
		except RuntimeError as e:
			print(f"Encountered RuntimeEerror with Entrez efetch on attempt #{attempt}. Retrying ...")
			attempt += 1

	if not success:
		raise ValueError

	print("Fetching and writing GenBank reference file ...")

	success = False
	attempt = 1
	while not success and attempt <= 3:
		try:
			with Entrez.efetch(db='nucleotide', id=best_reference_accno, rettype='gb', retmode='text') as infile, open(outfile_name + '.gb', 'w') as outfile:
				for line in infile.readlines():
					outfile.write(line)
			
			success = True
		except RuntimeError as e:
			print(f"Encountered RuntimeEerror with Entrez efetch on attempt #{attempt}. Retrying ...")
			attempt += 1

	if not success:
		raise ValueError

def get_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input', required=True, help='Referenceseeker results path')
	parser.add_argument('-o', '--outname', required=True, help='Prefix of output file')
	return parser.parse_args()

def main():

	args = get_args()

	refseeker_files = collect_refseeker_files(args.input)

	best_reference_id = find_best_reference(refseeker_files)

	write_reference_file(best_reference_id, args.outname)

	print("Complete.")


if __name__ == '__main__':
	main()

