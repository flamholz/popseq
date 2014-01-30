#!/usr/bin/python

"""Filter Reads.

Aligns reads to reference sequence (specifically entire plasmid sequence)
and then produces a FASTA file with only those reads that match

"""

import argparse
import itertools
import glob
import numpy as np
import os
import subprocess
import time

from os import path
from Bio import SearchIO
from Bio import SeqIO
from scripts.util import command_util
from scripts.util import filename_util

def ConvertFASTQToFASTA(fastq_fnames, output_dir):
	"""Calls out to compiled executable to convert FASTQ files to FASTA.
	
	Args:
		fastq_fnames: iterable of paths to FASTQ files.

	Returns:
		A list of paths to the FASTA files.
	"""
	fasta_fnames= []
	for fname in fastq_fnames:
		out_fname = filename_util.MakeFASTAFilename(fname, output_dir)
		fasta_fnames.append(out_fname)
		if path.exists(out_fname):
			print '\tSkipping conversion of %s as output exists' % fname
			continue

		print 'Writing FASTA output to', out_fname
		to_fasta_command = command_util.ToFASTACommand(fname, out_fname)
		ret = subprocess.call(to_fasta_command)
		assert ret == 0, '%s failed' % to_fasta_command[0]
	return fasta_fnames

def AlignReadsToDB(db_fname, fasta_fnames, output_dir,
				   output_filename_postfix=None,
				   blat_tile_size=10,
				   blat_step_size=3,
				   blat_min_score=10):
	"""Calls out to BLAT to align reads to the given sequences.

	Args:
		db_fname: the path to the file containing the target sequences to align to.
		fasta_fnames: the read-containing filenames to use as queries in alignment.
	
	Returns:
		A list of PSL filenames with the alignment output.
	"""
	psl_fnames = []
	for fasta_fname in fasta_fnames:
		psl_fname = filename_util.MakePSLFilename(
			fasta_fname, output_dir, postfix=output_filename_postfix)
		psl_fnames.append(psl_fname)
		if path.exists(psl_fname):
			print '\tSkipping alignment of %s as output exists' % fasta_fname
			continue

		print '\tWriting output to', psl_fname
		blat_command = command_util.BLATCommand(
			db_fname, fasta_fname, psl_fname,
			blat_tile_size=blat_tile_size,
			blat_step_size=blat_step_size,
			blat_min_score=blat_min_score)
		ret = subprocess.call(blat_command)
		assert ret == 0, 'blat failed for %s!' % fasta_fname
	return psl_fnames

def Main():
	parser = argparse.ArgumentParser(description='Filter reads.')
	parser.add_argument("-i", "--plasmid_db_filename", required=True,
      					help=("Path to FASTA file containing plasmid sequence."
							  "Will only retain reads that align well to this DB."))
	parser.add_argument("-r", "--read_filenames", nargs='+', required=True,
						help="Path to FASTQ files containing reads.")
	parser.add_argument("-t", "--tmp_dir",
						default="_read_filter_plasmid_data",
					    help="Path to use to store intermediate files.")
	parser.add_argument("--blat_tile_size", type=int, default=10,
						help="Tile size to use for BLAT search.")
	parser.add_argument("--blat_step_size", type=int, default=3,
						help="Step size to use for BLAT search.")
	parser.add_argument("--blat_min_score", type=int, default=10,
						help="Minimum score to retain a BLAT match.")
	args = parser.parse_args()

    # Check that everything we need exists.
	command_util.CheckAllInstalled(['fastq_to_fasta', 'blat'])

	# Get the filenames we are supposed to process.
	read_filenames = filename_util.ForceExpand(args.read_filenames)
	print read_filenames
	read_filenames = filter(lambda n: n.endswith('fq') or n.endswith('fastq'),
							read_filenames)
	print 'Read filenames:', ','.join(read_filenames)
	assert len(read_filenames) > 0, 'Must provide reads!'
	# Check that all the input files exist.
	filename_util.CheckAllExist(read_filenames)

	# Make the temporary directory if needed.
	if not path.exists(args.tmp_dir):
		os.makedirs(args.tmp_dir)

    # Convert the FASTQ input files to FASTA as BLAT doesn't seem
	# to take FASTQ input.
	print 'Converting FASTQ to FASTA'
	start_ts = time.time()
	fasta_fnames = ConvertFASTQToFASTA(read_filenames, args.tmp_dir)
	duration = time.time() - start_ts
	print 'Finished converting to FASTA, took %.3f seconds' % duration

    # Align the reads in FASTA format to the plasmid.
	print 'Aligning reads to plasmid sequence: %s' % args.plasmid_db_filename
	start_align_ts = time.time()
	insert_psl_fnames = AlignReadsToDB(
		args.plasmid_db_filename, fasta_fnames, args.tmp_dir,
		output_filename_postfix='plasmid_aligned',
		blat_tile_size=args.blat_tile_size,
		blat_step_size=args.blat_step_size,
		blat_min_score=args.blat_min_score)
	align_duration = time.time() - start_align_ts
	print 'Finished BLAT alignment to plasmid sequence, took %.3f seconds' % align_duration

	total_duration = time.time() - start_ts
	print 'Done, full script took %.2f minutes' % (total_duration / 60.0)
   

if __name__ == '__main__':
	Main()
