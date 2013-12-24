#!/usr/bin/python

"""Filters reads.

Aligns reads to reference sequences (e.g. insert) and then produces a FASTA file
with only those reads that match.

Requirements:
1) FASTX toolkit installed (for fastq_to_fasta).
	http://hannonlab.cshl.edu/fastx_toolkit/download.html
2) BLAT installed (for blat).
	http://genome.ucsc.edu/FAQ/FAQblat.html
3) Recent BioPython install with SearchIO for BLAT IO (1.61 or higher).
	http://biopython.org/wiki/SearchIO
4) NumPy installed.
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
from scripts.util import filename_util
from scripts.util import command_util


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


def AlignReadsToDB(insert_db_fname, fasta_fnames, output_dir,
				   output_filename_postfix=None,
				   blat_tile_size=10,
				   blat_step_size=3,
				   blat_min_score=10):
	"""Calls out to BLAT to align reads to the given sequences.

	Args:
		insert_db_fname: the path to the file containing the sequences to align to.
		fasta_fnames: the read filenames to align.
	
	Returns:
		A list of PSL filenames with the alignment output.
	"""
	insert_psl_fnames = []
	for fasta_fname in fasta_fnames:
		psl_fname = filename_util.MakePSLFilename(
			fasta_fname, output_dir, postfix=output_filename_postfix)
		insert_psl_fnames.append(psl_fname)
		if path.exists(psl_fname):
			print '\tSkipping alignment of %s as output exists' % fasta_fname
			continue

		print '\tWriting output to', psl_fname
		blat_command = command_util.BLATCommand(
			insert_db_fname, fasta_fname, psl_fname,
			blat_tile_size=blat_tile_size,
			blat_step_size=blat_step_size,
			blat_min_score=blat_min_score)
		ret = subprocess.call(blat_command)
		assert ret == 0, 'blat failed for %s!' % fasta_fname
	return insert_psl_fnames


def FilterFASTAByPSL(fasta_fnames, psl_fnames, output_dir,
					 output_filename_postfix=None):
	"""Produces a FASTA file with only those sequences in the PSL.

	Iterates over pairs of FASTA and PSL files and produces a new
	FASTA file containing only those sequences for which there was a 
	record in the PSL file.

	Args:
		fasta_fnames: the FASTA filenames.
		psl_fnames: the PSL (BLAT output) filenames in the same order.
	
	Returns:
		A list of paths to filtered FASTA files.
	"""
	filtered_fnames = []
	for fasta_fname, psl_fname in zip(fasta_fnames, psl_fnames):
		filtered_fname = filename_util.MakeFASTAFilename(
			fasta_fname, dest_dir=output_dir,
			postfix=output_filename_postfix)
		filtered_fnames.append(filtered_fname)
		if path.exists(filtered_fname):
			print 'Skipping filtering of %s as output exists' % psl_fname

		# Get the IDs of all the matching sequences.
		parsed = SearchIO.parse(psl_fname, 'blat-psl')
		ids_with_hits = set()
		for record in parsed:
			for hsp in record.hsps:
				ids_with_hits.add(hsp.query_id)

		parsed = SeqIO.parse(fasta_fname, 'fasta')
		retained = []
		n_seqs = 0
		for record in parsed:
			n_seqs += 1
			if record.id in ids_with_hits:
				retained.append(record)

		assert len(ids_with_hits) == len(retained), 'Some sequences missing!'
		pct_retained = 100 * float(len(retained)) / float(n_seqs)

		print '\tRetained %d of %d records (%.2f%%)' % (len(retained),
													   n_seqs, pct_retained)
		print '\tWriting output to', filtered_fname
		SeqIO.write(retained, filtered_fname, 'fasta')

		# Force delete these lists since they might be very big.
		del retained
		del ids_with_hits
	return filtered_fnames


def Main():
	parser = argparse.ArgumentParser(description='Filter reads.')
	parser.add_argument("-i", "--insert_db_filename", required=True,
      					help=("Path to FASTA file containing insert ends to align reads to. "
							  "Will only retain reads that align well to this DB."))
	parser.add_argument("-b", "--backbone_db_filename", required=True,
      					help=("Path to FASTA file containing backbone sequence. "
							  "Will bin reads by where they align to this sequence."))
	parser.add_argument("-r", "--read_filenames", nargs='+', required=True,
						help="Path to FASTQ files containing reads.")
	parser.add_argument("-t", "--tmp_dir",
						default="_read_filter_data",
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

	# Align the reads in FASTA format to the insert.
	print 'Aligning reads to insert database at %s' % args.insert_db_filename
	start_align_ts = time.time()
	insert_psl_fnames = AlignReadsToDB(
		args.insert_db_filename, fasta_fnames, args.tmp_dir,
		output_filename_postfix='insert_aligned',
		blat_tile_size=args.blat_tile_size,
		blat_step_size=args.blat_step_size,
		blat_min_score=args.blat_min_score)
	align_duration = time.time() - start_align_ts
	print 'Finished BLAT alignment to database, took %.3f seconds' % align_duration

	# Produce new FASTA files containing only those reads aligning with insert.
	start_filter_ts = time.time()
	filtered_fnames = FilterFASTAByPSL(
		fasta_fnames, insert_psl_fnames, args.tmp_dir,
		output_filename_postfix='insert_filtered')
	filter_duration = time.time() - start_filter_ts
	print 'Finished filtering output, took %.3f seconds' % filter_duration

	# Align the reads to the backbone.
	print 'Aligning reads to backbone database at %s' % args.backbone_db_filename
	start_align_ts = time.time()
	backbone_psl_fnames = AlignReadsToDB(
		args.backbone_db_filename, filtered_fnames, args.tmp_dir,
		output_filename_postfix='backbone_aligned',
		blat_tile_size=args.blat_tile_size,
		blat_step_size=args.blat_step_size,
		blat_min_score=args.blat_min_score)
	align_duration = time.time() - start_align_ts
	print 'Finished BLAT alignment to backbone, took %.3f seconds' % align_duration

	total_duration = time.time() - start_ts
	print 'Done, full script took %.2f minutes' % (total_duration / 60.0)


if __name__ == '__main__':
	Main()
