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


def MakeFASTAFilename(fastq_path, dest_dir):
	"""Makes a FASTA filename for the given FASTQ file.

	Args:
		fastq_path: path to the FASTQ file.
		dest_dir: where you want to save the FASTA file.
	
	Returns:
		The FASTA filename.
	"""
	head, tail = path.split(fastq_path)
	name, ext = path.splitext(tail)
	return path.join(dest_dir, '%s.fa' % name)


def MakeFilteredFASTAFilename(fasta_path):
	"""Makes a filename for a filtered version of the given file.

	Args:
		fasta_path: path to the FASTA file.
	
	Returns:
		The FASTA filename.
	"""
	name, ext = path.splitext(fasta_path)
	return '%s_filtered.fa' % name


def MakePSLFilename(fasta_path, dest_dir):
	"""Makes a BLAT PSL filename for the given FASTA file.

	Args:
		fasta_path: path to the FASTA file.
		dest_dir: where you want to save the PSL.
	
	Returns:
		The PSL filename.
	"""
	head, tail = path.split(fasta_path)
	name, ext = path.splitext(tail)
	return path.join(dest_dir, '%s_aligned.psl' % name)


def ToFASTACommand(fastq_path, fasta_path):
	"""Returns a shell command to convert the given FASTQ to FASTA."""
	return ['fastq_to_fasta',
		    '-i', fastq_path,
		    '-o', fasta_path,
			'-Q33']


def BLATCommand(db_fname, query_fname, output_fname,
				**kwargs):
	"""Returns a shell to BLAT search with."""
	command = ['blat', db_fname, query_fname, output_fname,
			   '-q=dna', '-t=dna']  # DB and Query are DNA.
	if 'blat_tile_size' in kwargs:
		command.append('-tileSize=%d' % kwargs['blat_tile_size'])
	if 'blat_step_size' in kwargs:
		command.append('-stepSize=%d' % kwargs['blat_step_size'])
	if 'blat_min_score' in kwargs:
		command.append('-minScore=%d' % kwargs['blat_min_score'])
	return command


def CheckInstalled(name):
	"""Assertion fails if the given shell command is not available."""
	test_command = ['which', name]
	ret = subprocess.check_output(test_command)
	assert ret, '%s is not installed!' % name


def CheckAllExist(paths):
	"""Assertion fails if any of the given paths does not exist."""
	for fname in paths:
		assert path.exists(fname), '%s does not exist!' % fname


def Main():
	parser = argparse.ArgumentParser(description='Filter reads.')
	parser.add_argument("database_filename",
      					help="Path to FASTA file containing DB to align reads to.")
	parser.add_argument("-r", "--read_filenames", nargs='+', help="Path to FASTQ file containing reads.")
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
	map(CheckInstalled, ['fastq_to_fasta', 'blat'])

	read_filenames = map(glob.glob, args.read_filenames)
	read_filenames = list(itertools.chain(*read_filenames))
	read_filenames = filter(lambda n: n.endswith('fq') or n.endswith('fastq'),
							read_filenames)
	print 'Read filenames:', read_filenames
	assert len(read_filenames) > 0, 'Must provide reads!'
	CheckAllExist(read_filenames)

	# Make the temporary directory if needed.
	if not path.exists(args.tmp_dir):
		os.makedirs(args.tmp_dir)

	print 'Converting FASTQ to FASTA'
	start_ts = time.time()

	fasta_fnames = []
	for fname in read_filenames:
		out_fname = MakeFASTAFilename(fname, args.tmp_dir)
		fasta_fnames.append(out_fname)
		if path.exists(out_fname):
			print '\tSkipping conversion of %s as output exists' % fname
			continue

		to_fasta_command = ToFASTACommand(fname, out_fname)
		ret = subprocess.call(to_fasta_command)
		assert ret == 0, '%s failed' % to_fasta_command[0]
	
	duration = time.time() - start_ts
	print 'Finished converting to FASTA, took %.3f seconds' % duration

	print 'Aligning reads to database at %s' % args.database_filename
	start_align_ts = time.time()
	psl_fnames = []
	blat_kwargs = {'blat_step_size': args.blat_step_size,
				   'blat_tile_size': args.blat_tile_size,
				   'blat_min_score': args.blat_min_score}

	for fasta_fname in fasta_fnames:
		psl_fname = MakePSLFilename(fasta_fname, args.tmp_dir)
		psl_fnames.append(psl_fname)
		if path.exists(psl_fname):
			print '\tSkipping alignment of %s as output exists' % fasta_fname
			continue

		blat_command = BLATCommand(
			args.database_filename, fasta_fname, psl_fname,
			**blat_kwargs)
		ret = subprocess.call(blat_command)
		assert ret == 0, 'blat failed for %s!' % fasta_fname
	align_duration = time.time() - start_align_ts
	print 'Finished BLAT alignment to database, took %.3f seconds' % duration

	start_filter_ts = time.time()
	for fasta_fname, psl_fname in zip(fasta_fnames, psl_fnames):
		filtered_fname = MakeFilteredFASTAFilename(fasta_fname)
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

		# Force the GCs hand since these may be very big.
		del retained
		del ids_with_hits

	duration = time.time() - start_filter_ts
	print 'Finished filtering output, took %.3f seconds' % duration

	total_duration = time.time() - start_ts
	print 'Done, full script took %.2f minutes' % (total_duration / 60.0)


if __name__ == '__main__':
	Main()
