#!/usr/bin/python

"""Filters reads.

Aligns reads to reference sequences (e.g. insert) and then produces a FASTA file
with only those reads that match.
"""

import argparse
import numpy as np
import os
import subprocess
import time

from os import path
from Bio import SearchIO


def MakeFASTAFilename(fastq_path, dest_dir):
	head, tail = path.split(fastq_path)
	name, ext = path.splitext(tail)
	return path.join(dest_dir, '%s.fa' % name)

def MakePSLFilename(fasta_path, dest_dir):
	head, tail = path.split(fasta_path)
	name, ext = path.splitext(tail)
	return path.join(dest_dir, '%s_aligned.psl' % name)

def CheckInstalled(name):
	test_command = ['which', name]
	ret = subprocess.check_output(test_command)
	assert ret, '%s is not installed!' % name

def CheckAllExist(files):
	for fname in files:
		assert path.exists(fname), '%s does not exist!' % fname

def Main():
	parser = argparse.ArgumentParser(description='Filter reads.')
	parser.add_argument("database_filename",
      					help="Path to FASTA file containing DB to align reads to.")
	parser.add_argument("read_filenames", nargs='+', help="Path to FASTQ file containing reads.")
	parser.add_argument("-t", "--tmp_dir",
					    help="Path to use to store intermediate files.",
						default="_tmp_data")
	args = parser.parse_args()

	# Check that everything we need exists.
	map(CheckInstalled, ['fastq_to_fasta', 'blat'])
	CheckAllExist(args.read_filenames)	

	# Make the temporary directory if needed.
	if not path.exists(args.tmp_dir):
		os.makedirs(args.tmp_dir)

	print 'Converting FASTQ to FASTA'
	start_ts = time.time()

	fasta_fnames = []
	for fname in args.read_filenames:
		out_fname = MakeFASTAFilename(fname, args.tmp_dir)
		fasta_fnames.append(out_fname)
		if path.exists(out_fname):
			print '\tSkipping conversion of %s as output exists' % fname
			continue

		to_fasta_command = ['fastq_to_fasta',
						    '-i', fname,
						    '-o', out_fname,
							'-Q33']
		ret = subprocess.call(to_fasta_command)
		assert ret == 0, '%s failed' % to_fasta_command[0]
	
	duration = time.time() - start_ts
	print 'Finished converting to FASTA, took %.3f seconds' % duration

	print 'Aligning reads to database at %s' % args.database_filename
	start_ts = time.time()
	psl_fnames = []
	for fasta_fname in fasta_fnames:
		psl_fname = MakePSLFilename(fasta_fname, args.tmp_dir)
		psl_fnames.append(psl_fname)
		if path.exists(psl_fname):
			print '\tSkipping alignment of %s as output exists' % fasta_fname
			# Assume reads haven't changed...
			continue

		blat_command = ['blat',
						'-tileSize=10',
						'-stepSize=4',
						'-minScore=12',
						'-q=dna', '-t=dna',
						args.database_filename,
						fasta_fname,
						psl_fname]
		
		ret = subprocess.call(blat_command)
		assert ret == 0, 'blat failed for %s!' % fasta_fname
	duration = time.time() - start_ts
	print 'Finished BLAT alignment to database, took %.3f seconds' % duration

	for fasta_fname, psl_fname in zip(fasta_fnames, psl_fnames):
		parsed = SearchIO.parse(psl_fname, 'blat-psl')
		ids_with_hits = set()
		for record in parsed:
			for hsp in record.hsps:
				ids_with_hits.add(hsp.query_id)

		print len(ids_with_hits), 'unique sequences with hits.'
	print 'Done'


if __name__ == '__main__':
	Main()
	

