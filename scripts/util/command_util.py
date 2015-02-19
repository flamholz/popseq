#!/usr/bin/python

import subprocess


def ToFASTACommand(fastq_path, fasta_path, keep_degenerate=True):
	"""Returns a shell command to convert the given FASTQ to FASTA."""
	l = ['fastq_to_fasta',
	    '-i', fastq_path,
	    '-o', fasta_path,
		'-Q33']
	if keep_degenerate:
		l.append('-n')
	return l

def BLATCommand(db_fname, query_fname, output_fname,
				blat_tile_size=10, blat_step_size=3,
				blat_min_score=10, blat_max_gap=0,
				output_type='psl'):
	"""Returns a shell to BLAT search with."""
	command = ['blat', db_fname, query_fname, output_fname,
			   '-q=dna', '-t=dna',  # DB and query are DNA
			   '-tileSize=%d' % blat_tile_size,
			   '-stepSize=%d' % blat_step_size,
			   '-minScore=%d' % blat_min_score,
               '-out=%s' % output_type] 
	return command


def CheckInstalled(name):
	"""Assertion fails if the given shell command is not available."""
	test_command = ['which', name]
	ret = subprocess.check_output(test_command)
	assert ret, '%s is not installed!' % name


def CheckAllInstalled(names):
	"""Assertion fails if any name is not installed."""
	map(CheckInstalled, names)
