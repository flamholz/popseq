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
				tile_size, step_size, min_score, min_match,
				max_gap, one_off, rep_match, output_type='pslx'):
	"""Returns a shell to BLAT search with."""
	command = ['blat', db_fname, query_fname, output_fname,
			   '-q=dna', '-t=dna',  # DB and query are DNA
			   '-tileSize=%d' % tile_size,
			   '-stepSize=%d' % step_size,
			   '-minScore=%d' % min_score,
			   '-minMatch=%d' % min_match,
			   '-maxGap=%d' % max_gap,
			   '-oneOff=%d' % one_off,
			   '-repMatch=%d' % rep_match,
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
