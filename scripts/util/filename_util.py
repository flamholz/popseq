#!/usr/bin/python

import glob
import itertools

from os import path


def ForceExpand(paths):
	"""Forcibly expand a list of (possibly-globbed) filenames."""
	filenames = map(glob.glob, paths)
	filenames = list(itertools.chain(*filenames))
	return filenames


def CheckAllExist(paths):
	"""Assertion fails if any of the given paths does not exist."""
	for fname in paths:
		assert path.exists(fname), '%s does not exist!' % fname


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


