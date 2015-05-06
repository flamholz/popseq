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


def MakeFname(in_fname, out_ext,
			  dest_dir=None,
			  postfix=None):
	"""Helper to make a new filename from an existing one.

	Args:
		in_fname: the existing filename.
		out_ext: the extension of the ouptut filename.
		dest_dir: directory to put the output filename into.
		postfix: a string to append to the filename.

	Returns:
		New filename for output.
	"""
	head, tail = path.split(in_fname)
	name, ext = path.splitext(tail)
	fname = '%s.%s' % (name, out_ext)
	if postfix:
		fname = '%s_%s.%s' % (name, postfix, out_ext)
	if dest_dir:
		return path.join(dest_dir, fname)
	return fname


def MakeFASTAFilename(in_fname,
					  dest_dir=None,
					  postfix=None):
	"""Makes a FASTA filename for the given input file.

	Args:
		in_fname: path to the existing file.
		dest_dir: where you want to save the FASTA file.
		postfix: a string to append to the filename (before extension).
	
	Returns:
		The FASTA filename.
	"""
	return MakeFname(in_fname, out_ext='fa',
					 dest_dir=dest_dir,
					 postfix=postfix)


def MakePSLFilename(fasta_path,
					dest_dir=None,
					postfix=None,
                    out_ext='psl'):
	"""Makes a BLAT PSL filename for the given FASTA file.

	Args:
		fasta_path: path to the FASTA file.
		dest_dir: where you want to save the PSL.
	
	Returns:
		The PSL filename.
	"""
	return MakeFname(fasta_path, out_ext=out_ext,
					 dest_dir=dest_dir,
					 postfix=postfix)

