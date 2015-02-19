#!/usr/bin/python

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


def ConvertFASTQToFASTA(fastq_fnames, output_dir, keep_degenerate=True):
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
		to_fasta_command = command_util.ToFASTACommand(
			fname, out_fname,
			keep_degenerate=keep_degenerate)
		ret = subprocess.call(to_fasta_command)
		assert ret == 0, '%s failed' % to_fasta_command[0]
	return fasta_fnames


def AlignReadsToDB(insert_db_fname, fasta_fnames, output_dir,
                   output_filename_postfix=None,
                   blat_tile_size=10,
                   blat_step_size=3,
                   blat_min_score=10,
                   blat_max_gap=0,
                   output_type='psl'):
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
            fasta_fname, output_dir, postfix=output_filename_postfix,
            out_ext=output_type)
        insert_psl_fnames.append(psl_fname)
        if path.exists(psl_fname):
            print '\tSkipping alignment of %s as output exists' % fasta_fname
            continue

        print '\tWriting output to', psl_fname
        blat_command = command_util.BLATCommand(
            insert_db_fname, fasta_fname, psl_fname,
            blat_tile_size=blat_tile_size,
            blat_step_size=blat_step_size,
            blat_min_score=blat_min_score,
            blat_max_gap=blat_max_gap,
            output_type=output_type)
        print 'BLAT command', ' '.join(blat_command)
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
