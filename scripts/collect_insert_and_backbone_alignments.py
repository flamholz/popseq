# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

import argparse
import numpy as np

from scripts.util import filename_util
from sequtils.read_alignment_data import ReadAlignmentData


parser = argparse.ArgumentParser(description='Filter reads.')
parser.add_argument("-i", "--insert_alignments", nargs='+', required=True,
                    help=("Path to pslx files containing alignments to backbone"))
parser.add_argument("-b", "--backbone_alignments", nargs='+', required=True,
                    help=("Path to pslx files containing alignments to insert"))
parser.add_argument("-r", "--read_filenames", nargs='+', required=True,
                    help="Path to FASTA files containing reads.")
parser.add_argument("-o", "--output_csv_filename", required=True,
                    help="Where to write CSV output to.")
args = parser.parse_args()

insert_aligned_fnames = args.insert_alignments
backbone_aligned_fnames = args.backbone_alignments
fasta_fnames = args.read_filenames

insert_aligned_fnames = filename_util.ForceExpand(insert_aligned_fnames)
backbone_aligned_fnames = filename_util.ForceExpand(backbone_aligned_fnames)
fasta_fnames = filename_util.ForceExpand(fasta_fnames)

"""
# For testing the script with less data
import glob
insert_aligned_fnames = glob.glob('sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_insert_aligned.pslx')
backbone_aligned_fnames = glob.glob('sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_backbone_aligned.pslx')
fasta_fnames = glob.glob('sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k.fa')
"""

# Make sure all the filenames are in the same order so we can zip them.
insert_aligned_fnames = sorted(insert_aligned_fnames)
backbone_aligned_fnames = sorted(backbone_aligned_fnames)
fasta_fnames = sorted(fasta_fnames)
assert insert_aligned_fnames
assert backbone_aligned_fnames
assert fasta_fnames
assert len(insert_aligned_fnames) == len(backbone_aligned_fnames)
assert len(fasta_fnames) == len(insert_aligned_fnames)
print 'Insert aligned filenames'
print insert_aligned_fnames
print 'Backbone aligned filenames'
print backbone_aligned_fnames
print 'FASTA fnames'
print fasta_fnames

# Gather all the reads information by read ID.
read_data_by_id = {}
for insert_fname, backbone_fname, fasta_fname in zip(insert_aligned_fnames,
                                                     backbone_aligned_fnames,
                                                     fasta_fnames):
    print 'Analyzing file set'
    print insert_fname
    print backbone_fname
    print fasta_fname
    read_data_by_id.update(ReadAlignmentData.DictFromFiles(insert_fname,
                                                           backbone_fname,
                                                           fasta_fname))

consistent = [r.has_insertion for r in read_data_by_id.itervalues()]

n_total_w_matches = len(read_data_by_id) 
n_consistent = np.sum(consistent)

print 'Total reads with any matches:', n_total_w_matches
print 'Reads with consistent matches to backbone and insert', n_consistent

import csv

print 'Writing insertion matches to', args.output_csv_filename
with open(args.output_csv_filename, 'w') as f:
    w = csv.DictWriter(f, ReadAlignmentData.DICT_FIELDNAMES)
    w.writeheader()
    for rd in read_data_by_id.itervalues():
        # Requires both matches they may not be forward or consistent.
        if rd.has_insert_backbone_matches:
            w.writerow(rd.AsDict())

