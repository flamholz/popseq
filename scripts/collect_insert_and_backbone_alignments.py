# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

import glob
import numpy as np

from sequtils.read_alignment_data import ReadAlignmentData

insert_aligned_fnames = glob.glob('_read_filter_data/DFS001_3_index7_CAGATC_L004_R1_*_insert_aligned.pslx')
backbone_aligned_fnames = glob.glob('_read_filter_data/DFS001_3_index7_CAGATC_L004_R1_*_backbone_aligned.pslx')
fasta_fnames = glob.glob('_read_filter_data/DFS001_3_index7_CAGATC_L004_R1_*.fa')

"""
insert_aligned_fnames = glob.glob('sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped_insert_aligned.pslx')
backbone_aligned_fnames = glob.glob('sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped_backbone_aligned.pslx')
fasta_fnames = glob.glob('sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped.fa')
"""

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
distances = np.array([r.BackboneInsertDistanceInRead() for r in read_data_by_id.itervalues()])

n_total_w_matches = len(read_data_by_id) 
n_consistent = np.sum(consistent)

print 'Total reads with any matches:', n_total_w_matches
print 'Reads with consistent matches to backbone and insert', n_consistent

positions = []
positions_3p = []
positions_5p = []
for r in read_data_by_id.itervalues():
    # NOTE: reverse insertions are when the insert and backbone matches are on opposite
    # strands of the read. We may want to quantify these later, but we are ignoring them
    # for now.
    if r.ConsistentInsertAndBackboneMatches():
        ip = r.insertion_site
        if ip < 0:
            print 'negative insertion site'
            print r.insert_hsp
            print r.backbone_hsp
        positions.append(ip)
        if r._insert_match_end == '3p':
            positions_3p.append(ip)
        else:
            positions_5p.append(ip)

import csv

output_fname = 'all_insertions_sample_3.csv'
print 'Writing insertion matches to', output_fname
with open(output_fname, 'w') as f:
    w = csv.DictWriter(f, ReadAlignmentData.DICT_FIELDNAMES)
    w.writeheader()
    for rd in read_data_by_id.itervalues():
        if rd.has_insertion:
            w.writerow(rd.AsDict())

import pylab
pylab.figure()
pylab.hist(positions_3p, bins=50, color='b', label='3p insertions')
pylab.hist(positions_5p, bins=50, color='g', label='5p insertions')
pylab.xlabel('Insertion Site (nt)')
pylab.xlim((-30, 4140))
pylab.legend()
pylab.savefig('3p_5p_insert_dist_sample_3.png')
pylab.savefig('3p_5p_insert_dist_sample_3.svg')

pylab.figure()
pylab.hist(positions, bins=50, color='b')
pylab.xlabel('Insertion Site (nt)')
pylab.xlim((-30, 4140))
pylab.savefig('insert_dist_sample_3.png')
pylab.savefig('insert_dist_sample_3.svg')

pylab.figure()
finite_idx = np.isfinite(distances)
pylab.hist(distances[finite_idx], bins=50)
pylab.xlabel('Insert - Backbone Distance (nt)')
pylab.savefig('insert_backbone_distance_sample_3.png')
pylab.savefig('insert_backbone_distance_sample_3.svg')
pylab.show()

