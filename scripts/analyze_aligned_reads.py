# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

import glob

# Filenames for the first end of the paired ends.
insert_aligned_fnames = glob.glob('../popseq/_read_filter_data/dana-*R1*_insert_aligned.psl')
backbone_aligned_fnames = glob.glob('../popseq/_read_filter_data/dana*R1*_backbone_aligned.psl')
insert_aligned_fnames = sorted(insert_aligned_fnames)
backbone_aligned_fnames = sorted(backbone_aligned_fnames)

# <codecell>

import numpy as np

from Bio import SearchIO

insert_aligned_ids = {}

for insert_aligned_fname in insert_aligned_fnames:
    insert_aligned_reader = SearchIO.parse(insert_aligned_fname, 'blat-psl')
    for record in insert_aligned_reader:
        for hsp in record.hsps:
            query_strands = np.array([f.query_strand for f in hsp.fragments])
            # NOTE: as far as I can tell, all fragments always come from the
            # same query strand...
            if (hsp.query_id in insert_aligned_ids and
                query_strands[0] != insert_aligned_ids[hsp.query_id][2]):
                print 'Already saw ID', hsp.query_id
                
            hit_end = '3p' if hsp.hit_id.endswith('_3p') else '5p'
            query_loc = hsp.query_start if hit_end == '5p' else hsp.query_end
            insert_aligned_ids[hsp.query_id] = (hit_end, query_loc, query_strands[0])

# <codecell>

backbone_aligned_ids = set()
hit_counts = {}
distances_from_insert = {}

for backbone_aligned_fname in backbone_aligned_fnames:
    backbone_aligned_reader = SearchIO.parse(backbone_aligned_fname, 'blat-psl') 
    for record in backbone_aligned_reader:
        for hsp in record.hsps:
            query_id = hsp.query_id
            query_strands = np.array([f.query_strand for f in hsp.fragments])
            query_strand = query_strands[0]
            
            insert_match_end, insert_match_loc, insert_query_strand = insert_aligned_ids[query_id]
            if (query_id not in insert_aligned_ids or
                query_strand != insert_query_strand):
                # Insert and backbone match on opposing strands.
                continue
            
            # It matters if we are considering a sequence that aligned to
            # the start or end of the insert. At the start we want to
            # use the end of the alignment and at the end we use the start...
            hit_position = hsp.hit_start if insert_match_end == '3p' else hsp.hit_end
            query_position = hsp.query_end if insert_match_end == '5p' else hsp.query_start
            distance_in_query = abs(query_position - insert_match_loc)
            distances_from_insert.setdefault(hit_position, []).append(distance_in_query)
            if distance_in_query > 20:
                # Insert and backbone are far apart in the read.
                continue
            
            backbone_aligned_ids.add(hsp.query_id)
            hit_counts[hit_position] = hit_counts.get(hit_position, 0) + 1
            

# <codecell>

n_insert_aligned = len(insert_aligned_ids)
n_backbone_aligned = len(backbone_aligned_ids)
frac = float(n_backbone_aligned) / float(n_insert_aligned)
print '%.3f%%' % (frac * 100.0)

# <codecell>

import pylab

pylab.figure()
positions = sorted(distances_from_insert.keys())
mean_displacement = [np.mean(distances_from_insert[k])
                     for k in positions]
pylab.plot(positions, mean_displacement, 'r.')
pylab.xlabel('position in MBP')
pylab.ylabel('mean displacement')

# <codecell>

pylab.figure(figsize=(11,11))
positions = np.array(sorted(hit_counts.keys()))
counts = np.array([hit_counts[k] for k in positions])
mean_count = np.mean(counts) 
stddev = np.std(counts)

total_count = np.sum(counts)
threshold = mean_count / 2.0
idx_above_thresh = np.where(counts > threshold)
my_counts = counts[idx_above_thresh]
my_positions = positions[idx_above_thresh]

pylab.yscale('log')
pylab.plot(positions, counts, 'r.')
pylab.plot(positions, mean_count * np.ones(positions.size),
           'b-', label='mean')
pylab.plot(positions, threshold * np.ones(positions.size),
           'g-', label='threshold')

for position, count in zip(my_positions, my_counts):
    pylab.text(position, count, '%s (%.1f)' % (position, position / 3.0))

pylab.xlabel('Insert position')
pylab.ylabel('Count')
pylab.legend()

# <codecell>

print 'Total count', total_count
print 'Mean count %.3f +/- %.3f ' % (mean_count, stddev)
print 'Count threshold', threshold
print
print 'Position => count:'
for position, count in zip(my_positions, my_counts):
    pct_of_total = (float(count) / float(total_count)) * 100.0
    print '\t', position, '=>', count, '(%.2f%% of total, AA %.2f)' % (pct_of_total, position / 3.0)

# <codecell>


pylab.show()
