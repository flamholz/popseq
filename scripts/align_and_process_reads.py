#!/usr/bin/python

"""Aligns reads to insert and backbone, finds matches.

One step script -- calls BLAT to perform alignment. Collects data from BLAT
to find those reads which match both insert and backbone. Uses match
statistics to calculate the position of the match and various other
data about it. Data is output to a CSV file which can be analyzed further.

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
import csv
import os
import time
import numpy as np

from os import path
from scripts.util import filename_util
from scripts.util import command_util
from scripts.util.fasta_stuff import ConvertFASTQToFASTA, AlignReadsToDB
from sequtils.read_alignment_data import ReadAlignmentData


def Main():
    parser = argparse.ArgumentParser(description='Filter reads.')
    parser.add_argument("-i", "--insert_db_filename", required=True,
                        help=("Path to FASTA file containing insert ends to align reads to. "
                              "Will only retain reads that align well to this DB."))
    parser.add_argument("-b", "--backbone_db_filename", required=True,
                        help=("Path to FASTA file containing backbone sequence. "
                              "Will bin reads by where they align to this sequence."))
    parser.add_argument("-r", "--read_filenames", nargs='+', required=True,
                        help="Path to FASTQ files containing reads.")
    parser.add_argument("-t", "--tmp_dir",
                        default="_read_filter_data",
                        help="Path to use to store intermediate files.")
    parser.add_argument("--fastq_ignore_degenerate", default=False, action='store_true',
                        help="If set, ignore reads with degenerate bases (N).")
    parser.add_argument("--blat_tile_size", type=int, default=6,
                        help="Tile size to use for BLAT search.")
    parser.add_argument("--blat_step_size", type=int, default=2,
                        help="Step size to use for BLAT search.")
    parser.add_argument("--blat_min_score", type=int, default=15,
                        help="Minimum score to retain a BLAT match.")
    parser.add_argument("--blat_max_gap", type=int, default=0,
                        help="Blat maximum number of gaps between tiles.")
    parser.add_argument("--blat_output_type", default="pslx",
                        help="Blat output format")
    parser.add_argument("-o", "--summary_output_csv_filename",
                        help="Where to write CSV output to.")
    parser.add_argument("-n", "--no_summary_output", dest='summary_output', action='store_false')
    parser.add_argument("-s", "--summary_output", dest='summary_output', action='store_true')
    parser.set_defaults(summary_output=True)
    args = parser.parse_args()

    # Check that everything we need exists.
    command_util.CheckAllInstalled(['fastq_to_fasta', 'blat'])

    # Get the filenames we are supposed to process.
    read_filenames = filename_util.ForceExpand(args.read_filenames)
    print 'Read filenames', read_filenames
    read_filenames = filter(lambda n: path.splitext(n)[1] in ['.fq', '.fastq', '.fa', '.fasta'],
                            read_filenames)
    print 'Read filenames:', ','.join(read_filenames)
    assert len(read_filenames) > 0, 'Must provide reads!'
    if args.summary_output:
        assert args.summary_output_csv_filename, 'Must provide output filename to write output.'
    
    # Check that all the input files exist.
    filename_util.CheckAllExist(read_filenames)

    # Make the temporary directory if needed.
    if not path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    # Convert the FASTQ input files to FASTA as BLAT doesn't seem
    # to take FASTQ input. NOTE: seems to ditch sequences containing uncalled bases
    # reported as "N" in this conversion. Revisit in future.
    all_fa = np.all(map(lambda fname:fname.endswith('.fa'), read_filenames))
    start_ts = time.time()
    if not all_fa:
        print 'Converting FASTQ to FASTA'
        keep_degenerate = not args.fastq_ignore_degenerate
        fasta_fnames = ConvertFASTQToFASTA(read_filenames, args.tmp_dir,
                                           keep_degenerate=keep_degenerate)
        duration = time.time() - start_ts
        print 'Finished converting to FASTA, took %.3f seconds' % duration
    else:
        fasta_fnames = read_filenames

    # Align the reads in FASTA format to the insert.
    print 'Aligning reads to insert database at %s' % args.insert_db_filename
    start_align_ts = time.time()
    insert_psl_fnames = AlignReadsToDB(
        args.insert_db_filename, fasta_fnames, args.tmp_dir,
        output_filename_postfix='insert_aligned',
        blat_tile_size=args.blat_tile_size,
        blat_step_size=args.blat_step_size,
        blat_min_score=args.blat_min_score,
        blat_max_gap=args.blat_max_gap,
        output_type='pslx')
    align_duration = time.time() - start_align_ts
    print 'Finished BLAT alignment to database, took %.3f seconds' % align_duration

    # Align the reads to the backbone.
    print 'Aligning reads to backbone database at %s' % args.backbone_db_filename
    start_align_ts = time.time()
    backbone_psl_fnames = AlignReadsToDB(
        args.backbone_db_filename, fasta_fnames, args.tmp_dir,
        output_filename_postfix='backbone_aligned',
        blat_tile_size=args.blat_tile_size,
        blat_step_size=args.blat_step_size,
        blat_min_score=args.blat_min_score,
        blat_max_gap=args.blat_max_gap,
        output_type='pslx')
    align_duration = time.time() - start_align_ts
    print 'Finished BLAT alignment to backbone, took %.3f seconds' % align_duration

    total_duration = time.time() - start_ts
    print 'Done aligning, took %.2f minutes' % (total_duration / 60.0)

    if not args.summary_output:
        print 'Not asked for summary output, bailing'
        return

    insert_aligned_fnames = insert_psl_fnames
    backbone_aligned_fnames = backbone_psl_fnames
    fasta_fnames = fasta_fnames

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
    start_ts = time.time()
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
    
    insertions = [r.has_insertion for r in read_data_by_id.itervalues()]
    fwd_insertions = [r.has_forward_insertion for r in read_data_by_id.itervalues()]
    
    n_total_w_matches = len(read_data_by_id) 
    n_insertions = np.sum(insertions)
    n_fwd_insertions = np.sum(fwd_insertions)
    
    print 'Total reads with any matches:', n_total_w_matches
    print 'Reads with insertions', n_insertions
    print 'Reads with forward insertions', n_fwd_insertions
    
    total_duration = time.time() - start_ts
    print 'Done collecting read statistics, took %.2f minutes' % (total_duration / 60.0)
    
    start_ts = time.time()
    out_fname = args.summary_output_csv_filename
    print 'Writing insertion matches to', out_fname
    with open(out_fname, 'w') as f:
        w = csv.DictWriter(f, ReadAlignmentData.DICT_FIELDNAMES)
        w.writeheader()
        for rd in read_data_by_id.itervalues():
            # Requires both matches they may not be forward or consistent.
            if rd.has_insert_backbone_matches:
                w.writerow(rd.AsDict())
    total_duration = time.time() - start_ts
    print 'Done writing read statistics, took %.2f minutes' % (total_duration / 60.0)

if __name__ == '__main__':
    Main()
