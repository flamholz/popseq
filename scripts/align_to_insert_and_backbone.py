#!/usr/bin/python

"""Filters reads.

Aligns reads to reference sequences (e.g. insert) and then produces a FASTA file
with only those reads that match.

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
import os
import time

from os import path
from scripts.util import filename_util
from scripts.util import command_util
from scripts.util.fasta_stuff import ConvertFASTQToFASTA, AlignReadsToDB


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
    parser.add_argument("--blat_tile_size", type=int, default=10,
                        help="Tile size to use for BLAT search.")
    parser.add_argument("--blat_step_size", type=int, default=3,
                        help="Step size to use for BLAT search.")
    parser.add_argument("--blat_min_score", type=int, default=10,
                        help="Minimum score to retain a BLAT match.")
    parser.add_argument("--blat_output_type", default="pslx",
                        help="Blat output format")
    args = parser.parse_args()

    # Check that everything we need exists.
    command_util.CheckAllInstalled(['fastq_to_fasta', 'blat'])

    # Get the filenames we are supposed to process.
    read_filenames = filename_util.ForceExpand(args.read_filenames)
    print 'Read filenames', read_filenames
    read_filenames = filter(lambda n: n.endswith('fq') or n.endswith('fastq'),
                            read_filenames)
    print 'Read filenames:', ','.join(read_filenames)
    assert len(read_filenames) > 0, 'Must provide reads!'
    # Check that all the input files exist.
    filename_util.CheckAllExist(read_filenames)

    # Make the temporary directory if needed.
    if not path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    # Convert the FASTQ input files to FASTA as BLAT doesn't seem
    # to take FASTQ input. NOTE: seems to ditch sequences containing uncalled bases
    # reported as "N" in this conversion. Revisit in future.
    print 'Converting FASTQ to FASTA'
    start_ts = time.time()
    fasta_fnames = ConvertFASTQToFASTA(read_filenames, args.tmp_dir)
    duration = time.time() - start_ts
    print 'Finished converting to FASTA, took %.3f seconds' % duration

    # Align the reads in FASTA format to the insert.
    print 'Aligning reads to insert database at %s' % args.insert_db_filename
    start_align_ts = time.time()
    insert_psl_fnames = AlignReadsToDB(
        args.insert_db_filename, fasta_fnames, args.tmp_dir,
        output_filename_postfix='insert_aligned',
        blat_tile_size=args.blat_tile_size,
        blat_step_size=args.blat_step_size,
        blat_min_score=args.blat_min_score,
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
        output_type='pslx')
    align_duration = time.time() - start_align_ts
    print 'Finished BLAT alignment to backbone, took %.3f seconds' % align_duration

    total_duration = time.time() - start_ts
    print 'Done, full script took %.2f minutes' % (total_duration / 60.0)


if __name__ == '__main__':
    Main()
