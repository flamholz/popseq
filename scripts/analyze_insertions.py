#!/usr/bin/env python

import argparse
import os
import subprocess
import time
import sequtils.read_alignment_data.factory as rad_factory
import logging

from Bio.Seq import Seq
from os import path
from scripts.util.filename_util import MakeFname, ForceExpand
from scripts.util import command_util
from sequtils.transposition_params import TranspositionParams

def ReadID(read_info):
    construct_num = read_info["cn"]
    read_num = read_info["rn"]
    read_id = 'c%d:r%d' % (construct_num, read_num)
    return read_id

def BBDukRetain(ref_fname, in_fname, out_fname, k):
    if path.exists(out_fname):
        print 'Output file %s exists. Bailing!' % out_fname
        return
    
    command = ['bbduk.sh',
               'ref=%s' % ref_fname,
               'in=%s' % in_fname,
               'outm=%s' % out_fname,
               'k=%d' % k,
               'hdist=1',
               'mm=f']
    print 'Running: \n\t%s' % ' '.join(command)
    ret = subprocess.call(command)
    assert ret == 0, 'bbduk.sh failed for %s!' % in_fname


def BBDukRetainMulti(ref_fname, in_fnames, out_fnames, k):
    """Encapsulate loop."""
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        BBDukRetain(ref_fname, in_fname, out_fname, k)


def BBDukMask(ref_fname, in_fname, out_fname, k, kmask):
    if path.exists(out_fname):
        print 'Output file %s exists. Bailing!' % out_fname
        return
    command = ['bbduk.sh',
               'ref=%s' % ref_fname,
               'in=%s' % in_fname,
               'out=%s' % out_fname,
               'k=%d' % k,
               'kmask=%s' % kmask,
               'hdist=1',
               'mm=f']
    ret = subprocess.call(command)
    assert ret == 0, 'bbduk.sh failed for %s!' % in_fname


def BBDukMaskMulti(ref_fname, in_fnames, out_fnames, k, kmask):
    """Encapsulate loop."""
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        BBDukMask(ref_fname, in_fname, out_fname, k, kmask)


def GrepTrim(pattern, in_fname, out_fname,
             trim_before, trim_after, trim_match):
    if path.exists(out_fname):
        print 'Output file %s exists. Bailing!' % out_fname
        return
    
    command = ['fastq-grep',
               pattern, in_fname]
    
    if trim_before:
        command.append('-b')
    if trim_after:
        command.append('-a')
    if trim_match:
        command.append('-t')
    
    command.extend(['| fastq-grep -v ZZZZ >', out_fname])
    ret = os.system(' '.join(command))
    assert ret == 0, 'fastq-grep failed for %s!' % in_fname

    
def GrepTrimMulti(pattern, in_fnames, out_fnames,
                  trim_before=False,
                  trim_after=False,
                  trim_match=False):
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        GrepTrim(pattern, in_fname, out_fname,
                 trim_before, trim_after, trim_match)

    
def BBMapAlign(index_name, in_fname, out_fname):
    if path.exists(out_fname):
        print 'Output file %s exists. Bailing!' % out_fname
        return
    
    command = ['bbmap.sh',
               'ref=%s' % index_name,
               'in=%s' % in_fname,
               'out=%s' % out_fname,
               'touppercase=t',
               'nodisk']
    ret = os.system(' '.join(command))
    assert ret == 0, 'bbmap.sh failed for %s!' % in_fname


def BBMapAlignMulti(index_name, in_fnames, out_fnames):
    for in_fname, out_fname in zip(in_fnames, out_fnames):
        BBMapAlign(index_name, in_fname, out_fname)
        
        
def SamtoolsIndex(in_fname):
    fname, unused_old_ext = os.path.splitext(in_fname)
    bam_fname = '%s.bam' % fname
    
    if path.exists(bam_fname):
        print 'Output file %s exists. Bailing!' % bam_fname
        return bam_fname
    
    # Make a sorted BAM file.
    command = ['samtools view -bS %s' % in_fname,
               '|'
               'samtools sort -o %s' % bam_fname]
    ret = os.system(' '.join(command)) # appends .bam
    assert ret == 0, 'samtools failed for %s!' % in_fname

    
    command = ['samtools index %s' % bam_fname]
    ret = os.system(' '.join(command))
    assert ret == 0, 'samtools failed to index %s!' % bam_fname

    return bam_fname


def SamtoolsIndexMutli(in_fnames):
    out_fnames = map(SamtoolsIndex, in_fnames)
    return out_fnames

def Main():
    logging.basicConfig(filename='out.log', level=logging.DEBUG)

    parser = argparse.ArgumentParser(description='Filter reads.', fromfile_prefix_chars='@')
    parser.add_argument("-i", "--insert_db_filename", required=True,
                        help=("Path to FASTA file containing insert ends to align reads to. "
                              "Will only retain reads that align well to this DB."))
    parser.add_argument("-r", "--read_filenames", nargs='+', required=True,
                        help="Path to FASTQ files containing reads.")
    parser.add_argument("-t", "--tmp_dir",
                        default="_read_filter_data",
                        help="Path to use to store intermediate files.")
    parser.add_argument("-o", "--output_fname", required=True,
                        help="Where to write output data (CSV).")
    TranspositionParams.AddArgs(parser)
    parser.set_defaults(summary_output=True)
    args = parser.parse_args()
    
    start_ts = time.time()
    command_util.CheckAllInstalled(['fastq-grep', 'bbduk.sh', 'bowtie2', 'samtools'])

    tn_params = TranspositionParams.FromArgs(args)
    print tn_params

    insert_db_fname = args.insert_db_filename
    bbone_db_fname = args.backbone_db_filename
    read_fnames = ForceExpand(args.read_filenames)
    
    assert len(read_fnames) > 0, 'There better be read files'

    print '##### Retaining reads that match insert #####'
    insert_filtered_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                        postfix='insert_filtered')
                              for i in read_fnames]
    BBDukRetainMulti(insert_db_fname,
                     read_fnames,
                     insert_filtered_fnames, k=12)
    
    print '##### Retaining reads also matching backbone #####'
    insert_bbone_filtered_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                              postfix='insert_bbone_filtered')
                                    for i in read_fnames] 
    BBDukRetainMulti(bbone_db_fname,
                     insert_filtered_fnames,
                     insert_bbone_filtered_fnames, k=12)
    
    print '##### Masking insert in reads #####'
    filtered_masked_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                        postfix='filtered_insert_masked')
                              for i in read_fnames]
    BBDukMaskMulti(insert_db_fname,
                   insert_bbone_filtered_fnames,
                   filtered_masked_fnames, k=12, kmask='Z')
    
    print '##### Trimming fixed sequence and insert from reads #####'
    fixed_5p = Seq(args.fixed_5p)
    fixed_3p = Seq(args.fixed_3p)
    pattern_5p = '%s[ATCG]{0,11}ZZZZZ' % args.fixed_5p
    pattern_3p = 'ZZZZZ[ATCG]{0,11}%s' % args.fixed_3p
    pattern_5p_rev = 'ZZZZZ[ATCG]{0,11}%s' % fixed_5p.reverse_complement()
    pattern_3p_rev = '%s[ATCG]{0,11}ZZZZZ' % fixed_3p.reverse_complement()
   
    trimmed_5p_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                   postfix='filtered_trimmed_5p')
                         for i in read_fnames]

    GrepTrimMulti('"{}"'.format(pattern_5p), filtered_masked_fnames, trimmed_5p_fnames,
                  trim_after=True, trim_match=True)
   
    trimmed_3p_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                   postfix='filtered_trimmed_3p')
                         for i in read_fnames]

    GrepTrimMulti('"{}"'.format(pattern_3p), filtered_masked_fnames, trimmed_3p_fnames,
                  trim_before=True, trim_match=True)
    
    trimmed_5p_rev_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                       postfix='filtered_trimmed_5p_rev')
                             for i in read_fnames]

    GrepTrimMulti('"{}"'.format(pattern_5p_rev), filtered_masked_fnames, trimmed_5p_rev_fnames,
                  trim_before=True, trim_match=True)
    
    trimmed_3p_rev_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                       postfix='filtered_trimmed_3p_rev')
                             for i in read_fnames]

    GrepTrimMulti('"{}"'.format(pattern_3p_rev), filtered_masked_fnames, trimmed_3p_rev_fnames,
                  trim_after=True, trim_match=True)

    print '##### Aligning to backbone #####'
    aligned_5p_fnames = [MakeFname(i, 'sam', dest_dir=args.tmp_dir,
                                   postfix='filtered_trimmed_5p_aligned')
                         for i in read_fnames]
    BBMapAlignMulti(bbone_db_fname, trimmed_5p_fnames, aligned_5p_fnames)
    
    aligned_3p_fnames = [MakeFname(i, 'sam', dest_dir=args.tmp_dir,
                                   postfix='filtered_trimmed_3p_aligned')
                         for i in read_fnames]
    BBMapAlignMulti(bbone_db_fname, trimmed_3p_fnames, aligned_3p_fnames)
    
    aligned_5p_rev_fnames = [MakeFname(i, 'sam', dest_dir=args.tmp_dir,
                                   postfix='filtered_trimmed_5p_rev_aligned')
                             for i in read_fnames]
    BBMapAlignMulti(bbone_db_fname, trimmed_5p_rev_fnames, aligned_5p_rev_fnames)
    
    aligned_3p_rev_fnames = [MakeFname(i, 'sam', dest_dir=args.tmp_dir,
                                       postfix='filtered_trimmed_3p_rev_aligned')
                             for i in read_fnames]
    BBMapAlignMulti(bbone_db_fname, trimmed_3p_rev_fnames, aligned_3p_rev_fnames)
    
    print '##### Indexing alignment output #####'
    logging.debug('This is aligned_5p_fnames: {}'.format(aligned_5p_fnames))
    aligned_5p_fnames_bam = SamtoolsIndexMutli(aligned_5p_fnames)
    aligned_3p_fnames_bam = SamtoolsIndexMutli(aligned_3p_fnames)
    aligned_5p_rev_fnames_bam = SamtoolsIndexMutli(aligned_5p_rev_fnames)
    aligned_3p_rev_fnames_bam = SamtoolsIndexMutli(aligned_3p_rev_fnames)
    
    print '##### Calculating insertions and writing output #####'
    
    # TODO: each of these calls to the factory is reading the same masked reads
    # again. If this is slow, should refactor into a container for the reads.
    # Otherwise, ignore.
    print 'Writing output to', args.output_fname
    out_f = open(args.output_fname, 'w')
    dict_writer = rad_factory.ReadAlignmentDataFactory.MakeDictWriter(out_f)
    
    total_matched_reads = 0 
    forward = 1
    reverse = -1
    factory = rad_factory.ReadAlignmentDataFactory(
        tn_params, '5p', forward)
    read_data_5p = factory.DictFromFileLists(
        insert_bbone_filtered_fnames, aligned_5p_fnames_bam)
    total_matched_reads += len(read_data_5p)
    factory.WriteToDictWriter(dict_writer, read_data_5p.itervalues())
    del read_data_5p  # hint to GC
    
    factory = rad_factory.ReadAlignmentDataFactory(
        tn_params, '3p', forward)
    read_data_3p = factory.DictFromFileLists(
        insert_bbone_filtered_fnames, aligned_3p_fnames_bam)
    total_matched_reads += len(read_data_3p)
    factory.WriteToDictWriter(dict_writer, read_data_3p.itervalues())
    del read_data_3p  # hint to GC
    
    factory = rad_factory.ReadAlignmentDataFactory(
        tn_params, '5p', reverse)
    read_data_5p_rev = factory.DictFromFileLists(
        insert_bbone_filtered_fnames, aligned_5p_rev_fnames_bam)
    total_matched_reads += len(read_data_5p_rev)
    factory.WriteToDictWriter(dict_writer, read_data_5p_rev.itervalues())
    del read_data_5p_rev  # hint to GC
    
    factory = rad_factory.ReadAlignmentDataFactory(
        tn_params, '3p', reverse)
    read_data_3p_rev = factory.DictFromFileLists(
        insert_bbone_filtered_fnames, aligned_3p_rev_fnames_bam)
    total_matched_reads += len(read_data_3p_rev)
    factory.WriteToDictWriter(dict_writer, read_data_3p_rev.itervalues())
    del read_data_3p_rev  # hint to GC
    out_f.close()
        
    print 'Saved', total_matched_reads, 'matching reads'
    
    duration = time.time() - start_ts
    print 'Running time: %.2f minutes' % (duration/60.0)
    
    
if __name__ == '__main__':
    Main()
