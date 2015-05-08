#!/usr/bin/python

import argparse
import itertools
import os
import subprocess
import time
import sequtils.read_alignment_data.factory as rad_factory

from Bio.Seq import Seq
from os import path
from scripts.util.filename_util import MakeFname
from scripts.util import command_util


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
               'samtools sort - %s' % fname]
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
    parser = argparse.ArgumentParser(description='Filter reads.', fromfile_prefix_chars='@')
    parser.add_argument("-i", "--insert_db_filename", required=True,
                        help=("Path to FASTA file containing insert ends to align reads to. "
                              "Will only retain reads that align well to this DB."))
    parser.add_argument("-b", "--backbone_db_filename", required=True,
                        help=("Path to FASTA file containing backbone sequence. "
                              "Will bin reads by where they align to this sequence."))
    parser.add_argument("-r", "--read_filenames", nargs='+', required=True,
                        help="Path to FASTQ files containing reads.")
    parser.add_argument("--start_offset", required=True, type=int,
                        help=("Offset of the start codon into the sequence used "
                              "to match the backbone (nt units)"))
    parser.add_argument("--tn_bp_duplicated", required=True, type=int,
                        help="Number of bases duplicated by transposition (nt units)")
    parser.add_argument("-t", "--tmp_dir",
                        default="_read_filter_data",
                        help="Path to use to store intermediate files.")
    parser.add_argument("-o", "--output_fname", required=True,
                        help="Where to write output data (CSV).")
    parser.add_argument("--fixed_5p",
                        default="TGCATC",
                        help="Fixed sequence found on 5' end of insert.")
    parser.add_argument("--fixed_3p",
                        default="GCGTCA",
                        help="Fixed sequence found on 3' end of insert.")
    parser.set_defaults(summary_output=True)
    args = parser.parse_args()
    
    start_ts = time.time()
    command_util.CheckAllInstalled(['fastq-grep', 'bbduk.sh', 'bowtie2', 'samtools'])

    insert_db_fname = args.insert_db_filename
    bbone_db_fname = args.backbone_db_filename
    read_fnames = args.read_filenames

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
    GrepTrimMulti(pattern_5p, filtered_masked_fnames, trimmed_5p_fnames,
                  trim_after=True, trim_match=True)
    
    trimmed_3p_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                   postfix='filtered_trimmed_3p')
                         for i in read_fnames]
    GrepTrimMulti(pattern_3p, filtered_masked_fnames, trimmed_3p_fnames,
                  trim_before=True, trim_match=True)
    
    trimmed_5p_rev_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                       postfix='filtered_trimmed_5p_rev')
                             for i in read_fnames]
    GrepTrimMulti(pattern_5p_rev, filtered_masked_fnames, trimmed_5p_rev_fnames,
                  trim_before=True, trim_match=True)
    
    trimmed_3p_rev_fnames = [MakeFname(i, 'fq', dest_dir=args.tmp_dir,
                                       postfix='filtered_trimmed_3p_rev')
                             for i in read_fnames]
    GrepTrimMulti(pattern_3p_rev, filtered_masked_fnames, trimmed_3p_rev_fnames,
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
    aligned_5p_fnames_bam = SamtoolsIndexMutli(aligned_5p_fnames)
    aligned_3p_fnames_bam = SamtoolsIndexMutli(aligned_3p_fnames)
    aligned_5p_rev_fnames_bam = SamtoolsIndexMutli(aligned_5p_rev_fnames)
    aligned_3p_rev_fnames_bam = SamtoolsIndexMutli(aligned_3p_rev_fnames)
    
    print '##### Evaluating backbone alignment #####'
    
    n_eq = 0
    total = 0
    
    reads_by_id = {}
    matching_ids = set()
    
    
    # TODO: each of these calls to the factory is reading the same masked reads
    # again. If this is slow, should refactor into a container for the reads.
    # Otherwise, ignore. 
    forward = 1
    reverse = -1
    factory = rad_factory.ReadAlignmentDataFactory(
        args.start_offset, fixed_5p, '5p', forward)
    read_data_5p = factory.DictFromFileLists(
        filtered_masked_fnames, aligned_5p_fnames_bam)
    
    factory = rad_factory.ReadAlignmentDataFactory(
        args.start_offset, fixed_3p, '3p', forward)
    read_data_3p = factory.DictFromFileLists(
        filtered_masked_fnames, aligned_3p_fnames_bam)
    
    factory = rad_factory.ReadAlignmentDataFactory(
        args.start_offset, fixed_5p, '5p', reverse)
    read_data_5p_rev = factory.DictFromFileLists(
        filtered_masked_fnames, aligned_5p_rev_fnames_bam)
    
    factory = rad_factory.ReadAlignmentDataFactory(
        args.start_offset, fixed_3p, '3p', reverse)
    read_data_3p_rev = factory.DictFromFileLists(
        filtered_masked_fnames, aligned_3p_rev_fnames_bam)
    
    out_fname = args.output_fname
    print 'Writing output to', out_fname
    all_matched_reads = itertools.chain(read_data_5p.itervalues(),
                                        read_data_3p.itervalues(),
                                        read_data_5p_rev.itervalues(),
                                        read_data_3p_rev.itervalues())
    factory.WriteCSVFilename(all_matched_reads, out_fname)
    
    """
    for bam_fname in aligned_5p_fnames_bam:
        alignedf = pysam.AlignmentFile(bam_fname, 'rb')
        aligned = alignedf.fetch()
        for i, it in enumerate(aligned):
            match_3p_end = it.reference_end
            if it.is_reverse:
                match_3p_end = it.reference_start + args.tn_bp_duplicated
        
            insertion_idx = match_3p_end
            insertion_site = insertion_idx - args.start_offset
            try:
                read_info = _parseReadInfo(it.query_name)
            except:
                print 'Error parsing', it
                continue
            
            expected_insertion_site = read_info["ins"]
            read_info["cins"] = insertion_site
            n_eq += (expected_insertion_site == insertion_site)            
            total += 1
            
            read_id = ReadID(read_info)
            read_info['rid'] = read_id
            if expected_insertion_site == insertion_site:
                if read_id in matching_ids:
                    print 'AAA 5p'
                matching_ids.add(read_id)
            reads_by_id[read_id] = read_info
            
    for bam_fname in aligned_3p_fnames_bam:
        alignedf = pysam.AlignmentFile(bam_fname, 'rb')
        aligned = alignedf.fetch()
        for i, it in enumerate(aligned):
            match_5p_end = it.reference_start + args.tn_bp_duplicated
            if it.is_reverse:
                match_5p_end = it.reference_end
            
            insertion_idx = match_5p_end
            insertion_site = insertion_idx - args.start_offset
            
            read_info = _parseReadInfo(it.query_name)
            expected_insertion_site = read_info["ins"]
            read_info["cins"] = insertion_site
            
            n_eq += (expected_insertion_site == insertion_site)            
            total += 1
            
            read_id = ReadID(read_info)
            read_info['rid'] = read_id
            if expected_insertion_site == insertion_site:
                if read_id in matching_ids:
                    print 'AAA 3p'
                matching_ids.add(read_id)
            reads_by_id[read_id] = read_info

    for bam_fname in aligned_5p_rev_fnames_bam:
        alignedf = pysam.AlignmentFile(bam_fname, 'rb')
        aligned = alignedf.fetch()
        for i, it in enumerate(aligned):
            match_5p_end = it.reference_start + args.tn_bp_duplicated
            if it.is_reverse:
                match_5p_end = it.reference_end
            
            insertion_idx = match_5p_end
            insertion_site = insertion_idx - args.start_offset
            
            read_info = _parseReadInfo(it.query_name)
            expected_insertion_site = read_info["ins"]
            read_info["cins"] = insertion_site
            n_eq += (expected_insertion_site == insertion_site)            
            total += 1
            
            read_id = ReadID(read_info)
            read_info['rid'] = read_id
            if expected_insertion_site == insertion_site:
                if read_id in matching_ids:
                    print 'AAA 5p rev'
                matching_ids.add(read_id)
            reads_by_id[read_id] = read_info

    for bam_fname in aligned_3p_rev_fnames_bam:
        alignedf = pysam.AlignmentFile(bam_fname, 'rb')
        aligned = alignedf.fetch()
        for i, it in enumerate(aligned):
            match_5p_end = it.reference_end
            if it.is_reverse:
                match_5p_end = it.reference_start + args.tn_bp_duplicated
            
            insertion_idx = match_5p_end
            insertion_site = insertion_idx - args.start_offset
            
            read_info = _parseReadInfo(it.query_name)
            expected_insertion_site = read_info["ins"]
            read_info["cins"] = insertion_site
            
            n_eq += (expected_insertion_site == insertion_site)            
            total += 1
            
            read_id = ReadID(read_info)
            read_info['rid'] = read_id
            if expected_insertion_site == insertion_site:
                if read_id in matching_ids:
                    print 'AAA 3p rev'
                matching_ids.add(read_id)
               
            reads_by_id[read_id] = read_info
    
    print '%d correct insertion sites' % n_eq
    print '%d total reads mapped' % total
    pct = 100.0 * float(n_eq) / float(total)
    print '%.2f%% pct of total' % pct
    
    
    for reads_fname in insert_bbone_filtered_fnames:
        reader = SeqIO.parse(reads_fname, 'fastq')
        for record in reader:
            read_info = _parseReadInfo(record.description)
            read_id = ReadID(read_info)
            read_info["cins"] = None
            read_info['rid'] = read_id
            if read_id not in reads_by_id:
                reads_by_id[read_id] = read_info
    
    print '##### Writing output file #####'
    fieldnames = sorted(read_info.keys())
    print 'total reads examined', len(reads_by_id)
    print 'total matching', len(matching_ids)
    out_csv_fname = os.path.join(args.tmp_dir, args.output_fname)
    with open(out_csv_fname, 'w') as outf:
        writer = csv.DictWriter(outf, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(reads_by_id.itervalues())
    """
    duration = time.time() - start_ts
    print 'Running time: %.2f minutes' % (duration/60.0)
    
    
if __name__ == '__main__':
    Main()