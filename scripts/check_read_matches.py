#!/usr/bin/python

import argparse
import numpy as np
import json
import pandas as pd
import pylab

from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from sequtils import read_alignment_data as rad
from sequtils.read_alignment_data import factory



def _parseReadInfo(rr):
    """Parses JSON encoded info stashed in read description."""
    desc = rr.description
    dict_start = desc.find('{')
    desc = desc[dict_start:]
    return json.loads(desc)


def Main():
    parser = argparse.ArgumentParser(description='Filter reads.', fromfile_prefix_chars='@')
    parser.add_argument("-i", "--insert_alignment_fname", required=True,
                        help=("Path to FASTA file containing insert ends to align reads to. "
                              "Will only retain reads that align well to this DB."))
    parser.add_argument("-b", "--backbone_alignment_fname", required=True,
                        help=("Path to FASTA file containing backbone sequence. "
                              "Will bin reads by where they align to this sequence."))
    parser.add_argument("-r", "--reads_fname", required=True,
                        help="Path to FASTQ files containing reads.")
    parser.add_argument("--start_offset", default=23, type=int,
                        help=("Offset of the start codon into the sequence used "
                              "to match the backbone (nt units)"))
    args = parser.parse_args()
    
    ofh = open('false_negative_reads.fa', 'w')
    writer = FastaIO.FastaWriter(ofh)
    writer.write_header()
    
    d_by_construct = {}
    
    FACTORY = factory.ReadAlignmentDataFactory(backbone_start_offset=args.start_offset,
                                               fixed_5p_seq=rad.DEFAULT_FIXED_5P_SEQ,
                                               fixed_3p_seq=rad.DEFAULT_FIXED_3P_SEQ)
    print args.reads_fname
    READ_DATA = FACTORY.DictFromFiles(args.insert_alignment_fname,
                                      args.backbone_alignment_fname,
                                      args.reads_fname)
    with open(args.reads_fname) as f:
        reader = SeqIO.parse(f, 'fasta')
        for record in reader:
            read_info = _parseReadInfo(record)
            
            cnum = read_info["construct_num"]
            d = d_by_construct.setdefault(cnum, read_info)
            
            found_match = (record.id in READ_DATA)
            should_match = read_info['should_match']
            true_pos = found_match and should_match
            false_pos = found_match and not should_match
            true_neg = not found_match and not should_match
            false_neg = not found_match and should_match
            
            d['total_reads'] = d.get('total_reads', 0) + 1 
            d['true_pos'] = d.get('true_pos', 0) + true_pos
            d['false_pos'] = d.get('false_pos', 0) + false_pos
            d['true_neg'] = d.get('true_neg', 0) + true_neg
            d['false_neg'] = d.get('false_neg', 0) + false_neg
            
            if false_neg:
                writer.write_record(record)
    writer.write_footer()
    ofh.close()

    df = pd.DataFrame.from_dict(d_by_construct, orient='index')
    df.to_csv('match_stats.csv')


if __name__ == '__main__':
    Main()