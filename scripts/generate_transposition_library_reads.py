#!/usr/bin/env python

import argparse

from Bio.SeqIO import FastaIO
from sequtils.insert_generator import InsertGenerator
from sequtils.synthetic_transposon import Transposition
from sequtils.transposition_params import TranspositionParams
from zipfile import ZipFile


def Main():
    parser = argparse.ArgumentParser(description='Generate synthetic reads.',
                                     fromfile_prefix_chars='@')
    parser.add_argument("-t", "--num_transpositions", type=int, default=50,
                        help="Number of transpositions to generate")
    parser.add_argument("-r", "--num_reads", type=int, default=10000,
                        help="Number of reads to generate per transposition.")
    parser.add_argument("-l", "--read_length", type=int, default=100,
                        help="Number of reads to generate per transposition.")
    parser.add_argument("-o", "--output_fname", default='generated_transposition_reads.fa',
                        help="Where to write FASTA output to.")
    TranspositionParams.AddArgs(parser)
    args = parser.parse_args()
    
    tn_params = TranspositionParams.FromArgs(args)
    insert_gen = InsertGenerator.FromTranspositionParams(tn_params)
    
    n_trans = args.num_transpositions
    n_reads = args.num_reads
    read_len = args.read_length
    print 'Generating %d random transpositions with %d random %d NT reads each' % (
                n_trans, n_reads, read_len)
    
    print 'Writing generated reads to FASTA'                                          
    x = 0
    with open(args.output_fname, 'w') as fh:
        writer = FastaIO.FastaWriter(fh)
        writer.write_header()
        for construct_num in xrange(n_trans):
            trans = Transposition(construct_num, insert_gen,
                                  tn_params.backbone_seq,
                                  tn_params.backbone_start_offset)
            for read_num in xrange(n_reads): 
                frag = trans.Shear(read_num, read_len)
                record = frag.ToSeqRecord()
                writer.write_record(record)
                
                x += 1
                if x % 1000000 == 0:
                    print "Created %d reads" % x
                    
        writer.write_footer()
    
    print 'Zipping generated reads.'
    zip_fname = '%s.zip' % args.output_fname
    with ZipFile(zip_fname, mode='w') as zipf:
        zipf.write(args.output_fname)


if __name__ == '__main__':
    Main()