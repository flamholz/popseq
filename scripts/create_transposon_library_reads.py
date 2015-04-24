#!/usr/bin/env python

import argparse
import json

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
from Bio.Alphabet import DNAAlphabet
from random import randrange
from random import random


# Target DNA region containing gene of iterest
target_cas9 = Seq('GATCTAAAGAGGAGAAAGGATCTATGGATAAGAAATACTCAATAGGCTTAGCTATCGGCACAAATAGCGTCGGATGGGCGGTGATCACTGATGAATATAAGGTTCCGTCTAAAAAGTTCAAGGTTCTGGGAAATACAGACCGCCACAGTATCAAAAAAAATCTTATAGGGGCTCTTTTATTTGACAGTGGAGAGACAGCGGAAGCGACTCGTCTCAAACGGACAGCTCGTAGAAGGTATACACGTCGGAAGAATCGTATTTGTTATCTACAGGAGATTTTTTCAAATGAGATGGCGAAAGTAGATGATAGTTTCTTTCATCGACTTGAAGAGTCTTTTTTGGTGGAAGAAGACAAGAAGCATGAACGTCATCCTATTTTTGGAAATATAGTAGATGAAGTTGCTTATCATGAGAAATATCCAACTATCTATCATCTGCGAAAAAAATTGGTAGATTCTACTGATAAAGCGGATTTGCGCTTAATCTATTTGGCCTTAGCGCATATGATTAAGTTTCGTGGTCATTTTTTGATTGAGGGAGATTTAAATCCTGATAATAGTGATGTGGACAAACTATTTATCCAGTTGGTACAAACCTACAATCAATTATTTGAAGAAAACCCTATTAACGCAAGTGGAGTAGATGCTAAAGCGATTCTTTCTGCACGATTGAGTAAATCAAGACGATTAGAAAATCTCATTGCTCAGCTCCCCGGTGAGAAGAAAAATGGCTTATTTGGGAATCTCATTGCTTTGTCATTGGGTTTGACCCCTAATTTTAAATCAAATTTTGATTTGGCAGAAGATGCTAAATTACAGCTTTCAAAAGATACTTACGATGATGATTTAGATAATTTATTGGCGCAAATTGGAGATCAATATGCTGATTTGTTTTTGGCAGCTAAGAATTTATCAGATGCTATTTTACTTTCAGATATCCTAAGAGTAAATACTGAAATAACTAAGGCTCCCCTATCAGCTTCAATGATTAAACGCTACGATGAACATCATCAAGACTTGACTCTTTTAAAAGCTTTAGTTCGACAACAACTTCCAGAAAAGTATAAAGAAATCTTTTTTGATCAATCAAAAAACGGATATGCAGGTTATATTGATGGGGGAGCTAGCCAAGAAGAATTTTATAAATTTATCAAACCAATTTTAGAAAAAATGGATGGTACTGAGGAATTATTGGTGAAACTAAATCGTGAAGATTTGCTGCGCAAGCAACGGACCTTTGACAACGGCTCTATTCCCCATCAAATTCACTTGGGTGAGCTGCATGCTATTTTGAGAAGACAAGAAGACTTTTATCCATTTTTAAAAGACAATCGTGAGAAGATTGAAAAAATCTTGACTTTTCGAATTCCTTATTATGTTGGTCCATTGGCGCGTGGCAATAGTCGTTTTGCATGGATGACTCGGAAGTCTGAAGAAACAATTACCCCATGGAATTTTGAAGAAGTTGTCGATAAAGGTGCTTCAGCTCAATCATTTATTGAACGCATGACAAACTTTGATAAAAATCTTCCAAATGAAAAAGTACTACCAAAACATAGTTTGCTTTATGAGTATTTTACGGTTTATAACGAATTGACAAAGGTCAAATATGTTACTGAAGGAATGCGAAAACCAGCATTTCTTTCAGGTGAACAGAAGAAAGCCATTGTTGATTTACTCTTCAAAACAAATCGAAAAGTAACCGTTAAGCAATTAAAAGAAGATTATTTCAAAAAAATAGAATGTTTTGATAGTGTTGAAATTTCAGGAGTTGAAGATAGATTTAATGCTTCATTAGGTACCTACCATGATTTGCTAAAAATTATTAAAGATAAAGATTTTTTGGATAATGAAGAAAATGAAGATATCTTAGAGGATATTGTTTTAACATTGACCTTATTTGAAGATAGGGAGATGATTGAGGAAAGACTTAAAACATATGCTCACCTCTTTGATGATAAGGTGATGAAACAGCTTAAACGTCGCCGTTATACTGGTTGGGGACGTTTGTCTCGAAAATTGATTAATGGTATTAGGGATAAGCAATCTGGCAAAACAATATTAGATTTTTTGAAATCAGATGGTTTTGCCAATCGCAATTTTATGCAGCTGATCCATGATGATAGTTTGACATTTAAAGAAGACATTCAAAAAGCACAAGTGTCTGGACAAGGCGATAGTTTACATGAACATATTGCAAATTTAGCTGGTAGCCCTGCTATTAAAAAAGGTATTTTACAGACTGTAAAAGTTGTTGATGAATTGGTCAAAGTAATGGGGCGGCATAAGCCAGAAAATATCGTTATTGAAATGGCACGTGAAAATCAGACAACTCAAAAGGGCCAGAAAAATTCGCGAGAGCGTATGAAACGAATCGAAGAAGGTATCAAAGAATTAGGAAGTCAGATTCTTAAAGAGCATCCTGTTGAAAATACTCAATTGCAAAATGAAAAGCTCTATCTCTATTATCTCCAAAATGGAAGAGACATGTATGTGGACCAAGAATTAGATATTAATCGTTTAAGTGATTATGATGTCGATGCCATTGTTCCACAAAGTTTCCTTAAAGACGATTCAATAGACAATAAGGTCTTAACGCGTTCTGATAAAAATCGTGGTAAATCGGATAACGTTCCAAGTGAAGAAGTAGTCAAAAAGATGAAAAACTATTGGAGACAACTTCTAAACGCCAAGTTAATCACTCAACGTAAGTTTGATAATTTAACGAAAGCTGAACGTGGAGGTTTGAGTGAACTTGATAAAGCTGGTTTTATCAAACGCCAATTGGTTGAAACTCGCCAAATCACTAAGCATGTGGCACAAATTTTGGATAGTCGCATGAATACTAAATACGATGAAAATGATAAACTTATTCGAGAGGTTAAAGTGATTACCTTAAAATCTAAATTAGTTTCTGACTTCCGAAAAGATTTCCAATTCTATAAAGTACGTGAGATTAACAATTACCATCATGCCCATGATGCGTATCTAAATGCCGTCGTTGGAACTGCTTTGATTAAGAAATATCCAAAACTTGAATCGGAGTTTGTCTATGGTGATTATAAAGTTTATGATGTTCGTAAAATGATTGCTAAGTCTGAGCAAGAAATAGGCAAAGCAACCGCAAAATATTTCTTTTACTCTAATATCATGAACTTCTTCAAAACAGAAATTACACTTGCAAATGGAGAGATTCGCAAACGCCCTCTAATCGAAACTAATGGGGAAACTGGAGAAATTGTCTGGGATAAAGGGCGAGATTTTGCCACAGTGCGCAAAGTATTGTCCATGCCCCAAGTCAATATTGTCAAGAAAACAGAAGTACAGACAGGCGGATTCTCCAAGGAGTCAATTTTACCAAAAAGAAATTCGGACAAGCTTATTGCTCGTAAAAAAGACTGGGATCCAAAAAAATATGGTGGTTTTGATAGTCCAACGGTAGCTTATTCAGTCCTAGTGGTTGCTAAGGTGGAAAAAGGGAAATCGAAGAAGTTAAAATCCGTTAAAGAGTTACTAGGGATCACAATTATGGAAAGAAGTTCCTTTGAAAAAAATCCGATTGACTTTTTAGAAGCTAAAGGATATAAGGAAGTTAAAAAAGACTTAATCATTAAACTACCTAAATATAGTCTTTTTGAGTTAGAAAACGGTCGTAAACGGATGCTGGCTAGTGCCGGAGAATTACAAAAAGGAAATGAGCTGGCTCTGCCAAGCAAATATGTGAATTTTTTATATTTAGCTAGTCATTATGAAAAGTTGAAGGGTAGTCCAGAAGATAACGAACAAAAACAATTGTTTGTGGAGCAGCATAAGCATTATTTAGATGAGATTATTGAGCAAATCAGTGAATTTTCTAAGCGTGTTATTTTAGCAGATGCCAATTTAGATAAAGTTCTTAGTGCATATAACAAACATAGAGACAAACCAATACGTGAACAAGCAGAAAATATTATTCATTTATTTACGTTGACGAATCTTGGAGCTCCCGCTGCTTTTAAATATTTTGATACAACAATTGATCGTAAACGATATACGTCTACAAAAGAAGTTTTAGATGCCACTCTTATCCATCAATCCATCACTGGTCTTTATGAAACACGCATTGATTTGAGTCAGCTAGGAGGTGACTAACTCGA', DNAAlphabet())

# First base of cas9 (A of ATG) is base 24 (1 indexed)
orf_start = 24

# Inserting transposon pre-cloned with domain of interest, no linker variation
trans_pdz = Seq('TGCATCTCAACGTCGGCGTGTGACGGTGCGCAAGGCCGACGCCGGCGGGCTGGGCATCAGCATCAAGGGGGGCCGGGAAAACAAGATGCCTATTCTCATTTCCAAAATCTTCAAGGGACTGGCAGCAGACCAGACGGAGGCCCTTTTTGTTGGGGATGCCATCCTGTCTGTGAATGGTGAAGATTTGTCCTCTGCCACCCACGATGAAGCGGTACAGGCCCTCAAGAAGACAGGCAAGGAGGTCGTGCTCGAAGTTAAGTACATGGCGTCA', DNAAlphabet())
pdz_len = len(trans_pdz)


class Fragment(object):
    
    def __init__(self, id_num, trans, frag_start, frag_end, forward):
        self.id = id_num
        self.trans = trans # Transposition from which derived.
        self.start = frag_start
        self.end = frag_end
        
        self.seq = trans.construct[frag_start:frag_end] # Fragment sequence
        self.forward = forward
        if not self.forward:
            self.seq = self.seq.reverse_complement()
        
        inlen = trans.insert_len
        ins = trans.insertion_site
        extbp = 20
        efs = frag_start + extbp
        efe = frag_end - extbp
        
        should_match_5p = efs <= ins <= efe
        should_match_3p = efs  <= ins + inlen <= efe
        self.should_match = (should_match_3p or should_match_5p)
        id_tuple = (self.trans.id, self.id, self.trans.insertion_site)
        self.id_str = "@C%06dR%04dINSPS%04d\n" % id_tuple
        self.info_dict = {"insertion_site": self.trans.expected_insertion_site,
                          "construct_num": self.trans.id,
                          "read_num": self.id,
                          "should_match": self.should_match}


class Transposition(object):
    
    def __init__(self, id_num, insert_seq, target_seq, target_orf_start):
        self.id = id_num
        self.insert = insert_seq
        self.insert_len = len(insert_seq)
        self.target = target_seq
        self.orf_start = target_orf_start
        self.forward_insertion = True
        
        # Random insert position within target sequence
        insert = self.insert
        target = self.target
        # Insert site relative to start of backbone sequence, not start codon.
        self.insertion_site = randrange(0, len(self.target)-5+1)
        
        # 50% of inserts are reverse
        if random() > 0.5:
            insert = insert.reverse_complement()
            self.forward_insertion = False
            
        # New sequence with insert tranposed and 5bp repeat created and ORF insert site calculated
        ins = self.insertion_site
        self.construct = target[:ins+5] + insert + target[ins:]
        # Expected insertion site 5' most base of insert sequence relative to start codon.
        self.expected_insertion_site = ins + 5 - (target_orf_start-1)
    
    def Shear(self, frag_id, fragment_length=100):
        shear_start = randrange(0, len(self.construct)-fragment_length+1)
        shear_end = shear_start + fragment_length
        fwd = random() <= 0.5
        return Fragment(frag_id, self, shear_start, shear_end, fwd)


def Main():
    parser = argparse.ArgumentParser(description='Filter reads.', fromfile_prefix_chars='@')
    parser.add_argument("-t", "--num_transpositions", type=int, default=50,
                        help="Number of transpositions to generate")
    parser.add_argument("-r", "--num_reads", type=int, default=10000,
                        help="Number of reads to generate per transposition.")
    parser.add_argument("-l", "--read_length", type=int, default=100,
                        help="Number of reads to generate per transposition.")
    parser.add_argument("-o", "--output_fname", default='generated_transposition_reads.fa',
                        help="Where to write FASTA output to.")
    
    parser.set_defaults(summary_output=True)
    args = parser.parse_args()
    x = 0
    with open(args.output_fname, 'w') as fh:
        records = []
        writer = FastaIO.FastaWriter(fh)
        writer.write_header()
        for construct_num in xrange(args.num_transpositions):
            trans = Transposition(construct_num, trans_pdz, target_cas9, orf_start)
            for read_num in xrange(args.num_reads): 
                frag = trans.Shear(read_num, args.read_length)
                construct_desc = json.dumps(frag.info_dict)
                record = SeqRecord(frag.seq, id=frag.id_str,
                                   name=frag.id_str,
                                   description=construct_desc)
                writer.write_record(record)
                
                x += 1
                if x % 1000000 == 0:
                    print "Created %d reads" % x
                    
        writer.write_footer()


if __name__ == '__main__':
    Main()