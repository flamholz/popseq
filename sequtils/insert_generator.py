#!/usr/bin/python

import random

from Bio.Seq import Seq


class InsertGenerator(object):
    def __init__(self, insert_seq, fixed_5p, fixed_3p,
                 extra_bp_5p=None,
                 linker_gen=None,
                 max_linker_repeats=2):
        self.insert_seq = insert_seq
        self.fixed_5p = fixed_5p
        self.fixed_3p = fixed_3p
        self.extra_bp_5p = extra_bp_5p or ''
        self.linker_gen = linker_gen
        self.max_linker_repeats = max_linker_repeats
    
    @staticmethod
    def FromTranspositionParams(tn_params):
        return InsertGenerator(tn_params.insert_seq,
                               tn_params.fixed_seq_5p,
                               tn_params.fixed_seq_3p,
                               tn_params.extra_bp_5p,
                               tn_params.linker_pattern,
                               tn_params.max_linker_repeats)
    
    def next(self):
        linker_5p = ''
        linker_3p = ''
        if self.linker_gen:
            linker_lens = range(self.max_linker_repeats + 1)
            linker_len_5p = random.choice(linker_lens)
            linker_len_3p = random.choice(linker_lens)
            linker_5p = self.linker_gen.Generate(n=linker_len_5p)
            linker_3p = self.linker_gen.Generate(n=linker_len_3p)
        
        output_seq = '%s%s%s%s%s%s' % (self.fixed_5p, self.extra_bp_5p,
                                       linker_5p, self.insert_seq, linker_3p,
                                       self.fixed_3p)
        
        return Seq(output_seq), linker_5p, linker_3p