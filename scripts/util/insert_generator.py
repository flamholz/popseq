#!/usr/bin/python

import random


class InsertGenerator(object):
    def __init__(self, insert_seq, fixed_5p, fixed_3p, linker_gen=None):
        self.insert_seq = insert_seq
        self.fixed_5p = fixed_5p
        self.fixed_3p = fixed_3p
        self.linker_gen = linker_gen
    
    def next(self):
        linker_5p = ''
        linker_3p = ''
        if self.linker_gen:
            linker_lens = range(4)
            linker_len_5p = random.choice(linker_lens)
            linker_len_3p = random.choice(linker_lens)
            linker_5p = self.linker_gen.Generate(n=linker_len_5p)
            linker_3p = self.linker_gen.Generate(n=linker_len_3p)
        
        output_seq = self.fixed_5p + linker_5p + self.insert_seq + linker_3p + self.fixed_3p
        
        return output_seq, linker_5p, linker_3p