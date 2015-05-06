#!/usr/bin/python

import random

from Bio.Data import IUPACData
from Bio.Seq import Seq


class AmbiguousSequenceGenerator(object):
    def __init__(self, pattern):
        self.pattern = pattern.upper()
    
    @staticmethod
    def _RandomNucleotide(b):
        possible_bases = IUPACData.ambiguous_dna_values[b]
        return random.choice(possible_bases)
    
    def Generate(self, n=1):
        out = []
        for _ in xrange(n):
            out.extend([self._RandomNucleotide(b) for b in self.pattern])
        return Seq(''.join(out))
