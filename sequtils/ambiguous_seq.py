#!/usr/bin/python

import random

from Bio.Data import IUPACData
from Bio.Seq import Seq


class AmbiguousSequenceGenerator(object):
    """Generates a concrete sequence at random from an ambiguous one."""
    def __init__(self, pattern):
        """Initialize.
        
        Args:
            pattern: a DNA sequence including ambiguous bases.
        """
        self.pattern = pattern.upper()
    
    @staticmethod
    def _RandomNucleotide(b):
        """Picks a random nucleotide from those b could code for.
        
        Args:
            b: a base which may or may not be ambiguous.
        
        Returns:
            A single base that could be coded for by b. If b is not
            ambiguous, will return b. Otherwise, picks one of the
            bases b could code for uniformly at random. 
        """
        possible_bases = IUPACData.ambiguous_dna_values[b]
        return random.choice(possible_bases)
    
    def Generate(self, n=1):
        """Generate a random sequence from the pattern with n repeats.
        
        Args:
            n: the number of times to repeat the pattern.
        
        Returns:
            A concrete (unambiguous) DNA sequence.
        """
        out = []
        for _ in xrange(n):
            out.extend([self._RandomNucleotide(b) for b in self.pattern])
        return Seq(''.join(out))
    
    def IsInstance(self, defined):
        """Returns > 0 if the defined sequence is an instance of this
            ambiguous one.
        
        Args:
            defined: a defined DNA sequence with no ambiguous nt.
        
        Returns:
            A number N equal to the number of repeats of the ambiguous template
            in the defined sequence.
        """
        defined = defined.upper()
        n = 0
        i = 0
        while i < len(defined):
            for ba in self.pattern:
                bd = defined[i]
                allowed_bases = IUPACData.ambiguous_dna_values[ba]
                if bd not in allowed_bases:
                    return 0
                i += 1
            n += 1
        return n
