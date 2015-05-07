#!/usr/bin/python

import unittest

from Bio.Data import IUPACData
from scripts.util.ambiguous_seq import AmbiguousSequenceGenerator


class AmbiguousSeqTest(unittest.TestCase):                                         
    
    def assertAmbiguousSeqInstance(self, ambiguous, defined, n):
        i = 0
        ambiguous = ambiguous.upper()
        # n is the number of repeats
        self.assertEquals(len(ambiguous) * n, len(defined))
        for _ in xrange(n):
            for ba in ambiguous:
                allowed_bases = IUPACData.ambiguous_dna_values[ba]
                actual_base = defined[i]
                self.assertIn(actual_base, allowed_bases)
                i += 1

    def testBasic(self):
        ambiguous_seqs = ['BCT',
                          'CGBNNAT',
                          'WVWVTCWVNN',
                          'bCaaN']
        for ambig in ambiguous_seqs:
            gen = AmbiguousSequenceGenerator(ambig)
            for n in xrange(10):
                seq = gen.Generate(n=n)
                # Check that the generated sequence could have been generated 
                # from the ambiguous sequence.
                self.assertAmbiguousSeqInstance(ambig, seq, n)


if __name__ == '__main__':
    unittest.main()