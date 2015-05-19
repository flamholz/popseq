#!/usr/bin/python

import unittest

from Bio.Data import IUPACData
from sequtils.ambiguous_seq import AmbiguousSequence


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

    def testEmpty(self):
        gen = AmbiguousSequence('')
        for n in xrange(20):
            self.assertEquals('', str(gen.Generate(n)))

    def testBasic(self):
        ambiguous_seqs = ['BCT',
                          'CGBNNAT',
                          'WVWVTCWVNN',
                          'bCaaN']
        for ambig in ambiguous_seqs:
            gen = AmbiguousSequence(ambig)
            for n in xrange(10):
                seq = gen.Generate(n=n)
                # Check that the generated sequence could have been generated 
                # from the ambiguous sequence.
                self.assertAmbiguousSeqInstance(ambig, seq, n)
                self.assertEquals(n, gen.IsInstance(seq))

    def testFailure(self):
        # These sequences are not instances of the ambiguous key sequence.
        ambiguous_seqs = {'BCT': ['TTT', 'AGA', 'CCTTCTACT', 'GCTCC'],
                          'CGBNNAT': ['cagtagt', 'CTCAAATC']}
        for ambig, defineds in ambiguous_seqs.iteritems():
            gen = AmbiguousSequence(ambig)
            self.assertEquals(0, gen.IsInstance(''))
            for should_fail in defineds:
                self.assertEquals(-1, gen.IsInstance(should_fail))


if __name__ == '__main__':
    unittest.main()