#!/usr/bin/python

import unittest

from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from scripts.util.insert_generator import InsertGenerator
from scripts.util.ambiguous_seq import AmbiguousSequenceGenerator


class AmbiguousSeqTest(unittest.TestCase):
    INSERT_SEQ = Seq('CAACGTCGGCGTGTGACGGTGCGCAGGTCGTGCTCGAAGTTAAGTACATG', DNAAlphabet())
    FIXED_5P = Seq('TGCATC' + 'T')
    FIXED_3P = Seq('GCGTCA')

    def testNoLinker(self):
        gen = InsertGenerator(self.INSERT_SEQ, self.FIXED_5P, self.FIXED_3P)
        expected = self.FIXED_5P + self.INSERT_SEQ + self.FIXED_3P
        expected_str = str(expected)
        for _ in xrange(50):
            generated, linker_5p, linker_3p = gen.next()
            # Need to compare as strings b/c metadata not equal.
            self.assertEquals(expected_str, str(generated))
            self.assertEquals('', str(linker_5p))
            self.assertEquals('', str(linker_3p))
            
    def testLinker(self):
        linker_seq = 'BCT'
        linker_gen = AmbiguousSequenceGenerator(linker_seq)
        gen = InsertGenerator(self.INSERT_SEQ, self.FIXED_5P, self.FIXED_3P, linker_gen)
        
        insert_str = str(self.INSERT_SEQ)
        fixed_5p_str = str(self.FIXED_5P)
        fixed_3p_str = str(self.FIXED_3P)
        
        for _ in xrange(50):
            generated, linker_5p, linker_3p = gen.next()
            generated_str = str(generated)
            
            self.assertTrue(generated_str.startswith(fixed_5p_str))
            self.assertTrue(generated_str.endswith(fixed_3p_str))
            
            start_idx = generated_str.find(insert_str[:12])
            end_idx = generated_str.find(fixed_3p_str)
            linker_5p_nt = start_idx - len(self.FIXED_5P)
            linker_3p_nt = end_idx - start_idx - len(insert_str)
             
            self.assertGreater(start_idx, -1)
            self.assertGreater(end_idx, start_idx)
            
            # Linkers should have integer number of codons
            self.assertEquals(linker_5p_nt % 3, 0)
            self.assertEquals(linker_3p_nt % 3, 0)
            self.assertEquals(len(linker_5p), linker_5p_nt)
            self.assertEquals(len(linker_3p), linker_3p_nt)
            

if __name__ == '__main__':
    unittest.main()