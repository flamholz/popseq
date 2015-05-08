#!/usr/bin/python

import unittest
import glob
import numpy as np
import sequtils.read_alignment_data.factory as rad_factory

from sequtils import read_alignment_data as rad
from sequtils.synthetic_transposon import Fragment



class ReadAligmentDataTest(unittest.TestCase):
    MASKED_READS_FNAME = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_insert_masked.fq'
    ALIGNED_3P = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_3p_aligned.bam'
    ALIGNED_5P = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_5p_aligned.bam'
    ALIGNED_3P_REV = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_3p_rev_aligned.bam'
    ALIGNED_5P_REV = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_5p_rev_aligned.bam'
    START_OFFSET = 23
    FIXED_5P = "TGCATC"
    FIXED_3P = "GCGTCA"
    FORWARD = 1
    REVERSE = -1

    def GenericTest(self, fixed_seq, end, orientation,
                    aligned_fname):
        factory = rad_factory.ReadAlignmentDataFactory(
            self.START_OFFSET, fixed_seq, end, orientation)
        rads_by_id = factory.DictFromFiles(self.MASKED_READS_FNAME,
                                           aligned_fname)
        for rad in rads_by_id.itervalues():
            read_info = Fragment.ParseInfoDict(rad.read_record)
            expected_insertion_site = read_info['ins']
            self.assertEquals(expected_insertion_site, rad.insertion_site)

    def testAligned5p(self):
        self.GenericTest(self.FIXED_5P, '5p', self.FORWARD, self.ALIGNED_5P)

    def testAligned3p(self):
        self.GenericTest(self.FIXED_3P, '3p', self.FORWARD, self.ALIGNED_3P)

    def testAligned5pReverse(self):
        self.GenericTest(self.FIXED_5P, '5p', self.REVERSE, self.ALIGNED_5P_REV)

    def testAligned3pReverse(self):
        self.GenericTest(self.FIXED_3P, '3p', self.REVERSE, self.ALIGNED_3P_REV)


if __name__ == '__main__':
    unittest.main()