#!/usr/bin/python

import unittest
import sequtils.read_alignment_data.factory as rad_factory

from StringIO import StringIO
from sequtils.transposition_params import TranspositionParams


class FactoryTest(unittest.TestCase):
    MASKED_READS_FNAME = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_filtered_insert_masked.fq'
    ALIGNED_3P = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_filtered_trimmed_3p_aligned.bam'
    ALIGNED_5P = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_filtered_trimmed_5p_aligned.bam'
    ALIGNED_3P_REV = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_filtered_trimmed_3p_rev_aligned.bam'
    ALIGNED_5P_REV = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_filtered_trimmed_5p_rev_aligned.bam'
    
    TN_PARAMS = TranspositionParams(insert_seq=None,
                                    backbone_seq=None,
                                    backbone_start_offset=23,
                                    fixed_seq_5p="TGCATC",
                                    fixed_seq_3p="GCGTCA",
                                    tn_bp_duplicated=5)
    FORWARD = 1
    REVERSE = -1

    def GenericTest(self, end, orientation,
                    aligned_fname):
        factory = rad_factory.ReadAlignmentDataFactory(
            self.TN_PARAMS, end, orientation)
        rads_by_id = factory.DictFromFiles(self.MASKED_READS_FNAME,
                                           aligned_fname)
        for read_id, rad in rads_by_id.iteritems():
            self.assertEquals(self.TN_PARAMS.backbone_start_offset, rad.start_offset)
            self.assertEquals(end, rad.fixed_seq_end)
            self.assertEquals(orientation, rad.fixed_seq_orientation)
            self.assertEquals(read_id, rad.read_record.id)

        outf = StringIO()
        factory.WriteCSVFile(rads_by_id.itervalues(), outf)

    def testAligned5p(self):
        self.GenericTest('5p', self.FORWARD, self.ALIGNED_5P)

    def testAligned3p(self):
        self.GenericTest('3p', self.FORWARD, self.ALIGNED_3P)

    def testAligned5pReverse(self):
        self.GenericTest('5p', self.REVERSE, self.ALIGNED_5P_REV)

    def testAligned3pReverse(self):
        self.GenericTest('3p', self.REVERSE, self.ALIGNED_3P_REV)


if __name__ == '__main__':
    unittest.main()