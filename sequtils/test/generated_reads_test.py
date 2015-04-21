#!/usr/bin/python

import unittest
import glob
import numpy as np
import json

from sequtils import read_alignment_data as rad

class ReadAligmentDataTest(unittest.TestCase):                                         

    INSERT_ALIGNMENT_PSL = 'sequtils/test/data/generated_transposition_reads_insert_aligned.pslx'    
    BACKBONE_ALIGNMENT_PSL = 'sequtils/test/data/generated_transposition_reads_backbone_aligned.pslx'
    READS_FASTA = 'sequtils/test/data/generated_transposition_reads.fa'
    
    FACTORY = rad.ReadAlignmentDataFactory(backbone_start_offset=23,
                                           fixed_5p_seq=rad.DEFAULT_FIXED_5P_SEQ,
                                           fixed_3p_seq=rad.DEFAULT_FIXED_3P_SEQ)
    READ_DATA = {}
    
    @classmethod
    def setUpClass(cls):
        cls.READ_DATA = cls.FACTORY.DictFromFiles(cls.INSERT_ALIGNMENT_PSL,
                                                  cls.BACKBONE_ALIGNMENT_PSL,
                                                  cls.READS_FASTA)
        print 'Finished reading test data.'
        
    def testCorrectInsertionSite(self):
        for read_id, rd in self.READ_DATA.iteritems():
            desc = rd.read_record.description
            dict_start = desc.find('{')
            desc = desc[dict_start:]
            
            read_info = json.loads(desc)
            insertion_site = read_info["insertion_site"]
            self.assertEquals(insertion_site, rd.insertion_site)
            self.assertTrue('read_id' in rd.AsDict())
            

if __name__ == '__main__':
    unittest.main()