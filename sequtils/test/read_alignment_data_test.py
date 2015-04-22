#!/usr/bin/python

import unittest
import glob
import numpy as np

from sequtils import read_alignment_data as rad
from sequtils.read_alignment_data import factory


class ReadAligmentDataTest(unittest.TestCase):

    INSERT_ALIGNMENT_PSL = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_insert_aligned.pslx'    
    BACKBONE_ALIGNMENT_PSL = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k_backbone_aligned.pslx'
    READS_FASTA = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped100k.fa'
    
    FACTORY = factory.ReadAlignmentDataFactory(backbone_start_offset=23,
                                               fixed_5p_seq=rad.DEFAULT_FIXED_5P_SEQ,
                                               fixed_3p_seq=rad.DEFAULT_FIXED_3P_SEQ)
    READ_DATA = {}
    
    @classmethod
    def setUpClass(cls):
        cls.READ_DATA = cls.FACTORY.DictFromFiles(cls.INSERT_ALIGNMENT_PSL,
                                                  cls.BACKBONE_ALIGNMENT_PSL,
                                                  cls.READS_FASTA)
        print 'Finished reading test data.'
    
    def test5pPlusStrandPerfectBboneAlignment(self):
        print 'test5pPlusStrandPerfectBboneAlignment'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:14592:10486'
        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('CCTTCT', rd.linker_seq.tostring())
        self.assertEquals(693, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals('5p', rd._insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
    
    def test5pPlusStrandAlignmentOverrun(self):
        print 'test5pPlusStrandAlignmentOverrun'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:16999:7093'
        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(253, rd.insertion_site)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('TCTCCT', rd.linker_seq.tostring())
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
        
    def test5pMinusStrand(self):
        print 'test5pMinusStrand'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:7152:10611'
        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(84, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(3, rd.linker_length)
        self.assertEquals('GCT', rd.linker_seq.tostring())
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
    
    def test3pPlusStrand(self):
        print 'test3pPlusStrand'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:13100:2818'
        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(3516, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('TCTGCT', rd.linker_seq.tostring())
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
        
    def test3pMinusStrand(self):
        print 'test3pMinusStrand'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:9758:3066'
        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(24, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('CCTCCT', rd.linker_seq.tostring())
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
        
    def test5BPOverrun1(self):
        print "test5BPOverrun1"
        ex_key = "HWI-700460R:435:C5EWEACXX:4:1101:16901:9478"

        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(352, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('CCTTCT', rd.linker_seq.tostring())
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
    
    def testSuperWeirdRead(self):
        ex_key = "HWI-700460R:435:C5EWEACXX:4:1101:11878:10018"

        rd = self.READ_DATA[ex_key]
        
        self.assertFalse(rd.has_insertion)
        self.assertIsNone(rd.insertion_site)
        self.assertIsNone(rd.linker_length)
        self.assertIsNone(rd.linker_seq)
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
    
    def testGeneratedReads(self):
        """Reads generated from real sequencing data."""
        print 'testGeneratedReads'
        reads = sorted(glob.glob('sequtils/test/data/cas9_pdz_inserts/Cas9_*.fa'))
        insert_aligned = sorted(glob.glob('sequtils/test/data/cas9_pdz_inserts/Cas9_*insert_aligned.pslx'))
        backbone_aligned = sorted(glob.glob('sequtils/test/data/cas9_pdz_inserts/Cas9_*backbone_aligned.pslx'))
        
        for insert_fname, bbone_fname, reads_fname in zip(insert_aligned, backbone_aligned, reads):
            read_data = self.FACTORY.DictFromFiles(insert_fname, bbone_fname, reads_fname)
            self.assertEquals(4, len(read_data), msg='should have 4 reads per example that matched')
            
            insert_sites = []
            for rd in read_data.itervalues():
                self.assertTrue(rd.has_insertion)
                # All these examples are in frame.
                self.assertTrue(rd.insertion_in_frame)
                self.assertIsNotNone(rd.linker_seq)
                self.assertIsNotNone(rd.linker_length)
                self.assertIsNotNone(rd.fixed_start)
                self.assertIsNotNone(rd.fixed_start)
                self.assertIsNotNone(rd.fixed_seq)
                
                # Check the insertion site.
                expected_insert_site = int(rd.read_record.description.split('_')[-3])
                self.assertEquals(expected_insert_site, rd.insertion_site)
                
                insert_sites.append(rd.insertion_site)
            
            insert_sites = np.array(insert_sites)
            self.assertTrue(np.all(insert_sites[0] == insert_sites)) 
        
    def testReverseInsertion5pInsertPlusBackboneMinus(self):
        print "testReverseInsertion5pInsertPlusBackboneMinus"
        ex_key = "HWI-700460R:435:C5EWEACXX:4:1101:15409:7214"

        rd = self.READ_DATA[ex_key]
        rd.PrettyPrint()
        
        self.assertTrue(rd.has_insertion)
        self.assertFalse(rd.has_forward_insertion)
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
        self.assertEquals(3477, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('GCTGCT', rd.linker_seq.tostring())

    def testReverseInsertion3pInsertPlusBackboneMinus(self):
        print "testReverseInsertion3pInsertPlusBackboneMinus"
        ex_key = "HWI-700460R:435:C5EWEACXX:4:1101:7921:3089"

        rd = self.READ_DATA[ex_key]
        
        self.assertTrue(rd.has_insertion)
        self.assertFalse(rd.has_forward_insertion)
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
        self.assertEquals(1962, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length) 
        self.assertEquals('TCTCCT', rd.linker_seq.tostring())

    def testReverseInsertion3pInsertMinusBackbonePlus(self):
        print "testReverseInsertion3pInsertMinusBackbonePlus"
        ex_key = "HWI-700460R:435:C5EWEACXX:4:1101:1835:4416"

        rd = self.READ_DATA[ex_key]
        
        self.assertTrue(rd.has_insertion)
        self.assertFalse(rd.has_forward_insertion)
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
        self.assertEquals(300, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('CCTCCT', rd.linker_seq.tostring())

    def testReverseInsertion5pInsertMinusBackbonePlus(self):
        print "testReverseInsertion5pInsertMinusBackbonePlus"
        ex_key = "HWI-700460R:435:C5EWEACXX:4:1101:15525:6413"

        rd = self.READ_DATA[ex_key]
        
        self.assertTrue(rd.has_insertion)
        self.assertFalse(rd.has_forward_insertion)
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
        self.assertEquals(3224, rd.insertion_site)
        self.assertFalse(rd.insertion_in_frame)
        self.assertEquals(6, rd.linker_length)
        self.assertEquals('GCTGCT', rd.linker_seq.tostring())
        
    def testRunAll(self):
        """Rename to run."""
        for key, read_data in self.READ_DATA.iteritems():
            if read_data.has_insertion:
                read_data.insertion_site
                if not read_data.has_forward_insertion:
                    print 
                    print key
                    read_data.PrettyPrint()
                    print


if __name__ == '__main__':
    unittest.main()