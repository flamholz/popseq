#!/usr/bin/python

import unittest

from sequtils.read_alignment_data import ReadAlignmentData


class ReadAligmentDataTest(unittest.TestCase):

    INSERT_ALIGNMENT_PSL = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped_insert_aligned.pslx'    
    BACKBONE_ALIGNMENT_PSL = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped_backbone_aligned.pslx'
    READS_FASTA = 'sequtils/test/data/DFS001_2_index6_GCCAAT_L004_R1_001_clipped.fa'
    
    READ_DATA = {}
    
    @classmethod
    def setUpClass(cls):
        cls.READ_DATA = ReadAlignmentData.DictFromFiles(cls.INSERT_ALIGNMENT_PSL,
                                                        cls.BACKBONE_ALIGNMENT_PSL,
                                                        cls.READS_FASTA)
        print 'Finished reading test data.'
    
    def test5pPlusStrandPerfectBboneAlignment(self):
        print 'test5pPlusStrandPerfectBboneAlignment'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:14592:10486'
        rd = self.READ_DATA[ex_key]
        print rd.read_seq
        print rd.insertion_site
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(693, rd.insertion_site)
        self.assertEquals('5p', rd._insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
    
    def test5pPlusStrandAlignmentOverrun(self):
        print 'test5pPlusStrandAlignmentOverrun'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:16999:7093'
        rd = self.READ_DATA[ex_key]
        print rd.read_seq
        print rd.insertion_site
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(253, rd.insertion_site)
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
        
    def test5pMinusStrand(self):
        print 'test5pMinusStrand'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:7152:10611'
        rd = self.READ_DATA[ex_key]
        print rd.read_seq
        print rd.insertion_site
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(84, rd.insertion_site)
        self.assertEquals('5p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
    
    def test3pPlusStrand(self):
        print 'test3pPlusStrand'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:13100:2818'
        rd = self.READ_DATA[ex_key]
        print rd.read_seq
        print rd.insertion_site
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(3516, rd.insertion_site)
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(1, rd.insert_match_strand)
        self.assertEquals(1, rd.backbone_match_strand)
        
    def test3pMinusStrand(self):
        print 'test3pMinusStrand'
        ex_key = 'HWI-700460R:435:C5EWEACXX:4:1101:9758:3066'
        rd = self.READ_DATA[ex_key]
        print rd.read_seq
        print rd.insertion_site
        
        self.assertTrue(rd.has_insertion)
        self.assertEquals(24, rd.insertion_site)
        self.assertEquals('3p', rd.insert_match_end)
        self.assertEquals(-1, rd.insert_match_strand)
        self.assertEquals(-1, rd.backbone_match_strand)
    
    """
    def testRunAll(self):
        for key, read_data in self.READ_DATA.iteritems():
            if read_data.has_insertion:
                read_data.insertion_site
    """ 