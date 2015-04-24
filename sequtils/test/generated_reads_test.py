#!/usr/bin/python

import unittest
import glob
import numpy as np
import json
import pylab

from Bio import SeqIO
from sequtils import read_alignment_data as rad
from sequtils.read_alignment_data import factory


class ReadAligmentDataTest(unittest.TestCase):                                         

    INSERT_ALIGNMENT_PSL = 'sequtils/test/data/generated_transposition_reads_insert_aligned.pslx'    
    BACKBONE_ALIGNMENT_PSL = 'sequtils/test/data/generated_transposition_reads_backbone_aligned.pslx'
    READS_FASTA = 'sequtils/test/data/generated_transposition_reads.fa'
    
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
    
    @staticmethod
    def _parseReadInfo(rr):
        """Parses JSON encoded info stashed in read description."""
        desc = rr.description
        dict_start = desc.find('{')
        desc = desc[dict_start:]
        return json.loads(desc)
    
    def testCorrectInsertionSite(self):
        """Tests that all the reads matched yield the correct insertion site."""
        total_matches = len(self.READ_DATA)
        matches_expected = 0
        for read_id, rd in self.READ_DATA.iteritems():
            read_info = self._parseReadInfo(rd.read_record)
            insertion_site = read_info["insertion_site"]
            self.assertEquals(insertion_site, rd.insertion_site)
            self.assertEquals(0, rd.linker_length, msg='generated reads have no linkers')
            self.assertTrue('read_id' in rd.AsDict())
            matches_expected += read_info["should_match"]
        frac_expected = float(matches_expected) / float(total_matches)
        print 'Expected %d of %d (%.2f%%)' % (matches_expected, total_matches, 100*frac_expected)
        self.assertGreater(frac_expected, 0.6)
    
    def testCoverage(self):
        should_match_ids = set()
        did_match_ids = set()
        with open(self.READS_FASTA) as f:
            reader = SeqIO.parse(f, 'fasta')
            for record in reader:
                read_info = self._parseReadInfo(record)
                found_match = (record.id in self.READ_DATA)
                if found_match:
                    did_match_ids.add(record.id)
                if read_info["should_match"]:
                    should_match_ids.add(record.id)
        
        should_did = should_match_ids.intersection(did_match_ids)
        n_did_match = len(did_match_ids)
        n_found = len(should_did)
        n_expected = len(should_match_ids)
        recall = float(n_found) / float(n_expected)
        precision = float(n_found) / float(n_did_match)
        print 'precision: %.2f' % precision
        print 'recall: %.2f' % recall
        print 'total expected: %d' % n_expected
        print 'total found: %d' % n_did_match
        print 'num expected found: %d' % n_found
        self.assertGreater(precision, 0.85)
        self.assertGreater(recall, 0.85)
    
    def testEndCorrelation(self):
        matches = np.zeros((4200, 2))
        for read_id, rd in self.READ_DATA.iteritems():
            idx = rd.insertion_site
            col = 0 if rd.insert_match_end == '5p' else 1
            matches[idx,col] += 1
        corr = np.corrcoef(matches[:,0], matches[:,1])[1,0]
        # Reads at ends should correlate well even with random noise and low total read count.
        self.assertGreater(corr, 0.8)    

if __name__ == '__main__':
    unittest.main()