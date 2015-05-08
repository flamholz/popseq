#!/usr/bin/python

import unittest
import json
import pandas as pd
import pylab

from Bio import SeqIO
from Bio.SeqIO import FastaIO
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
        #self.assertGreater(frac_expected, 0.6)
    
    def testCoverage(self):
        should_match_ids = set()
        did_match_ids = set()
        not_found_sites = []
        ofh = open('false_negative_reads.fa', 'w')
        d_by_construct = {}
        writer = FastaIO.FastaWriter(ofh)
        writer.write_header()
        with open(self.READS_FASTA) as f:
            reader = SeqIO.parse(f, 'fasta')
            for record in reader:
                read_info = self._parseReadInfo(record)
                cnum = read_info["construct_num"]
                d = d_by_construct.setdefault(cnum, read_info)
                found_match = (record.id in self.READ_DATA)
                should_match = read_info['should_match']
                true_pos = found_match and should_match
                false_pos = found_match and not should_match
                true_neg = not found_match and not should_match
                false_neg = not found_match and should_match
                d['true_pos'] = d.get('true_pos', 0) + true_pos
                d['false_pos'] = d.get('false_pos', 0) + false_pos
                d['true_neg'] = d.get('true_neg', 0) + true_neg
                d['false_neg'] = d.get('false_neg', 0) + false_neg
                
                if found_match:
                    did_match_ids.add(record.id)
                if read_info["should_match"]:
                    should_match_ids.add(record.id)
                    if not found_match:
                        not_found_sites.append(read_info["insertion_site"])
                        writer.write_record(record)
        writer.write_footer()
        ofh.close()
        
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
        #self.assertGreater(precision, 0.85)
        #self.assertGreater(recall, 0.85)
        
        df = pd.DataFrame.from_dict(d, orient='index')
        df.to_csv('match_stats.csv')
        
        pylab.figure()
        pylab.hist(not_found_sites, bins=200)
        pylab.show()

if __name__ == '__main__':
    unittest.main()