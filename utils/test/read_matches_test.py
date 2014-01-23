#!/usr/bin/python

import unittest

from Bio import SearchIO
from utils.read_matches import ReadMatches


class ReadMatchesTest(unittest.TestCase):

    def setUp(self):
        psl_fname = 'utils/test/sample.psl'
        parsed = SearchIO.parse(psl_fname, 'blat-psl')
        self.records = [r for r in parsed]
        self.first_record = self.records[0]
        self.sample_hsp = self.first_record.hsps[0]

    def testInit(self):
        read_id = 'foobar'
        matches = ReadMatches(read_id)

        # Prior to setting the matches, object is incomplete.
        self.assertFalse(matches.Complete())
        self.assertFalse(matches.Consistent())
        self.assertEquals(read_id, matches.read_id)
        self.assertIsNone(matches.insert_match)
        self.assertIsNone(matches.insert_location_in_query)
        self.assertIsNone(matches.insert_match_range)
        self.assertIsNone(matches.insert_match_strand)
        self.assertIsNone(matches.backbone_match)
        self.assertIsNone(matches.backbone_location_in_query)
        self.assertIsNone(matches.match_location_in_backbone)
        self.assertIsNone(matches.backbone_match_range)
        self.assertIsNone(matches.backbone_match_strand)

        # Still incomplete - only the insert match set.
        matches.insert_match = self.sample_hsp
        self.assertFalse(matches.Complete())
        self.assertFalse(matches.Consistent())
        self.assertIsNotNone(matches.insert_match)
        self.assertIsNotNone(matches.insert_location_in_query)
        self.assertIsNotNone(matches.insert_match_range)
        self.assertIsNotNone(matches.insert_match_strand)
        self.assertIsNone(matches.backbone_match)
        self.assertIsNone(matches.backbone_location_in_query)
        self.assertIsNone(matches.match_location_in_backbone)
        self.assertIsNone(matches.backbone_match_range)
        self.assertIsNone(matches.backbone_match_strand)

        # Now complete - both matches.
        matches.backbone_match = self.sample_hsp
        self.assertTrue(matches.Complete())
		# I happen to know that this sample match is inconsistent.
        self.assertFalse(matches.Consistent())
        self.assertIsNotNone(matches.insert_match)
        self.assertIsNotNone(matches.insert_location_in_query)
        self.assertIsNotNone(matches.insert_match_range)
        self.assertIsNotNone(matches.insert_match_strand)
        self.assertIsNotNone(matches.backbone_match)
        self.assertIsNotNone(matches.backbone_location_in_query)
        self.assertIsNotNone(matches.match_location_in_backbone)
        self.assertIsNotNone(matches.backbone_match_range)
        self.assertIsNotNone(matches.backbone_match_strand)


if __name__ == '__main__':
    unittest.main()
        
