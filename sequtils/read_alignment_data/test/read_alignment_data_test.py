#!/usr/bin/python

import gzip
import sequtils.read_alignment_data.factory as rad_factory
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from sequtils import read_alignment_data as rad
from sequtils.ambiguous_seq import AmbiguousSequence
from sequtils.synthetic_transposon import Fragment
from sequtils.transposition_params import TranspositionParams


class ReadAlignmentDataTest(unittest.TestCase):
    RAW_READS = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6.fa.gz'
    FILTERED_READS = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_insert_bbone_filtered.fq'
    MASKED_READS_FNAME = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_insert_masked.fq'
    ALIGNED_3P = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_3p_aligned.bam'
    ALIGNED_5P = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_5p_aligned.bam'
    ALIGNED_3P_REV = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_3p_rev_aligned.bam'
    ALIGNED_5P_REV = 'sequtils/test/data/generated_PDZ_insertions_linker_1e6_filtered_trimmed_5p_rev_aligned.bam'
    
    TN_PARAMS = TranspositionParams(
        insert_seq=TranspositionParams.LoadFASTA('data/sequences/cas9/pdz_insert.fa'),
        backbone_seq=TranspositionParams.LoadFASTA('data/sequences/cas9/dCas9_cloning_ends.fa'),
        backbone_start_offset=23,
        fixed_seq_5p=Seq("TGCATC"),
        fixed_seq_3p=Seq("GCGTCA"),
        linker_pattern=AmbiguousSequence(Seq("BCT")),
        extra_bp_5p=Seq('T'),
        max_linker_repeats=3,
        tn_bp_duplicated=5)
    FORWARD = 1
    REVERSE = -1

    INSERT_SEQ_FNAME = 'data/sequences/cas9/pdz_insert.fa'
    INSERT_SEQ = None

    @classmethod
    def setUpClass(cls):
        parsed = SeqIO.parse(cls.INSERT_SEQ_FNAME, 'fasta')
        cls.INSERT_SEQ = parsed.next()

    # TODO: test precision and recall from raw read baseline.
    def GetReadAlignementData(self, end, orientation, aligned_fname):
        factory = rad_factory.ReadAlignmentDataFactory(
            self.TN_PARAMS, end, orientation)
        rads_by_id = factory.DictFromFiles(self.FILTERED_READS,
                                           aligned_fname)
        return rads_by_id
    
    def GetAllReadData(self):
        all_data = [
            self.GetReadAlignementData('5p', self.FORWARD, self.ALIGNED_5P),
            self.GetReadAlignementData('3p', self.FORWARD, self.ALIGNED_3P),
            self.GetReadAlignementData('5p', self.REVERSE, self.ALIGNED_5P_REV),
            self.GetReadAlignementData('3p', self.REVERSE, self.ALIGNED_3P_REV)]
        
        out_dict = dict()
        for d in all_data:
            out_dict.update(d)
        return out_dict
    
    def testCoverage(self):
        all_recovered = self.GetAllReadData()
        insert_len = len(self.INSERT_SEQ.seq)
        
        raw_reads_gzf = gzip.GzipFile(self.RAW_READS)
        raw_reads = SeqIO.parse(raw_reads_gzf, 'fasta')
        masked_reads = SeqIO.parse(self.MASKED_READS_FNAME, 'fastq')
        n_guessed = 0
        n_guessed_correctly = 0
        n_guessed_incorrectly = 0
        n_should_have_guessed = 0
        n_total_raw_reads = 0
        
        positions_did_not_guess = []
        false_negative_ids = set()
        for read_record in raw_reads:
            read_id = read_record.id
            read_info = Fragment.ParseInfoDict(read_record)
            len_5p = len(read_info['l5p'])
            len_3p = len(read_info['l3p'])
            insert_idx = read_info['inx']
            frag_start = read_info['fst']
            frag_end = read_info['fen']
            backbone_bp_in_read_5p = insert_idx - frag_start
            total_insert_len = (insert_len + len_3p + len_5p +
                                len(self.TN_PARAMS.fixed_seq_5p) +
                                len(self.TN_PARAMS.fixed_seq_3p))
            backbone_bp_in_read_3p = frag_end - insert_idx - total_insert_len
            
            should_guess = ((20 <= backbone_bp_in_read_5p <= 80) or
                            (20 <= backbone_bp_in_read_3p <= 80))
            rad = all_recovered.get(read_id)
            true_site = read_info['ins']
            did_guess = rad is not None and rad.insertion_site is not None
            did_match = rad is not None and rad.insertion_site == true_site
            
            n_guessed += did_guess
            n_guessed_correctly += did_match
            n_guessed_incorrectly += should_guess and not did_match
            n_should_have_guessed += should_guess and not did_guess
            n_total_raw_reads += 1
            if should_guess and not did_match:
                false_negative_ids.add(read_id)
                positions_did_not_guess.append(true_site)
        raw_reads_gzf.close() ## close zipfile
        
        # Perfect precision - all guesses were correct.
        self.assertEquals(n_guessed, n_guessed_correctly)
        
        # Perfect precision --
        # All the cases where we guessed "incorrectly" were actually 
        # cases where we did not guess at all.
        self.assertEquals(n_should_have_guessed, n_guessed_incorrectly)

        false_negative_rate = (float(n_should_have_guessed) /
                               (float(n_should_have_guessed) + float(n_guessed_correctly)))
        self.assertLessEqual(false_negative_rate, 0.07)
        
        print 'Number of total reads', n_total_raw_reads
        print 'Number guessed', n_guessed
        print 'Number guessed correctly', n_guessed_correctly
        print 'Number guessed incorrectly', n_guessed_incorrectly
        print 'Number should have guessed', n_should_have_guessed
        print 'False negative rate', false_negative_rate
        
        masked_reads_by_id = dict((r.id, r) for r in masked_reads)
        masked_ids = set(masked_reads_by_id.iterkeys())
        intersect = false_negative_ids.intersection(masked_ids)
        print 'NUmber of false negative reads', len(false_negative_ids)
        print 'Intersection size with filtered reads', len(intersect)
        #for read_id in intersect:
        #    print masked_reads_by_id[read_id].seq

    def GenericTest(self, end, orientation, aligned_fname):
        rads_by_id = self.GetReadAlignementData(
            end, orientation, aligned_fname)
        for rad in rads_by_id.itervalues():
            read_info = Fragment.ParseInfoDict(rad.read_record)
            expected_insertion_site = read_info['ins']
            expected_forward_ins = read_info['fin']
            expected_forward_read = read_info['frd']
            expected_linker_seq = read_info['l5p']
            if rad.fixed_seq_end == '3p':
                expected_linker_seq = read_info['l3p']

            # NOTE: not testing insertion index because I carelessly defined
            # then differently between the read generator and the read analysis.            
            self.assertEquals(expected_insertion_site, rad.insertion_site)
            self.assertEquals(expected_forward_ins, rad.forward_insertion)
            self.assertTrue(rad.valid_linker)
            self.assertEquals(expected_linker_seq, str(rad.linker_seq))
            self.assertEquals(expected_forward_read, rad.backbone_match_strand > 0)
            
            rad_dict = rad.AsDict()
            self.assertSetEqual(set(rad.DICT_FIELDNAMES),
                                set(rad_dict.keys()))
            # Data should come right out of the object.
            for k,v in rad_dict.iteritems():
                self.assertEquals(getattr(rad, k), v)

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