# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

import numpy as np

from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO

"""
TODO: extract linker sequence and length.
TODO: make sure that reads on the 3' and 5' end of the same 
      insertion will yield the same insertion site.
"""


class ReadAlignmentData(object):
    
    # Make these configurable.
    FIXED_5P_SEQ = Seq.Seq('TGCATC')
    FIXED_3P_SEQ = Seq.Seq('GCGTCA')
    
    def __init__(self, read_id):
        self.read_id = read_id
        self._read_record = None
        self._backbone_hsp = None
        self._insert_hsp = None
        
        # 6 BP sequence fixed based on cloning.
        # Differs on 3' and 5' ends of insert.
        # Start and end positions are in query coordinates.
        self._fixed_seq = None
        self._fixed_start = None
        self._fixed_end = None
        
        # If set, this is the position of the insertion.
        # The read may have no insertion in it.
        self._has_insertion = None
        self._insertion_site = None
        self._insert_match_end = None
        self._insert_match_strand = None
        self._backbone_match_strand = None
        # Insertion site is the DNA coordinate in the backbone of
        # of the insertion. Calculated to give the same value regardless
        # of whether the read matched the 3' of 5' end of the insert.
        # Position is calculated to be the site of the 5' end onf the 
        # insert in the backbone just prior to the fixed sequence introduced
        # by the transposon.
        
        # If set, the sequence of the linker between the fixed
        # sequence and the 
        self._linker_seq = None
    
    DICT_FIELDNAMES = ['read_id', 'insert_match_end', 'insert_match_strand',
                       'backbone_match_strand', 'fixed_seq_start', 'fixed_seq_end',
                       'insertion_site']
    def AsDict(self):
        return {'read_id': self.read_id,
                'insert_match_end': self.insert_match_end,
                'insert_match_strand': self.insert_match_strand,
                'backbone_match_strand': self.backbone_match_strand,
                'fixed_seq_start': self.fixed_start,
                'fixed_seq_end': self.fixed_end,
                'insertion_site': self.insertion_site}
            

    @staticmethod
    def DictFromFiles(insert_alignment_fname,
                      backbone_alignment_fname,
                      reads_fasta_fname):
        """Makes a dictionarof initialized ReadAlignmentData.
        
        Args:
            backbone_alignment_fname: the path to the BLAT output for the backbone.
            insert_alignment_fname: the path to the BLAT output for the insert.
            reads_fasta_fname: the path to the reads in FASTA format.
        
        Returns:
            A dictionary mapping read ID to ReadAlignmentData.
            Retains data only for those reads mapping to both the insert and the backbone 
            according to the alignment data passed in.
        """
        read_data = {}
        
        # Fetch all insert alignments
        # Do the insert first because there are many fewer matches to the insert.
        # Can avoid memory blowup of making records for all the backbone matches this way.
        insert_aligned_reader = SearchIO.parse(insert_alignment_fname,
                                               'blat-psl', pslx=True)
        duplicated_insert_ids = set()
        for i, record in enumerate(insert_aligned_reader):
            if i % 10000 == 0:
                print 'processing insert match %d' % i
                
            for hsp in record.hsps:
                
                read_id = hsp.query_id
                if read_id in read_data:
                    duplicated_insert_ids.add(read_id)
                else:
                    # Create a new one if not there yet.
                    rd = ReadAlignmentData(read_id)
                    rd.insert_hsp = hsp
                    read_data[read_id] = rd
                
        # Fetch all the backbone alignments
        backbone_aligned_ids = set()
        duplicated_backbone_ids = set()
        backbone_aligned_reader = SearchIO.parse(backbone_alignment_fname,
                                                 'blat-psl', pslx=True) 
        for i, record in enumerate(backbone_aligned_reader):
            if i % 10000 == 0:
                print 'processing backbone match %d' % i
                
            for hsp in record.hsps:
            
                read_id = hsp.query_id
                if read_id in backbone_aligned_ids:
                    duplicated_backbone_ids.add(read_id)
                backbone_aligned_ids.add(read_id)

                # Only save if we had an insert match there.
                rd = read_data.get(read_id)
                if rd is not None:
                    # If there was no previous backbone match or this one is 
                    # longer, replace it. 
                    if (rd.backbone_hsp is None or
                        rd.backbone_hsp.query_span < hsp.query_span):
                        rd.backbone_hsp = hsp
        
        print len(duplicated_insert_ids), 'duplicated insert ids'
        print len(duplicated_backbone_ids), 'duplicated backbone ids'
        rids_to_remove = [rid for rid in read_data
                          if not read_data[rid].HasBothHSPs()]
        print 'Removing', len(rids_to_remove), 'reads with no insertions'
        for rid in rids_to_remove:
            read_data.pop(rid)
        print len(read_data), 'with matches to insert and backbone'
        
        # Fetch the sequences of the reads.
        with open(reads_fasta_fname) as f:
            reader = SeqIO.parse(f, 'fasta')
            for record in reader:
                read_id = record.id
                if read_id in read_data:
                    rd = read_data[read_id]
                    rd.read_record = record

        return read_data
    
    def HasBackboneHSP(self):
        return self._backbone_hsp is not None

    def HasInsertHSP(self):
        return self._insert_hsp is not None

    def HasBothHSPs(self):
        return self.HasBackboneHSP() and self.HasInsertHSP()

    def BackboneInsertMatchesOverlap(self):
        """Test the case that they overlap at exactly one bp.

        NOTE: this code appears to be wrong....
        """
        b_hsp = self._backbone_hsp
        i_hsp = self._insert_hsp
        diff_start_end = i_hsp.query_start - b_hsp.query_end
        diff_end_start = i_hsp.query_end - b_hsp.query_start
        return np.sign(diff_start_end) == np.sign(diff_end_start)

    def ConsistentInsertAndBackboneMatches(self):
        return (self.HasBothHSPs() and
                self._backbone_hsp.query_id == self._insert_hsp.query_id and 
                self._n_insert_match_fragments == 1 and
                self._n_backbone_match_fragments == 1 and
                self._insert_match_strand == self._backbone_match_strand)
        
    def BackboneFirstInReadCoords(self):
        """Returns True if the backbone matched before the insert in 5' to 3' direction of the read."""
        return self._backbone_hsp.query_start < self._insert_hsp.query_start

    def BackboneInsertDistanceInRead(self):
        if not self.ConsistentInsertAndBackboneMatches():
            return np.NaN
        b_hsp = self._backbone_hsp
        i_hsp = self._insert_hsp
        if self.BackboneFirstInReadCoords():
            return i_hsp.query_start - b_hsp.query_end
        return b_hsp.query_start - i_hsp.query_end

    def CalculateInsertion(self):
        self._has_insertion = self.ConsistentInsertAndBackboneMatches()
        
        fixed_5p_seq = self.FIXED_5P_SEQ
        fixed_3p_seq = self.FIXED_3P_SEQ
        read = self.read_seq

        fixed_seq = fixed_5p_seq
        if self._insert_match_end == '3p':
            fixed_seq = fixed_3p_seq
        self._fixed_seq = fixed_seq

        # fixed position will be the coordinated in the query of the 5'
        # end of the fixed sequence. 
        # TODO: deal with case of multiple matches to fixed seq.
        fixed_pos = -1
        if self._insert_match_strand > 0:
            # Matched in the fwd direction.
            # Only need to test the + strand for fixed seq.
            fixed_pos = read.find(fixed_seq)
        else:
            # Matched in the rvs direction
            # Only need to test the - strand for fixed seq.
            fixed_pos = read.find(fixed_seq.reverse_complement())
        # Note: start is always 5' end of the fixed sequence.
        self._fixed_start = fixed_pos
        self._fixed_end = fixed_pos + len(fixed_seq)
            
        b_hsp = self._backbone_hsp
        
        insert_position = None
        if self._insert_match_end == '5p':
            if self._insert_match_strand > 0:
                # + strand of read, matched 5' end of insert.
                bb_match_end = b_hsp.hit_end
                bb_end_inq = b_hsp.query_end
                diff = fixed_pos - bb_end_inq
                insert_position = bb_match_end + diff
            else:
                # - strand of read, matched 5' end of insert
                bb_match_end = b_hsp.hit_end
                bb_start_inq = b_hsp.query_start
                diff = bb_start_inq - (fixed_pos + len(fixed_seq))
                insert_position = bb_match_end + diff
        else:
            if self._insert_match_strand > 0:
                # + strand of read, matched 3' end of insert.
                bb_match_start = b_hsp.hit_start
                bb_start_inq = b_hsp.query_start
                diff = bb_start_inq - (fixed_pos + len(fixed_seq))
                insert_position = bb_match_start - diff + 5
            else:
                # - strand of read, matched 3' end of insert.
                bb_match_start = b_hsp.hit_start
                bb_end_inq = b_hsp.query_end
                diff = fixed_pos - bb_end_inq
                insert_position = bb_match_start - diff + 5 
        self._insertion_site = insert_position
        
        # NOTE: +5 bp on the 3' end of the insert b/c of the
        # 5 bp duplicated by the transposase

        """
        i_hsp = self._insert_hsp
        print self.read_id, read.tostring()
        print 'Insert position', insert_position
        print 'Insert end', self._insert_match_end
        print 'bbone hit pos', b_hsp.hit_start, b_hsp.hit_end
        print 'bbone hit in q', b_hsp.query_start, b_hsp.query_end
        print 'insert hit in q', i_hsp.query_start, i_hsp.query_end
        print 'fixed pos', fixed_pos, self._insert_match_strand
        """
        return insert_position
    
    def _get_read_record(self):
        return self._read_record
    def _set_read_record(self, read_record):
        self._read_record = read_record
        
    def _get_backbone_hsp(self):
        return self._backbone_hsp
    def _set_backbone_hsp(self, hsp):        
        self._backbone_hsp = hsp
        self._n_backbone_match_fragments = len(hsp.fragments)
        self._backbone_match_strand = hsp.fragments[0].query_strand
    
    def _get_insert_hsp(self):
        return self._insert_hsp
    def _set_insert_hsp(self, hsp):
        self._insert_hsp = hsp
        self._n_insert_match_fragments = len(hsp.fragments)
        self._insert_match_strand = hsp.fragments[0].query_strand
        self._insert_match_end = '3p' if hsp.hit_id.endswith('_3p') else '5p'
    
    def _calc_if_none(self, key):
        val = self.__getattribute__(key)
        if val is not None:
            return val
        self.CalculateInsertion()
        return self.__getattribute__(key)
    
    # Properties        
    backbone_hsp = property(_get_backbone_hsp, _set_backbone_hsp)
    insert_hsp = property(_get_insert_hsp, _set_insert_hsp)
    read_record = property(_get_read_record, _set_read_record)
    read_seq = property(lambda self: self._read_record.seq)
    linker_seq = property(lambda self: self._calc_if_none('_linker_seq'))
    linker_len = property(lambda self: len(self._calc_if_none('_linker_seq')))
    has_insertion = property(lambda self: self._calc_if_none('_has_insertion'))
    insertion_site = property(lambda self: self._calc_if_none('_insertion_site'))
    insert_match_end = property(lambda self: self._insert_match_end)
    insert_match_strand = property(lambda self: self._insert_match_strand)
    backbone_match_strand = property(lambda self: self._backbone_match_strand)
    fixed_start = property(lambda self: self._calc_if_none('_fixed_start'))
    fixed_end = property(lambda self: self._calc_if_none('_fixed_end'))
    fixed_seq = property(lambda self: self._calc_if_none('_fixed_seq'))
    
    
    
    