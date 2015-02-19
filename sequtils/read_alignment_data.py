# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

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
        self._has_fixed = None
        self._fixed_seq = None
        self._fixed_start = None
        self._fixed_end = None
        
        # If set, this is the position of the insertion.
        # The read may have no insertion in it.
        self._has_insert_backbone_matches = None
        self._has_insertion = None
        self._has_forward_insertion = None
        self._insertion_site = None
        self._insertion_in_frame = None
        self._linker_length = None
        self._insert_match_end = None
        self._insert_match_strand = None
        self._backbone_match_strand = None
        # Insertion site is the DNA coordinate in the backbone of
        # of the insertion. Calculated to give the same value regardless
        # of whether the read matched the 3' of 5' end of the insert.
        # Position is calculated to be the site of the 5' end of the 
        # insert in the backbone just prior to the fixed sequence introduced
        # by the transposon.
        
        # If set, the sequence of the linker between the fixed
        # sequence and the 
        self._linker_seq = None
    
    DICT_FIELDNAMES = ['read_id', 'insert_match_end', 'insert_match_strand',
                       'backbone_match_strand', 'fixed_seq_start', 'fixed_seq_end',
                       'insertion_site', 'linker_length', 'linker_seq', 'insertion_in_frame']
    def AsDict(self):
        return {'read_id': self.read_id,
                'insert_match_end': self.insert_match_end,
                'insert_match_strand': self.insert_match_strand,
                'backbone_match_strand': self.backbone_match_strand,
                'fixed_seq_start': self.fixed_start,
                'fixed_seq_end': self.fixed_end,
                'insertion_site': self.insertion_site,
                'linker_length': self.linker_length,
                'linker_seq': self.linker_seq,
                'insertion_in_frame': self.insertion_in_frame}
    
    def PrettyPrint(self):
        read_seq = str(self.read_record.seq)
        if self._has_insertion is None:
            self.CalculateInsertion()
            
        insert_start = self.insert_hsp.query_start
        insert_end = self.insert_hsp.query_end
        backbone_start = self.backbone_hsp.query_start
        backbone_end = self.backbone_hsp.query_end
        annotation_l = ['-'] * len(read_seq)
        for i in xrange(insert_start, insert_end):
            if annotation_l[i] != '-':
                annotation_l[i] = 'X'
            else: 
                annotation_l[i] = 'I'
        for i in xrange(backbone_start, backbone_end):
            if annotation_l[i] != '-':
                annotation_l[i] = 'X'
            else:
                annotation_l[i] = 'B'
        fs = self._fixed_start or 0
        fe = self._fixed_end or 0
        for i in xrange(fs, fe):
            if annotation_l[i] != '-':
                annotation_l[i] = 'X'
            else:
                annotation_l[i] = 'F'
        
        position_l = ['-'] * len(read_seq)
        if self.insert_match_strand < 0:
            position_l[insert_start] = self.insert_hsp.hit_end
            position_l[insert_end - 1] = self.insert_hsp.hit_start
        else:
            position_l[insert_start] = self.insert_hsp.hit_start
            position_l[insert_end - 1] = self.insert_hsp.hit_end
        if self.backbone_match_strand < 0:
            position_l[backbone_start] = self.backbone_hsp.hit_end
            position_l[backbone_end - 1] = self.backbone_hsp.hit_start
        else:
            position_l[backbone_start] = self.backbone_hsp.hit_start
            position_l[backbone_end - 1] = self.backbone_hsp.hit_end
        print read_seq
        print ''.join(annotation_l)
        print ''.join(map(str, position_l))
        
        i_hsp = self.insert_hsp
        b_hsp = self.backbone_hsp
        print 'Insert query position %d:%d' % (i_hsp.query_start, i_hsp.query_end)
        print 'Backbone query position %d:%d' % (b_hsp.query_start, b_hsp.query_end)
        print 'Insert position %d:%d' % (i_hsp.hit_start, i_hsp.hit_end)
        print 'Backbone position %d:%d' % (b_hsp.hit_start, b_hsp.hit_end) 
        print 'Read length', len(self._read_record.seq)
        print 'Insert end matched', self.insert_match_end
        print 'Insert matches strand', self.insert_match_strand
        print 'Backbone matches strand', self.backbone_match_strand
        print 'Calculated insert position', self.insertion_site
        print 'Calculated linker length', self.linker_length
        

    @staticmethod
    def DictFromFiles(insert_alignment_fname,
                      backbone_alignment_fname,
                      reads_fasta_fname):
        """Makes a dictionary of initialized ReadAlignmentData.
        
        Args:
            insert_alignment_fname: the path to the BLAT output for the insert.
            backbone_alignment_fname: the path to the BLAT output for the backbone.
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
                # Don't even care if the match isn't contiguous.
                # TODO: test if this removed real matches.
                if hsp.is_fragmented:
                    continue
                
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
                # Don't even care if the match isn't contiguous.
                # TODO: test if this removed real matches.
                if hsp.is_fragmented:
                    continue
                
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

    def ConsistentInsertAndBackboneMatches(self):
        return (self.HasBothHSPs() and
                self._backbone_hsp.query_id == self._insert_hsp.query_id and 
                self._n_insert_match_fragments == 1 and
                self._n_backbone_match_fragments == 1 and
                self._insert_match_strand == self._backbone_match_strand)
        
    def CalculateInsertion(self):
        self._has_insert_backbone_matches = self.HasBothHSPs()
        self._has_insertion = self.HasBothHSPs()
        self._has_forward_insertion = self.ConsistentInsertAndBackboneMatches()
        if not self._has_insertion:
            self._has_fixed = False
            return
        
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
        self._has_fixed = (fixed_pos != -1)
        if not self._has_fixed:
            self._has_insertion = False
            self._has_forward_insertion = False
            return
        
        self._fixed_start = fixed_pos
        self._fixed_end = fixed_pos + len(fixed_seq)
            
        b_hsp = self._backbone_hsp
        i_hsp = self._insert_hsp
        
        insert_position = None
        linker_length = None
        linker_seq = None
        same_strand = self._insert_match_strand == self._backbone_match_strand
        diff_sign = 1 if same_strand else -1
        if self._insert_match_end == '5p':
            if self._insert_match_strand > 0:
                # + strand of read for insert, matched 5' end of insert.
                bb_match_end = b_hsp.hit_end if same_strand else b_hsp.hit_start
                bb_end_inq = b_hsp.query_end
                diff = (fixed_pos - bb_end_inq) #* diff_sign 
                insert_offset = 0 if same_strand else 5
                insert_position = bb_match_end + diff + insert_offset
                linker_length = i_hsp.query_start - self._fixed_end
                linker_seq = self.read_seq[self._fixed_end+1:i_hsp.query_start]
            else:
                # - strand of read, matched 5' end of insert
                bb_match_end = b_hsp.hit_end if same_strand else b_hsp.hit_start
                bb_start_inq = b_hsp.query_start
                diff = diff_sign * (bb_start_inq - (fixed_pos + len(fixed_seq)))
                insert_offset = 0 if same_strand else 5
                insert_position = bb_match_end + diff + insert_offset
                linker_length = self._fixed_start - i_hsp.query_end
                linker_seq = self.read_seq[i_hsp.query_end:self._fixed_start-1]
                linker_seq = linker_seq.reverse_complement()
        else:
            if self._insert_match_strand > 0:
                # + strand of read, matched 3' end of insert.
                bb_match_start = b_hsp.hit_start if same_strand else b_hsp.hit_end
                bb_start_inq = b_hsp.query_start
                diff = (bb_start_inq - (fixed_pos + len(fixed_seq))) #* diff_sign 
                insert_offset = 5 if same_strand else 0
                insert_position = bb_match_start - diff + insert_offset
                linker_length = self._fixed_start - i_hsp.query_end
                linker_seq = self.read_seq[i_hsp.query_end:self._fixed_start]
            else:
                # - strand of read, matched 3' end of insert.
                bb_match_start = b_hsp.hit_start if same_strand else b_hsp.hit_end
                bb_end_inq = b_hsp.query_end
                diff = diff_sign * (fixed_pos - bb_end_inq)
                insert_offset = 5 if same_strand else 0
                insert_position = bb_match_start - diff + insert_offset
                linker_length = i_hsp.query_start - self._fixed_end
                linker_seq = self.read_seq[self._fixed_end:i_hsp.query_start]
                linker_seq = linker_seq.reverse_complement()
        self._insertion_site = insert_position
        if self._insert_match_end == '5p':
            linker_length -= 1
        self._linker_length = linker_length
        self._linker_seq = linker_seq
        
        insert_in_frame = ((insert_position + 1) % 3 == 0)
        if not self._has_forward_insertion:
            insert_in_frame = False
        linker_in_frame = (linker_length % 3) == 0
        self._insertion_in_frame = (insert_in_frame and linker_in_frame)
        
        # NOTE: +5 bp on the 3' end of the insert b/c of the
        # 5 bp duplicated by the transposase
    
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
    has_insert_backbone_matches = property(lambda self: len(self._calc_if_none('_has_insert_backbone_matches')))
    has_forward_insertion = property(lambda self: self._calc_if_none('_has_forward_insertion'))
    has_insertion = property(lambda self: self._calc_if_none('_has_insertion'))
    insertion_site = property(lambda self: self._calc_if_none('_insertion_site'))
    insertion_in_frame = property(lambda self: self._calc_if_none('_insertion_in_frame'))
    linker_length = property(lambda self: self._calc_if_none('_linker_length'))
    insert_match_end = property(lambda self: self._insert_match_end)
    insert_match_strand = property(lambda self: self._insert_match_strand)
    backbone_match_strand = property(lambda self: self._backbone_match_strand)
    fixed_start = property(lambda self: self._calc_if_none('_fixed_start'))
    fixed_end = property(lambda self: self._calc_if_none('_fixed_end'))
    fixed_seq = property(lambda self: self._calc_if_none('_fixed_seq'))
    
    
    
    