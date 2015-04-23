# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

from Bio import Seq

DEFAULT_FIXED_5P_SEQ = Seq.Seq('TGCATC')
DEFAULT_FIXED_3P_SEQ = Seq.Seq('GCGTCA')


class ReadAlignmentData(object):
    """Contains alignments of a transposon mutagenesis library.
    
    Class performs all DNA math required to locate the site at which
    a transposon inserted a DNA sequence (the insert) into a backbone
    sequence. This particular class is designed to do the math for the 
    case where the inserts 5' and 3' ends are not reverse complements
    of each other. 
    
    Class assumes that the read matches both insert and backbone. This 
    will not be the case for all reads. Therefore, Builder subclass can
    be used to accumulate read-specific data before it is clear whether
    the read matched both insert and backbone (i.e. while reading files).
    
    TODO: rename this class appropriate to its function.
    """ 
    
    class Builder(object):
        """Builder for ReadAlignmentData.
        
        Used because we need a vessel to collect read-specific data before
        analysis. Allows construction of ReadAlignmentData objects only once
        all data available.
        """
        def __init__(self, read_id, backbone_start_offset,
                     fixed_5p_seq=DEFAULT_FIXED_5P_SEQ,
                     fixed_3p_seq=DEFAULT_FIXED_3P_SEQ):
            self.read_id = read_id
            self.read_record = None
            self.backbone_hsp = None
            self.insert_hsp = None
        
            self._start_offset = backbone_start_offset
            self._fixed_5p_seq = fixed_5p_seq
            self._fixed_3p_seq = fixed_3p_seq
        
        def HasInsertAndBackboneHSPs(self):
            return (self.insert_hsp is not None and
                    self.backbone_hsp is not None)
        
        def Complete(self):
            """Has all required data."""
            return (self.read_id and
                    self.read_record and 
                    self.insert_hsp and
                    self.backbone_hsp)
        
        def Build(self):
            """Build a ReadAlignmentData object.
            
            Returns an instance of ReadAlignmentData from accumulated data.
            """
            assert self.Complete()
            rd = ReadAlignmentData(self.read_id, self.read_record,
                                   self.insert_hsp, self.backbone_hsp, self._start_offset,
                                   self._fixed_5p_seq, self._fixed_3p_seq)
            return rd
    
    def __init__(self,
                 read_id,
                 read_record,
                 insert_hsp,
                 backbone_hsp,
                 backbone_start_offset,
                 fixed_5p_seq=DEFAULT_FIXED_5P_SEQ,
                 fixed_3p_seq=DEFAULT_FIXED_3P_SEQ):
        """Object encapsulating variant calling from reads.
        
        TODO: rework all this. Very strange pattern we are using here. 
        
        Args:
            read_id: ID of the read.
            backbone_start_offset: offset of the start codon into the
                sequence used to match the backbone (nt units). Used
                to make insertion sites relative to the start codon.
            fixed_5p_seq: fixed sequence verifying 5' end.
            fixed_3p_seq: fixed sequence verifying 3' end.
        """
        self.read_id = read_id
        self._fixed_5p_seq = fixed_5p_seq
        self._fixed_3p_seq = fixed_3p_seq
        self._read_record = read_record
        self._backbone_hsp = backbone_hsp
        self._insert_hsp = insert_hsp
        self._start_offset = backbone_start_offset
        
        # Fetch some data right out of the HSPs
        self._insert_match_strand = insert_hsp.fragments[0].query_strand
        self._insert_match_end = '3p' if insert_hsp.hit_id.endswith('_3p') else '5p'
        self._n_insert_match_fragments = len(insert_hsp.fragments)
        self._backbone_match_strand = backbone_hsp.fragments[0].query_strand
        self._n_backbone_match_fragments = len(backbone_hsp.fragments)
        
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
        
        # Insertion site is the DNA coordinate in the backbone of
        # of the insertion. Calculated to give the same value regardless
        # of whether the read matched the 3' of 5' end of the insert.
        # Position is calculated to be the site of the 5' end of the 
        # insert in the backbone just prior to the fixed sequence introduced
        # by the transposon.
        
        # If set, the sequence of the linker between the fixed
        # sequence and the 
        self._linker_seq = None
        
        # Calculate insertion parameters on construction!
        self._CalculateInsertion()
    
    DICT_FIELDNAMES = ['read_id', 'forward_insertion',
                       'insert_match_end', 'insert_match_strand',
                       'backbone_match_strand', 'fixed_seq_start', 'fixed_seq_end',
                       'insertion_site', 'linker_length', 'linker_seq', 'insertion_in_frame']
    def AsDict(self):
        return {'read_id': self.read_id,
                'forward_insertion': self.has_forward_insertion,
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
        """Prints the read data nicely to the commandline."""
        read_seq = str(self.read_record.seq)
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
    
    def HasBothHSPs(self):
        return (self._backbone_hsp is not None and
                self._insert_hsp is not None)

    def ConsistentInsertAndBackboneMatches(self):
        return (self.HasBothHSPs() and
                self._backbone_hsp.query_id == self._insert_hsp.query_id and 
                self._n_insert_match_fragments == 1 and
                self._n_backbone_match_fragments == 1 and
                self._insert_match_strand == self._backbone_match_strand)
        
    def _CalculateInsertion(self):
        """Calculates the position of the insertion in the backbone."""
        self._has_insert_backbone_matches = self.HasBothHSPs()
        self._has_insertion = self.HasBothHSPs()
        self._has_forward_insertion = self.ConsistentInsertAndBackboneMatches()
        if not self._has_insertion:
            self._has_fixed = False
            return
        
        fixed_5p_seq = self._fixed_5p_seq
        fixed_3p_seq = self._fixed_3p_seq
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
                diff = (fixed_pos - bb_end_inq) * diff_sign 
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
                diff = (bb_start_inq - (fixed_pos + len(fixed_seq))) * diff_sign 
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
        self._insertion_site = insert_position - self._start_offset
        if self._insert_match_end == '5p':
            linker_length -= 1
        self._linker_length = linker_length
        self._linker_seq = linker_seq
        
        insert_in_frame = ((self._insertion_site + 1) % 3 == 0)
        if not self._has_forward_insertion:
            insert_in_frame = False
        linker_in_frame = (linker_length % 3) == 0
        self._insertion_in_frame = (insert_in_frame and linker_in_frame)
        
        # NOTE: +5 bp on the 3' end of the insert b/c of the
        # 5 bp duplicated by the transposase
    
    @property
    def read_record(self):
        return self._read_record
    @property
    def backbone_hsp(self):
        return self._backbone_hsp
    @property
    def insert_hsp(self):
        return self._insert_hsp
    
    read_seq = property(lambda self: self._read_record.seq)
    linker_seq = property(lambda self: self._linker_seq)
    linker_len = property(lambda self: len(self._linker_seq))
    has_insert_backbone_matches = property(lambda self: self._has_insert_backbone_matches)
    has_forward_insertion = property(lambda self: self._has_forward_insertion)
    has_insertion = property(lambda self: self._has_insertion)
    insertion_site = property(lambda self: self._insertion_site)
    insertion_in_frame = property(lambda self: self._insertion_in_frame)
    linker_length = property(lambda self: self._linker_length)
    insert_match_end = property(lambda self: self._insert_match_end)
    insert_match_strand = property(lambda self: self._insert_match_strand)
    backbone_match_strand = property(lambda self: self._backbone_match_strand)
    fixed_start = property(lambda self: self._fixed_start)
    fixed_end = property(lambda self: self._fixed_end)
    fixed_seq = property(lambda self: self._fixed_seq)