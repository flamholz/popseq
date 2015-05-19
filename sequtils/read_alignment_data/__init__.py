# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

from Bio import Seq


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
        def __init__(self, transposition_params,
                     fixed_seq_end, fixed_seq_orientation,
                     read_record=None, backbone_alignment=None):
            # These values may be set after construction.
            self.read_record = read_record
            self.backbone_alignment = backbone_alignment
        
            self.tn_params = transposition_params
            self.fixed_seq_end = fixed_seq_end
            self.fixed_seq_orientation = fixed_seq_orientation

        def Complete(self):
            """Has all required data."""
            return (self.read_record and 
                    self.backbone_alignment and 
                    self.tn_params and 
                    self.fixed_seq_end and
                    self.fixed_seq_orientation)
        
        def Build(self):
            """Build a ReadAlignmentData object.
            
            Returns an instance of ReadAlignmentData from accumulated data.
            """
            assert self.Complete()
            rad = ReadAlignmentData(
                self.read_record,
                self.backbone_alignment,
                self.tn_params,
                self.fixed_seq_end,
                self.fixed_seq_orientation)
            return rad
    
    def __init__(self,
                 read_record,
                 backbone_alignment,
                 transposition_params,
                 fixed_seq_end, fixed_seq_orientation):
        """Object encapsulating variant calling from reads.
        
        TODO: rework all this. Very strange pattern we are using here. 
        
        Args:
            read_record: ID of the read.
            backbone_alignment: alignment of backbone. 
            backbone_start_offset: offset of the start codon into the
                sequence used to match the backbone (nt units). Used
                to make insertion sites relative to the start codon.
            fixed_seq: fixed sequence in this read.
            fixed_seq_orientation: the orientation of the fixed sequence (+/-1).
            fixed_seq_end: whether this is the 3' or 5' end.
        """
        self.tn_params = transposition_params
        self.backbone_alignment = backbone_alignment
        self.read_record = read_record
        self.fixed_seq_end = fixed_seq_end
        self.fixed_seq_orientation = fixed_seq_orientation
                
        self.insertion_index = None
        self.insertion_site = None
        self.linker_seq = None
        self.expected_insert_end_seq = None
        
        # True if the linker matches the configured pattern
        self.valid_linker = False  
        
        # Strands of read to which matches of insert/backbone were found
        self.insert_match_strand = None
        self.backbone_match_strand = None
        
        # Indices into the read of the part that matches the fixed sequence
        self.fixed_seq_start_idx = None
        self.fixed_seq_end_idx = None
        
        # Indices into the read of the part that matches the fixed sequence
        self.linker_start_idx = None
        self.linker_end_idx = None
        
        # Indices into the read of the part that aligned to backbone
        self.bbone_start_idx = None
        self.bbone_end_idx = None

        # Indices into the read of the part that aligned to insert
        self.insert_start_idx = None
        self.insert_end_idx = None
        
        # Whether the insertion is in-frame with the backbone sequence.
        # Note: this can be true even if the insert is reversed since
        # it will still allow for in-frame translation of the 3' end
        # of the backbone. 
        self.in_frame_insertion = False
        # True if the insert was inserted in the same orientation as
        # the backbone (5' - 3')
        self.forward_insertion = False
        
        
        # Calculate insertion parameters on construction!
        self._CalculateInsertion()
        
    """
    TODO: more complete CSV ouput
    DICT_FIELDNAMES = ['read_id', 'forward_insertion',
                       'insert_match_end', 'insert_match_strand',
                       'backbone_match_strand', 'fixed_seq_start', 'fixed_seq_end',
                       'insertion_site', 'linker_length', 'linker_seq', 'insertion_in_frame']
    """
    DICT_FIELDNAMES = ['read_id', 'insertion_index', 'insertion_site']
    
    def AsDict(self):
        return {'read_id': self.read_record.id,
                'insertion_index': self.insertion_index,
                'insertion_site': self.insertion_site}
        """
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
        """
    def PrettyPrint(self):
        """Prints the read data nicely to the commandline."""
        read_seq = str(self.read_record.seq)
        backbone_start = self.bbone_start_idx
        backbone_end = self.bbone_end_idx
        insert_start = self.insert_start_idx
        insert_end = self.insert_end_idx
        linker_start = self.linker_start_idx
        linker_end = self.linker_end_idx
        
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
        fs = self.fixed_seq_start_idx
        fe = self.fixed_seq_end_idx
        for i in xrange(fs, fe):
            if annotation_l[i] != '-':
                annotation_l[i] = 'X'
            else:
                annotation_l[i] = 'F'
        for i in xrange(linker_start, linker_end):
            if annotation_l[i] != '-':
                annotation_l[i] = 'X'
            else:
                annotation_l[i] = 'L'
        
        print '##### Read %s #####' % self.read_record.id 
        print ''.join(annotation_l)
        print self.read_record.seq
        print 'Expected insert seq', self.expected_insert_end_seq
        print 'Fixed sequence end', self.fixed_seq_end
        print 'Fixed sequence orientation', self.fixed_seq_orientation
        print 'Backbone matches strand', self.backbone_match_strand
        print 'Calculated insert start', self.insert_start_idx
        print 'Calculated insert end', self.insert_end_idx
        print 'Calculated insertion index', self.insertion_index
        print 'Calculated insertion site', self.insertion_site
        print 'Calculated linker seq', self.linker_seq
        print 'Valid linker seq', self.valid_linker
        print 'Calculated linker start', self.linker_start_idx
        print 'Calculated linker end', self.linker_end_idx
        print 'Forward insertion:', self.forward_insertion
        print 'In frame insertion:', self.in_frame_insertion
        print 
    
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
        # Shorthand
        bbone_al = self.backbone_alignment
        bp_dup = self.tn_params.tn_bp_duplicated
        start_offset = self.tn_params.backbone_start_offset
        
        self.backbone_match_strand = -1 if bbone_al.is_reverse else 1
        
        aligned_seq = Seq.Seq(bbone_al.query_alignment_sequence)
        if self.backbone_match_strand < 0:
            aligned_seq = aligned_seq.reverse_complement()
            
        self.bbone_start_idx = self.read_record.seq.find(aligned_seq)
        self.bbone_end_idx = self.bbone_start_idx + len(aligned_seq)
        
        # Position of insertion in nt sequence of target
        # not adjusted for start codon
        insertion_index = None
        if self.fixed_seq_end == '5p':
            if self.fixed_seq_orientation > 0:
                # Matched forward complement of 5' fixed sequence
                insertion_index = bbone_al.reference_end
                if self.backbone_alignment.is_reverse:
                    insertion_index = bbone_al.reference_start + bp_dup
            else:
                # Matched reverse complement of 5' fixed sequence
                insertion_index = bbone_al.reference_start + bp_dup
                if self.backbone_alignment.is_reverse:
                    insertion_index = bbone_al.reference_end
        elif self.fixed_seq_end == '3p':
            if self.fixed_seq_orientation > 0:
                # Matched forward complement of 5' fixed sequence
                insertion_index = bbone_al.reference_start + bp_dup
                if self.backbone_alignment.is_reverse:
                    insertion_index = bbone_al.reference_end
            else:
                # Matched reverse complement of 5' fixed sequence
                insertion_index = bbone_al.reference_end
                if self.backbone_alignment.is_reverse:
                    insertion_index = bbone_al.reference_start + bp_dup
        else:
            raise ValueError('Illegal value for fixed end %s' % self.fixed_end)
        
        self.insertion_index = insertion_index
        self.insertion_site = insertion_index - start_offset
        
        fixed_seq = self.tn_params.GetFixedSequence(self.fixed_seq_end,
                                                    self.fixed_seq_orientation)
        fixed_l = len(fixed_seq)
        fixed_seq_start = self.read_record.seq.find(fixed_seq)
        self.fixed_seq_start_idx = fixed_seq_start
        self.fixed_seq_end_idx = fixed_seq_start + fixed_l
        
        n_insert_bp = 10
        insert_end_seq = self.tn_params.GetInsertEndSequence(
            self.fixed_seq_end, self.fixed_seq_orientation, n_insert_bp)
        self.expected_insert_end_seq = insert_end_seq
        self.insert_start_idx = self.read_record.seq.find(insert_end_seq)
        self.insert_end_idx = self.insert_start_idx + n_insert_bp
        # TODO: handle the case that we don't find the insert...
        
        # Bases added on the 3' end of 5' linker to keep stuff in frame
        n_extra_bp_5p = self.tn_params.n_extra_bp_5p
        if self.fixed_seq_end == '5p':
            if self.fixed_seq_orientation > 0:
                fixed_end = fixed_seq_start + fixed_l
                linker_start = fixed_end + n_extra_bp_5p
                linker_end = self.insert_start_idx 
            else:
                fixed_end = fixed_seq_start - fixed_l
                linker_start = self.insert_end_idx
                linker_end = fixed_seq_start - n_extra_bp_5p
        elif self.fixed_seq_end == '3p':
            # Note: the offsetting nt is only on the 5' end.
            # Therefore we only need that offset correction above.
            # TODO: should make that a part of TranspositionParams.
            if self.fixed_seq_orientation > 0:
                fixed_end = fixed_seq_start - fixed_l
                linker_start = self.insert_end_idx
                linker_end = fixed_seq_start
            else:
                fixed_end = fixed_seq_start + fixed_l
                linker_start = fixed_end
                linker_end = self.insert_start_idx
        
        self.linker_start_idx = linker_start
        self.linker_end_idx = linker_end
        linker = self.read_record.seq[linker_start:linker_end]
        if self.fixed_seq_orientation < 0:
            linker = linker.reverse_complement()
        self.linker_seq = linker
        self.valid_linker = self.tn_params.ValidLinker(linker)
        
        self.in_frame_insertion = ((self.insertion_site + 1) % 3 == 0)
        self.forward_insertion = (self.backbone_match_strand == self.fixed_seq_orientation)