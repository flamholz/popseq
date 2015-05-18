#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
from sequtils.ambiguous_seq import AmbiguousSequence


class TranspositionParams(object):
    """Parameters associated with a transposition library."""
    
    def __init__(self, insert_seq, backbone_seq, backbone_start_offset,
                 fixed_seq_5p, fixed_seq_3p, linker_pattern,
                 max_linker_repeats, tn_bp_duplicated):
        """Initialize.
        
        Args:
            insert_seq: Bio.Seq insert sequence not including linkers
                or fixed sequences.
            backbone_seq: Bio.Seq backbone sequence. Might include some
                sequence flanking the ORF.
            backbone_start_offset: index of first base of start codon
                in backbone_seq.
            fixed_seq_5p: the fixed sequence found on the 5' end of
                the insert.
            fixed_seq_3p: the fixed sequence found on the 3' end of
                the insert.
            linker_pattern: AmbiguousSequence object defining the pattern.
            max_linker_repeats: the maximum number of repeats of the
                linker pattern.
            tn_bp_duplicated: the number of bases duplicated during
                transposition.
        """
        assert insert_seq
        assert backbone_seq
        assert fixed_seq_5p
        assert fixed_seq_3p
        
        self.insert_seq = insert_seq
        self.backbone_seq = backbone_seq
        self.backbone_start_offset = backbone_start_offset
        self.fixed_seq_5p = fixed_seq_5p
        self.fixed_seq_3p = fixed_seq_3p
        self.linker_pattern = linker_pattern
        self.max_linker_repeats = max_linker_repeats
        self.tn_bp_duplicated = tn_bp_duplicated
        
        self.fixed_seq_by_end = {('3p', 1): self.fixed_seq_3p,
                                 ('5p', 1): self.fixed_seq_5p,
                                 ('3p', -1): self.fixed_seq_3p.reverse_complement(),
                                 ('5p', -1): self.fixed_seq_5p.reverse_complement()}
    
    def GetFixedSequence(self, fixed_end, orientation):
        return self.fixed_seq_by_end[(fixed_end, orientation)]
    
    def GetInsertEndSequence(self, fixed_end, orientation, n):
        seq = self.insert_seq[:n]
        if fixed_end == '3p':
            seq = self.insert_seq[-n:]
        if orientation == -1:
            seq = seq.reverse_complement()
        return seq
    
    def ValidLinker(self, linker_seq):
        """Returns True if this linker sequence is valid.
        
        i.e. does the linker match the known pattern and fall within the
        maximum number of repeats.
        
        Args:
            linker_seq: sequence observed for the linker.
        """
        n = self.linker_pattern.IsInstance(linker_seq)
        valid = n >= 0 and n <= self.max_linker_repeats
        return valid
            
    @classmethod
    def LoadFASTA(cls, fname):
        return SeqIO.parse(fname, 'fasta').next().seq
    
    @classmethod
    def FromArgs(cls, args):
        """Parse from commandline arguments in a standard form."""
        insert_seq = cls.LoadFASTA(args.insert_seq_filename)
        bbone_seq = cls.LoadFASTA(args.backbone_db_filename)
        linker_pattern = AmbiguousSequence(args.linker_pattern)
        return TranspositionParams(insert_seq, bbone_seq, args.start_offset,
                                   Seq(args.fixed_5p), Seq(args.fixed_3p),
                                   linker_pattern, args.max_linker_repeats,
                                   args.tn_bp_duplicated)
    
    def __str__(self):
        l = ['<TranspositionParams backbone_start_offset=%d' % self.backbone_start_offset,
             'fixed_seq_5p="%s"' % self.fixed_seq_5p,
             'fixed_seq_3p="%s"' % self.fixed_seq_3p,
             'tn_bp_duplicated=%d>' % self.tn_bp_duplicated]
        return ' '.join(l)
        