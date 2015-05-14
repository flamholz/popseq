#!/usr/bin/python

from Bio import SeqIO


class TranspositionParams(object):
    """Parameters associated with a transposition library."""
    
    def __init__(self, insert_seq, backbone_seq, backbone_start_offset,
                 fixed_seq_5p, fixed_seq_3p, tn_bp_duplicated):
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
            tn_bp_duplicated: the number of bases duplicated during
                transposition.
        """
        self.insert_seq = insert_seq
        self.backbone_seq = backbone_seq
        self.backbone_start_offset = backbone_start_offset
        self.fixed_seq_5p = fixed_seq_5p
        self.fixed_seq_3p = fixed_seq_3p
        self.tn_bp_duplicated = tn_bp_duplicated
        
        self.fixed_seq_by_end = {'3p': self.fixed_seq_3p,
                                 '5p': self.fixed_seq_5p}
        
    @staticmethod
    def FromArgs(args):
        """Parse from commandline arguments in a standard form."""
        insert_seq = SeqIO.parse(args.insert_seq_filename, 'fasta').next()
        bbone_seq = SeqIO.parse(args.backbone_db_filename, 'fasta').next()
        return TranspositionParams(insert_seq, bbone_seq, args.start_offset,
                                   args.fixed_5p, args.fixed_3p, args.tn_bp_duplicated)
        
    def __str__(self):
        l = ['<TranspositionParams backbone_start_offset=%d' % self.backbone_start_offset,
             'fixed_seq_5p="%s"' % self.fixed_seq_5p,
             'fixed_seq_3p="%s"' % self.fixed_seq_3p,
             'tn_bp_duplicated=%d>' % self.tn_bp_duplicated]
        return ' '.join(l)
        