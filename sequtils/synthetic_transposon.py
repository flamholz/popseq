#!/usr/bin/python

import json
import random

from Bio.SeqRecord import SeqRecord


class Fragment(object):
    """Fragment of a larger sequence, usually representing a read."""

    INFO_DICT_KEYS = [
        'ins', 'inx', 'l5p', 'l3p', 'fst', 'fen', 'cn', 'rn', 'frd', 'fin'
        ]
    
    def __init__(self, id_num, trans, frag_start, frag_end, forward):
        self.id = id_num
        self.trans = trans # Transposition from which derived.
        self.start = frag_start
        self.end = frag_end
        
        self.seq = trans.construct[frag_start:frag_end] # Fragment sequence
        self.forward = forward
        if not self.forward:
            self.seq = self.seq.reverse_complement()
        
        inlen = trans.insert_len
        ins_start = trans.insertion_site + 5
        ins_end = trans.insertion_site + 5 + inlen
        extbp = 25
        efs = frag_start + extbp
        efe = frag_end - extbp
        
        should_match_5p = efs <= ins_start <= efe
        should_match_3p = efs  <= ins_end <= efe
        self.should_match = (should_match_3p or should_match_5p)
        id_tuple = (self.trans.id, self.id, self.trans.insertion_site)
        self.id_str = "C%06dR%04dINSPS%04d\n" % id_tuple
        
        # Longer names leads to a longer JSON string that gets mangled
        # by the aligner for some reason...
        self.info_dict = {
            # Site of insertion into coding sequence, NT units
            "ins": self.trans.expected_insertion_site,
            # Site of insertion into complete target sequence, i.e.
            # without accounting for ORF start. NT units.
            "inx": self.trans.insertion_site,
            # Linker on 5' side of insert.
            "l5p": str(self.trans.linker_5p),
            # Linker on 3' side of insert.
            "l3p": str(self.trans.linker_3p),
            # Fragment start in target coords, NT units/
            "fst": self.start,
            # Fragment end in target coords, NT units/
            "fen": self.end,
            # Construct ID number
            "cn": self.trans.id,
            # Read ID number
            "rn": self.id,
            # True if a forward read.
            "frd": self.forward,
            # True if a forward insertion.
            "fin": self.trans.forward_insertion}
        self.desc = json.dumps(self.info_dict)
        
    def ToSeqRecord(self):
        record = SeqRecord(self.seq, id=self.id_str,
                           name=self.id_str,
                           description=self.desc)
        return record
    
    @staticmethod
    def ParseInfoDict(read_record):
        desc = read_record.description
        dict_start = desc.find('{')
        desc = desc[dict_start:]
        return json.loads(desc)


class Transposition(object):
    """Simulates transposition into a defined backbone."""
    
    def __init__(self, id_num, insert_gen, target_seq, target_orf_start):
        """Init.
        
        Args:
            id_num: identifier of this transposition.
            insert_gen: generates the insert, potentially with random linkers.
            target_seq: sequence that is being transposed into.
            target_orf_start: start of the ORF in the target sequence.
        """
        self.id = id_num
        self.insert_gen = insert_gen
        self.insert, self.linker_5p, self.linker_3p = insert_gen.next()
        self.insert_len = len(self.insert)
        self.target = target_seq
        self.orf_start = target_orf_start
        self.forward_insertion = True
        
        # Random insert position within target sequence
        insert = self.insert
        target = self.target
        # Insert site relative to start of backbone sequence, not start codon.
        self.insertion_site = random.randrange(0, len(self.target)-5+1)
        
        # 50% of inserts are reverse
        if random.random() > 0.5:
            insert = insert.reverse_complement()
            self.forward_insertion = False
            
        # New sequence with insert tranposed and 5bp repeat created and ORF insert site calculated
        ins = self.insertion_site
        self.construct = target[:ins+5] + insert + target[ins:]
        # Expected insertion site 5' most base of insert sequence relative to start codon.
        self.expected_insertion_site = ins + 5 - (target_orf_start-1)
    
    def Shear(self, frag_id, fragment_length=100):
        """Generate a random subsequence of given length.
        
        Returns:
            Fragment object.
        """
        shear_start = random.randrange(0, len(self.construct)-fragment_length+1)
        shear_end = shear_start + fragment_length
        fwd = (random.random() <= 0.5)
        return Fragment(frag_id, self, shear_start, shear_end, fwd)