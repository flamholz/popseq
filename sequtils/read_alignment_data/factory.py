# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

from Bio import SearchIO
from Bio import Seq
from Bio import SeqIO

import sequtils.read_alignment_data as rad


class ReadAlignmentDataFactory(object):
    def __init__(self, backbone_start_offset,
                 fixed_5p_seq, fixed_3p_seq):
        """Object that makes ReadAlignmentData instances.
        
        Args:
            backbone_start_offset: offset of the start codon into the
                sequence used to match the backbone (nt units). Used
                to make insertion sites relative to the start codon.
            fixed_5p_seq: fixed sequence verifying 5' end.
            fixed_3p_seq: fixed sequence verifying 3' end.
        """
        self._start_offset = backbone_start_offset 
        self._fixed_5p_seq = fixed_5p_seq
        self._fixed_3p_seq = fixed_3p_seq
    
    def NewBuilder(self, read_id):
        """Creates a ReadAlignmentBuilder object."""
        return rad.ReadAlignmentData.Builder(read_id, self._start_offset,
                                             self._fixed_5p_seq, self._fixed_3p_seq)

    def DictFromFiles(self,
                      insert_alignment_fname,
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
        read_data_builders = {}
        
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
                if read_id in read_data_builders:
                    duplicated_insert_ids.add(read_id)
                else:
                    # Create a new one if not there yet.
                    rb = self.NewBuilder(read_id)
                    rb.insert_hsp = hsp
                    read_data_builders[read_id] = rb
                
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
                rb = read_data_builders.get(read_id)
                if rb is not None:
                    # If there was no previous backbone match or this one is 
                    # longer, replace it. 
                    if (rb.backbone_hsp is None or
                        rb.backbone_hsp.query_span < hsp.query_span):
                        rb.backbone_hsp = hsp
        
        print len(duplicated_insert_ids), 'duplicated insert ids'
        print len(duplicated_backbone_ids), 'duplicated backbone ids'
        rids_to_remove = [rid for rid in read_data_builders
                          if not read_data_builders[rid].HasInsertAndBackboneHSPs()]
        print 'Removing', len(rids_to_remove), 'reads with no insertions'
        for rid in rids_to_remove:
            read_data_builders.pop(rid)
        print len(read_data_builders), 'with matches to insert and backbone'
        
        # Fetch the sequences of the reads.
        with open(reads_fasta_fname) as f:
            reader = SeqIO.parse(f, 'fasta')
            for record in reader:
                read_id = record.id
                if read_id in read_data_builders:
                    rb = read_data_builders[read_id]
                    rb.read_record = record

        # Generate ReadAlignmentData from builder for every remaining read.
        read_data = dict((k, v.Build()) for k,v in read_data_builders.iteritems())
        return read_data
    
    def DictFromFileLists(self,
                          insert_alignment_fnames,
                          backbone_alignment_fnames,
                          reads_fasta_fnames):
        """Generate a dictionary of ReadAlignmentData from lists of filenames.
        
        Lists assumed to be in the same order. Iterates over tuples.
        
        Args:
            insert_alignment_fnames: list of paths to BLAT output for the insert.
            backbone_alignment_fnames: list of paths to BLAT output for the backbone.
            reads_fasta_fnames: list of paths to reads in FASTA format.
        """
        read_data_by_id = {}
        zipped = zip(insert_alignment_fnames, backbone_alignment_fnames, reads_fasta_fnames)
        for insert_fname, backbone_fname, fasta_fname in zipped:
            read_data_by_id.update(self.DictFromFiles(insert_fname,
                                                      backbone_fname,
                                                      fasta_fname))
        return read_data_by_id
    