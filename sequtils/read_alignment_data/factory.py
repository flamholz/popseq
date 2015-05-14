# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

import csv
import pysam
import sequtils.read_alignment_data as rad

from Bio import SeqIO


class ReadAlignmentDataFactory(object):
    
    ALLOWED_FIXED_SEQ_ENDS = ('3p', '5p')
    
    def __init__(self, tn_params,
                 fixed_seq_end, fixed_seq_orientation):
        """Object that makes ReadAlignmentData instances.
        
        Args:
            tn_params: TranspositionParams object.
            fixed_seq_end: '3p' or '5p' - whether this sequence represents
                the 3' or 5' end of the transposon insert.
            orientation: +/-1. +1 indicating fixed sequence found as
                forward complement, -1 for reverse.
        """
        assert fixed_seq_end in self.ALLOWED_FIXED_SEQ_ENDS
        assert abs(fixed_seq_orientation) == 1
        self.tn_params = tn_params 
        self.fixed_seq_end = fixed_seq_end
        self.fixed_seq_orientation = fixed_seq_orientation
    
    def NewBuilder(self, read_record):
        """Creates a ReadAlignmentBuilder object."""
        return rad.ReadAlignmentData.Builder(
            self.tn_params, self.fixed_seq_end,
            self.fixed_seq_orientation,
            read_record=read_record)

    def DictFromFiles(self,
                      masked_reads_fname,
                      aligned_reads_fname):
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
        
        print masked_reads_fname
        print aligned_reads_fname
        masked_reads_reader = SeqIO.parse(masked_reads_fname, 'fastq')
        for i, record in enumerate(masked_reads_reader):
            if i % 10000 == 0:
                print 'Processing masked read %d' % i
            
            read_id = record.id
            if read_id in read_data_builders:
                print 'Already saw ID', read_id
            builder = self.NewBuilder(record)
            read_data_builders[read_id] = builder
            
        alignedf = pysam.AlignmentFile(aligned_reads_fname, 'rb')
        for i, it in enumerate(alignedf.fetch()):
            # Note that the aligner includes the read description in the name.
            # They are separated by whitespace.
            read_id = it.query_name.split()[0]
            if read_id not in read_data_builders:
                print 'Never saw read ID', read_id
                continue
            
            builder = read_data_builders[read_id]
            builder.backbone_alignment = it

        ids_to_remove = []
        for read_id, builder in read_data_builders.iteritems():
            if not builder.Complete():
                ids_to_remove.append(read_id)
        for read_id in ids_to_remove:
            # Can't modify dict while iterating over it.
            read_data_builders.pop(read_id)

        # Generate ReadAlignmentData from builder for every remaining read.
        read_data = dict((k, v.Build()) for k,v in read_data_builders.iteritems())
        return read_data
    
    def DictFromFileLists(self,
                          masked_reads_fnames,
                          aligned_reads_fnames):
        """Generate a dictionary of ReadAlignmentData from lists of filenames.
        
        Lists assumed to be in the same order. Iterates over tuples.
        
        Args:
            insert_alignment_fnames: list of paths to BLAT output for the insert.
            backbone_alignment_fnames: list of paths to BLAT output for the backbone.
            reads_fasta_fnames: list of paths to reads in FASTA format.
        """
        read_data_by_id = {}
        zipped = zip(masked_reads_fnames, aligned_reads_fnames)
        for masked_fname, aligned_fname in zipped:
            read_data_by_id.update(self.DictFromFiles(masked_fname,
                                                      aligned_fname))
        return read_data_by_id
    
    def WriteCSVFile(self, iterable, outfile):
        """Writes ReadAlignmentData to a CSV file.
        
        Args:
            iterable: iterable of ReadAlignmentData objects.
            outfile: file like object.
        """
        w = csv.DictWriter(outfile, rad.ReadAlignmentData.DICT_FIELDNAMES)
        w.writeheader()
        for rd in iterable:
            w.writerow(rd.AsDict())
                
    def WriteCSVFilename(self, iterable, out_fname):
        """Writes ReadAlignmentData to a CSV file.
        
        Args:
            iterable: iterable of ReadAlignmentData objects.
            out_fname: filename to write to.
        """
        with open(out_fname, 'w') as outf:
            self.WriteCSVFile(iterable, outf)