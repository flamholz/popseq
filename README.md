popseq
======

Code for handling data from library sequencing.

# Data Formats
<<<<<<< HEAD
* Read data is expected in FASTQ format.
* Alignments are output in SAM and converted to BAM for analysis.
* Final output of variant calling is a CSV file.
=======
* Read data is expected in FASTQ or FASTA format.
* BLAT output is produced in PSLX format.
>>>>>>> 5b8ee12b612564e3cc56a1d82280a479a6f0342a

# Python Dependencies
* numpy
* scipy
* pandas
* biopython
* rpy2 (for R-based analysis)

# External Dependencies
<<<<<<< HEAD
* BBMap toolkit for fast filtering and alignment of reads.
* fastq-tools for regex trimming of reads.
* R, DESeq: http://bioconductor.org/packages/release/bioc/html/DESeq.html

# Running tests
* Install the nose test runner. 
** sudo pip install nose
* Run all tests
** nosetests ./
=======
* BLAT: https://genome.ucsc.edu/FAQ/FAQblat.html
* FASTX toolkit: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
* R, DESeq: http://bioconductor.org/packages/release/bioc/html/DESeq.html

>>>>>>> 5b8ee12b612564e3cc56a1d82280a479a6f0342a
