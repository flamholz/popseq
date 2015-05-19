popseq
======

Code for handling data from library sequencing.

# Data Formats
* Read data is expected in FASTQ format.
* Alignments are output in SAM and converted to BAM for analysis.
* Final output of variant calling is a CSV file.

# Python Dependencies
* numpy
* scipy
* pandas
* biopython
* rpy2 (for R-based analysis)

# External Dependencies
* BBMap toolkit for fast filtering and alignment of reads. 
 * http://sourceforge.net/projects/bbmap/
* fastq-tools for regex trimming of reads.
 * note: you need to install our custom version of fastq-tools.
 * https://github.com/SavageLab/fastq-tools
* R, DESeq
 * http://bioconductor.org/packages/release/bioc/html/DESeq.html

# Running tests
* Install the nose test runner. 
 * code(sudo pip install nose)
* Run all tests
 * code(nosetests ./)
