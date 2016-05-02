dipseq
======

Code for handling data from domain insertion library sequencing.

# Data Formats
* Read data is expected in FASTQ format.
* Alignments are output in SAM and converted to BAM for analysis.
* Final output of variant calling is a CSV file.

# Python Dependencies
* Basic dependencies
 * scipy numpy pandas biopython
* Dependencies for enrichment analysis
 * rpy2

# External Dependencies
* Java runtime (for BBmap).
* C compiler (GCC, for fastq-grep installation).
* BBMap toolkit for fast filtering and alignment of reads. 
 * http://sourceforge.net/projects/bbmap/
* fastq-tools for regex trimming of reads.
 * You may need to install from source for most recent version.
 * https://github.com/dcjones/fastq-tools
* R, DESeq
 * http://bioconductor.org/packages/release/bioc/html/DESeq.html

# Installation
0. Make sure that you have GCC and the Java runtime installed.
 * If you are on a Macintosh, the easiest way to do this is to install the developer tools.
  * You can do this by running `xcode-select --install` from the command line.
 * You can also install the Java runtime independently by following these instructions.
  * https://www.java.com/en/download/help/mac_install.xml
1. Install the BBMap toolkit. 
 * Download the project zip from sourceforge (link above).
 * Unzip and move the unzipped folder somewhere permanent (e.g. make a ~/bin directory)
 * Add the BBMap folder to your PATH.
 * Test by running bbduk.sh from another folder.
2. Install the latest version of fastq-tools.
 * Clone the project into your local workspace (link above).  
 * Run `./autogen.sh` from inside the fastq-tools directory to generate the the makefile.
 * Run `./configure.sh && sudo make install` to configure, build and install the tools.
 * Note: if you don't have PCRE installed this will fail. It seems like re-installing the developer tools is the easiest solution to this problem.

If you follow the above steps, you have enough installed to run the insertion site identification script `analyze_insertions.py`. If this is all you want, then you are done. If you also want to use the DESeq analysis script `calculate_enrichment.py` then you also need to install the R runtime, the BioConductor and DESeq. You can follow instructions for installing R on their website (http://cran.r-project.org/bin/macosx/). You can install DESeq through BioConductor. BioConductor is available here: http://www.bioconductor.org/install/.  

# Running tests
* Install the nose test runner. 
 * sudo pip install nose
* Run all tests
 * nosetests ./
