#!/usr/bin/python

"""Analysis of relative enrichment of a variant between samples."""

import argparse
import numpy as np
import os
import pylab
import pandas as pd
import pandas.rpy.common as com

from Bio import SeqIO
from os import path
from rpy2 import robjects
from rpy2.robjects.packages import importr


def Main():
    parser = argparse.ArgumentParser(description='Filter reads.')
    parser.add_argument("--data_pre_csv", required=True,
                        help="Path to CSV summarizing variants in reference sample.")
    parser.add_argument("--data_post_csv", required=True,
                        help="Path to CSV summarizing variants in experimental sample.")
    parser.add_argument("-o", "--output_fname", required=True,
                        help="Where to write CSV output.")
    args = parser.parse_args()

    # Read in data as data frames
    pre_fname = args.data_pre_csv
    post_fname = args.data_post_csv
    assert path.exists(pre_fname), '%s does not exist' % pre_fname
    assert path.exists(post_fname), '%s does not exist' % post_fname
    print 'Reading input data'
    df_pre = pd.read_csv(pre_fname)
    df_post = pd.read_csv(post_fname)

    # Grab condition names from the filenames.
    # TODO(flamholz): allow these to supplied via cmdline. 
    pre_name, _ = path.splitext(path.split(pre_fname)[-1])
    post_name, _ = path.splitext(path.split(post_fname)[-1])

    dfs = [df_pre, df_post]
    df_labels = [pre_name, post_name]

    # Stratify on the end matched.
    print 'Stratifying data based on end matched.'
    end_dfs = []
    end_labels = []
    for end in ['3p', '5p']:
        for adf, label in zip(dfs, df_labels):
            end_dfs.append(adf[adf.insert_match_end == end])
            end_labels.append('%s_matched_%s' % (label, end))
    
    min_site = np.min([adf.insertion_site.min() for adf in end_dfs]) - 1
    max_site = np.max([adf.insertion_site.max() for adf in end_dfs])
    
    # Merge stratified count data into 1 dataframe
    # A column per condition and a row for every insertion site.
    counts_df = pd.DataFrame()
    for adf, df_label in zip(end_dfs, end_labels):
        myhist = np.histogram(adf.insertion_site, bins=np.arange(min_site, max_site))
        counts, sites = myhist
        counts_df[df_label] = pd.Series(counts, sites[:-1])
    counts_df.index.name = 'insertion_site'

    # Make an R dataframe and import R functionality.
    print 'Running DESeq enrichment analysis'
    r_df = com.convert_to_r_dataframe(counts_df)
    deseq = importr("DESeq")
    estimateSizeFactors = robjects.r["estimateSizeFactors"]
    sizeFactors = robjects.r["sizeFactors"]
    estimateDispersions = robjects.r["estimateDispersions"]
    fitInfo = robjects.r["fitInfo"]
    nbinomTest = robjects.r["nbinomTest"]
    plotDispEsts = robjects.r["plotDispEsts"]

    condition = robjects.FactorVector([pre_name, post_name] * 2)
    cds = deseq.newCountDataSet(r_df, condition)
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    deseq_res = nbinomTest(cds, pre_name, post_name)

    # Write output to CSV.
    print 'Writing DESeq results to', args.output_fname
    res_df = com.convert_robj(deseq_res)
    res_df.to_csv(args.output_fname)
    
    

if __name__ == '__main__':
    Main()