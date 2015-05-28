#!/usr/bin/python

"""Analysis of enrichment of a variant between samples.

Uses DESeq R package to calculate normalized counts from replicated samples 
and, from the counts, calculate fold changes (post/pre) and attach P-values
to the fold change calculation.  
"""

import argparse
import pandas as pd
import pandas.rpy.common as com

from os import path
from rpy2 import robjects
from rpy2.robjects.packages import importr


def Main():
    parser = argparse.ArgumentParser(description='Calculate fold enrichment and P-values.')
    parser.add_argument("--data_pre_csv", required=True,
                        help="Path to CSV summarizing variants in reference sample.")
    parser.add_argument("--data_post_csv", required=True,
                        help="Path to CSV summarizing variants in experimental sample.")
    parser.add_argument("--pre_condition", required=True,
                        help="Name of the pre condition.")
    parser.add_argument("--post_condition", required=True,
                        help="Name of the post condition")
    parser.add_argument("-o", "--output_fname", required=True,
                        help="Where to write CSV output.")
    args = parser.parse_args()

    # Read in data as data frames
    pre_fname = args.data_pre_csv
    post_fname = args.data_post_csv
    assert path.exists(pre_fname), '%s does not exist' % pre_fname
    assert path.exists(post_fname), '%s does not exist' % post_fname
    print 'Reading input data'
    df_pre = pd.read_csv(pre_fname, index_col=0)
    df_post = pd.read_csv(post_fname, index_col=0)

    # Grab condition names from the filenames.
    pre_name = args.pre_condition
    post_name = args.post_condition

    # Merge the data frames
    merged_df = df_pre.reset_index().merge(df_post, how='outer').set_index('insertion_site_name')
    
    # Drop the non-relevant data
    cols_to_drop = ['insertion_site', 'forward_insertion', 'in_frame']
    merged_df_dropped = merged_df.drop(cols_to_drop, axis=1)
    merged_df_dropped = merged_df_dropped.fillna(0)
    merged_df_dropped = merged_df_dropped.astype(int)
    
    # Set up the condition labels
    condition_names = []
    for colname in merged_df_dropped.columns.tolist():
        if colname in df_pre.columns.tolist():
            condition_names.append(pre_name)
        elif colname in df_post.columns.tolist():
            condition_names.append(post_name)
        else:
            assert False, 'found unrecognized column name!'

    print 'Column names'
    print '\t%s' % ', '.join(merged_df_dropped.columns.tolist())
    print 'Condition names'
    print '\t%s' % ', '.join(condition_names)

    # Make an R dataframe and import R functionality.
    print 'Running DESeq enrichment analysis'
    r_df = com.convert_to_r_dataframe(merged_df_dropped)
    deseq = importr("DESeq")
    estimateSizeFactors = robjects.r["estimateSizeFactors"]
    sizeFactors = robjects.r["sizeFactors"]
    estimateDispersions = robjects.r["estimateDispersions"]
    fitInfo = robjects.r["fitInfo"]
    nbinomTest = robjects.r["nbinomTest"]
    plotDispEsts = robjects.r["plotDispEsts"]

    condition = robjects.FactorVector(condition_names)
    cds = deseq.newCountDataSet(r_df, condition)
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    deseq_res = nbinomTest(cds, pre_name, post_name)

    # Write output to CSV.
    print 'Writing DESeq results to', args.output_fname
    res_df = com.convert_robj(deseq_res)

    # Put the dropped columns back
    dropped_df = merged_df[cols_to_drop]
    res_df = res_df.merge(dropped_df, how='left', left_on='id', right_index=True)
    res_df.to_csv(args.output_fname)
    
    

if __name__ == '__main__':
    Main()