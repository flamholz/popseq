#!/usr/bin/python

"""Analysis of relative enrichment of a variant between samples."""

import argparse
import numpy as np
import pandas as pd

from os import path
from scripts.util.filename_util import ForceExpand


def Main():
    parser = argparse.ArgumentParser(description='Calculate fold enrichment and P-values.')
    parser.add_argument("-i", "--input_filenames", nargs='+', required=True,
                        help="Paths to CSV summarizing variants in reference sample.")
    parser.add_argument("-o", "--output_fname", required=True,
                        help="Where to write CSV output.")
    args = parser.parse_args()

    # Read in data as data frames
    all_fnames = ForceExpand(args.input_filenames)
    for fname in all_fnames:
        assert path.exists(fname)

    print 'Reading input data'
    dfs = [pd.read_csv(fname) for fname in all_fnames]
    df_labels = [path.splitext(path.split(fname)[-1])[0]
                 for fname in all_fnames]

    print 'Filtering reads where insert was not located'
    dfs = [df[df.insert_start_idx >= 0] for df in dfs]

    print 'Stratifying data based on end matched.'
    end_dfs = []
    end_labels = []
    for adf, label in zip(dfs, df_labels):
        for end in ['3p', '5p']:
            end_dfs.append(adf[adf.fixed_seq_end == end])
            end_labels.append('%s_matched_%s' % (label, end))
    
    print 'Calculating per-site counts'
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
    
    print 'Writing merged counts to', args.output_fname
    counts_df.to_csv(args.output_fname)
    
    
if __name__ == '__main__':
    Main()