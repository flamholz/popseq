#!/usr/bin/python

"""Merges counts per insertion site.

Ignores any per-site linker variability.
"""

import argparse
import numpy as np
import pandas as pd

from os import path
from scripts.util.filename_util import ForceExpand


def Main():
    parser = argparse.ArgumentParser(description='Merges counts per insertion site.')
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
        adf_fwd = adf[adf.forward_insertion == True]
        adf_rev = adf[adf.forward_insertion == False]
        
        myhist_fwd = np.histogram(adf_fwd.insertion_site, bins=np.arange(min_site, max_site))
        myhist_rev = np.histogram(adf_rev.insertion_site, bins=np.arange(min_site, max_site))
        counts_fwd, sites_fwd = myhist_fwd
        counts_rev, sites_rev = myhist_rev
        
        index_fwd = ['%d_fwd' % i for i in sites_fwd[:-1]]
        index_rev = ['%d_rev' % i for i in sites_rev[:-1]]
        
        series = pd.Series(counts_fwd, index_fwd)
        series = series.append(pd.Series(counts_rev, index_rev))
        counts_df[df_label] = series
        
    actual_site = np.array([int(ins.split('_')[0]) for ins in counts_df.index])
    fwd = [(ins.split('_')[1] == 'fwd') for ins in counts_df.index]
    in_frame = (actual_site + 1) % 3 == 0
    counts_df['insertion_site'] = actual_site
    counts_df['forward_insertion'] = fwd
    counts_df['in_frame'] = in_frame
    counts_df.index.name = 'insertion_site_name'
        
    print 'Writing merged counts to', args.output_fname
    counts_df.to_csv(args.output_fname)
    
    
if __name__ == '__main__':
    Main()