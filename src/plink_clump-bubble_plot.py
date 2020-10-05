#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''A script to plot heatmap'''

import argparse

import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap


def getopts():
    '''Get command line interface options.'''
    parser = argparse.ArgumentParser(description='A script to plot heatmap of SNPs')
    parser.add_argument('-i', '--from-file', dest='meta_file', default=None, help='A file from which read file path.')
    parser.add_argument('-F', '--output-fmt', dest='output_fmt', default=['png'], nargs='*', choices=['png', 'pdf', 'svg'], help='The output format. Default: %(default)s')
    parser.add_argument('-o', '--output-pref', dest='output_pref', default='./', help='The output prefix. Default: %(default)s')

    return parser.parse_args()


def load_meta_file(meta_file):
    '''Load file that including meta information for pre-/post-clump files.
    '''
    file_pool = []
    with open(meta_file, 'r') as ipf:
        for line in ipf:
            if line.startswith('#'):
                header = line.lstrip('#').strip('\n').split(',')
            else:
                record = line.strip('\n').split(',')
                file_pool.append(dict(zip(header, record)))

    return file_pool


def read_plink_clumped_file(fp):
    records = []
    with open(fp, 'r') as pst_clmp_fh:
        for line in pst_clmp_fh:
            if "CHR" in line:
                header = line.split()
            else:
                record = line.split()
                if record:
                    records.append(dict(zip(header, record)))

    return pd.DataFrame.from_records(records)


def load_file(meta_info, **kwargs):
    cohort = meta_info.get('cohort')
    mthd = meta_info.get('mthd')
    celltype = meta_info.get('celltype')
    pre_clmp_fp = meta_info.get('pre_clump')
    pst_clmp_fp = meta_info.get('post_clump')

    pre_clmp_df = pd.read_csv(pre_clmp_fp, **kwargs)
    pre_clmp_df['cohort'] = cohort
    pre_clmp_df['mthd'] = mthd
    pre_clmp_df['celltype'] = celltype

    pst_clmp_df = read_plink_clumped_file(pst_clmp_fp)
    #pre_clmp_df.loc[:, 'chrom'] = pst_clmp_df.loc[0, 'CHR']
    idx_snps = pst_clmp_df['SNP'].to_list()
    pre_clmp_df.loc[:, 'Index SNPs'] = 0
    pre_clmp_df.loc[pre_clmp_df['snps'].isin(idx_snps), 'Index SNPs'] = 'Y'

    return pre_clmp_df


def load_files(meta_info_pool, **kwargs):
    return pd.concat([load_file(x, **kwargs) for x in meta_info_pool])


def prepare_data(dtfm):
    dtfm.loc[:, 'celltype'] = dtfm['celltype'].apply(lambda x: x.split('_', 2)[-1])
    dtfm.loc[:, 'cohort'] = dtfm['cohort'].str.upper()
    dtfm.loc[:, 'mthd'] = dtfm['mthd'].replace({'gm': 'MFI', 'pc': 'PCP'})
    dtfm.loc[:, 'run'] = dtfm[['celltype', 'mthd', 'cohort']].apply(', '.join, axis=1)
    dtfm.loc[:, 'pvalue'] = dtfm['pvalue'].apply(lambda x: -np.log10(x))
    dtfm.loc[dtfm['pvalue'] >= 16, 'pvalue'] = 16
    dtfm.loc[dtfm['beta'] <= -2, 'beta'] = -2
    dtfm.loc[dtfm['beta'] >= 2, 'beta'] = 2
    idx_snps = dtfm.loc[dtfm['Index SNPs'] == 'Y', 'snps'].drop_duplicates().to_list()
    dtfm.loc[((dtfm['snps'].isin(idx_snps)) & (dtfm['Index SNPs'] == 0)), 'Index SNPs'] = 'N'
    #dtfm.loc[:, 'snps'] = (dtfm.loc[:, ['snps', 'chrom']]
    #.apply(lambda x: '{}, Chr{}'.format(*x.to_list()), axis=1))

    # dtfm.sort_values(['mthd', 'cohort', 'celltype'], inplace=True)
    return dtfm


def bubble_plot(dtfm):
    mpl.rcParams['legend.fontsize'] = 'large'
    fig, axes = plt.subplots()

    beta_max, beta_min = dtfm['beta'].max(), dtfm['beta'].min()
    colors = ["lightseagreen", 'white', "darkorange"]
    nodes = ([beta_min, 0, beta_max] - beta_min) / (beta_max - beta_min)
    cmap = LinearSegmentedColormap.from_list("mycap", list(zip(nodes, colors)))

    pmax, pmin = dtfm['pvalue'].max(), dtfm['pvalue'].min()

    arguments = dict(x='run', y='snps',
                     size='pvalue', sizes=(50, 250), size_norm=(pmin, pmax),
                     hue='beta', palette=cmap, hue_norm=(beta_min, beta_max),
                     style='Index SNPs', markers=['o', '*'], linewidth=.2,
                     edgecolor='0', axes=axes)

    idx_snps = dtfm['Index SNPs'] != 0
    sb.scatterplot(data=dtfm.loc[idx_snps, :], **arguments)

    axes.set_xlabel('Experiments')
    axes.set_ylabel('SNPs')

    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)

    axes.legend(bbox_to_anchor=(1, 1))
    plt.setp(axes.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    axes.grid(linestyle='--', linewidth=0.5, zorder=0)
    axes.set_axisbelow(True)

    fig.set_figwidth(12)
    fig.set_figheight(8)
    fig.set_tight_layout(True)

    return fig


def main():
    '''The main entry.'''
    opts = getopts()
    meta_file = opts.meta_file
    output_fmt = opts.output_fmt
    output_pref = opts.output_pref

    input_file_pool = load_meta_file(meta_file)
    dtfm = load_files(input_file_pool, sep='\t')
    dtfm = prepare_data(dtfm)
    fig = bubble_plot(dtfm)

    if isinstance(output_fmt, str):
        output_fmt = [output_fmt]

    flnm = 'bubble_plot' if output_pref.endswith('/') else '-bubble_plot'
    for fmt in output_fmt:
        optp = '{}{}.{}'.format(output_pref, flnm, fmt)
        fig.savefig(optp)


if __name__ == '__main__':
    main()
