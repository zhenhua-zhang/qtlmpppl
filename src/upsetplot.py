#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''A script to draw upset plot.
'''

import argparse

import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt


def getargs():
    '''Get CLI arguments'''
    parser = argparse.ArgumentParser(description="A scripts to draw upsetplot.")
    parser.add_argument('-i', '--intput-file', dest='input_file', required=True,
                        help='The input file.')
    parser.add_argument('-s', '--intput-sep', dest='input_sep', default='\t',
                        help='The field splitter of input file. Default: %(default)s')

    parser.add_argument('--y-col-idx', dest='y_col_idx', default='SNPS',
                        help='Index of columns used for Y axis. Default: %(default)s') # SNPs
    parser.add_argument('--left-x-col-idx', dest='left_x_col_idx', default='HIV_MODEL',
                        help='Index of columns uses for X axis of left panel. default: %(default)s') # Public HIV Traits
    parser.add_argument('--right-x-col-idx', dest='right_x_col_idx', default='MAPPED_TRAIT',
                        help='Index of columns used for X axis of right panel. Default: %(default)s') # Project HIV traits
    parser.add_argument('--bubble-size-col-idx', dest='bubble_size_col_idx', default='HIV_PVAL',
                        help='Index of columns used for the size of bubbles in right panel. Default: %(default)s') #
    parser.add_argument('--bubble-color-col-idx', dest='bubble_color_col_idx', default='HIV_BETA',
                        help='Index of columns used for the color of bubbles in right panel. Default: %(default)s') #

    parser.add_argument('-S', '--fig-size', default=(12, 16), dest='figsize', nargs=2, help='The figsize. Default: %(default)s')

    parser.add_argument('-F', '--output-fmt', dest='output_fmt', default=['svg'], choices=('png', 'svg', 'pdf'), metavar="FMT", nargs='*',
                        help='No more than three. Default: %(default)s')
    parser.add_argument('-o', '--output-pref', dest='output_pref', default='./test',
                        help='Output pref. Default: %(default)s')

    return parser


def load_input_file(iptf, **kwargs):
    '''A very Jilei function to load input file.
    '''
    return pd.read_csv(iptf, **kwargs)


def mk_itsdf(data, x_col='MAPPED_TRAIT', y_col='SNPS', class_least_count=3,
             idx_order=None):
    '''Make data frame for the its_plot().
    '''
    intersected_col = (data.loc[:, [x_col, y_col]].drop_duplicates()
                       .loc[:, x_col].str.split(', ').explode())

    tmp_dict = {}
    for idx, val in zip(intersected_col.index, intersected_col.values):
        if val in tmp_dict:
            tmp_dict[val] += [idx]
        else:
            tmp_dict[val] = [idx]

    if idx_order is None:
        idx_order = data.loc[:, y_col].index

    exp_tmp_dict = {key: [1 if _idx in val else 0 for _idx in idx_order]
                    for key, val in tmp_dict.items()}

    dtfm = pd.DataFrame(exp_tmp_dict).stack().reset_index()
    dtfm.columns = [y_col, x_col, 'plot']
    dtfm = dtfm.loc[dtfm['plot'] == 1, ]

    transform_dict = dict(zip(data.index, data.loc[:, y_col]))
    dtfm.loc[:, y_col] = dtfm.loc[:, y_col].apply(lambda x: transform_dict[x])

    x_col_vals = dtfm.loc[:, x_col].value_counts()
    x_col_vals = x_col_vals[x_col_vals >= class_least_count].index
    dtfm = dtfm.query('{} in @x_col_vals'.format(x_col))

    return dtfm


def its_plot(axes, data, x_col, y_col, y_col_order: dict, **kwargs):
    '''Plot intersect plot.
    '''
    x_vec, y_vec = data.loc[:, [x_col, y_col]].T.values
    x_vec_len = {x: list(x_vec).count(x) for x in set(x_vec)}

    x_y_coord = sorted(zip(x_vec, y_vec), key=lambda i: x_vec_len[i[0]])
    x_vec = [i[0] for i in x_y_coord]
    y_vec = [y_col_order[i[1]] for i in x_y_coord]

    axes.scatter(x_vec, y_vec, **kwargs)

    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_visible(False)

    axes.invert_xaxis()
    axes.set_xlabel(x_col)
    axes.tick_params(axis='y', which='both', length=0)

    axes.grid(True, linewidth=0.25, linestyle='--')

    plt.setp(axes.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    return axes


def bbl_plot(axes, data, **kwargs):
    '''Plot bubble plot.
    '''
    axes = sb.scatterplot(data=data, ax=axes, **kwargs)

    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    plt.setp(axes.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    artist, labels = axes.get_legend_handles_labels()
    axes.legend(artist, labels, fontsize='large')

    return axes


def draw_upsetplot(data, y_col, left_x_col, right_x_col=None,
                   bcolor_col=None, bsize_col=None, figsize=None,
                   output_path='./test.png', colormap='bwr',
                   bsize_func=None, bcolor_func=None):
    '''Draw upsetplot.
    '''
    itsdf = None if right_x_col is None else mk_itsdf(data, right_x_col, y_col)

    if itsdf is None:
        (fig, l_axes), r_axes = plt.subplots(), None
    else:
        left_right_ratio = [len(data.loc[:, left_x_col].drop_duplicates()) * 1.5,
                            len(data.loc[:, right_x_col].drop_duplicates())]
        gridspec = {'width_ratios': left_right_ratio}
        fig, (l_axes, r_axes) = plt.subplots(ncols=2, sharey=True,
                                             gridspec_kw=gridspec)

    if figsize is None:
        figwidth = sum(left_right_ratio) / 2
        figheight = len(data.loc[:, y_col].drop_duplicates()) / 3
    else:
        figwidth, figheight = figsize

    fig.set_figwidth(figwidth)
    fig.set_figheight(figheight)
    fig.set_tight_layout(True)

    if bsize_func:
        data.loc[:, bsize_col] = data.loc[:, bsize_col].apply(bsize_func)

    if bcolor_func:
        data.loc[:, bcolor_col] = data.loc[:, bcolor_col].apply(bcolor_func)

    bubble_size_range = (figwidth * 10, figwidth * 20)
    bbl_ax = bbl_plot(l_axes, data, x=left_x_col, y=y_col, hue=bcolor_col,
                      size=bsize_col, sizes=bubble_size_range,
                      palette=colormap)
    bbl_ax.set(title='Enrichment plot')

    if l_axes is not None:
        y_col_order = {v:i for i, v in enumerate(data[y_col].drop_duplicates())}
        its_ax = its_plot(r_axes, itsdf, right_x_col, y_col, c='0.1',
                          y_col_order=y_col_order)
        its_ax.set(title='Intersection plot')

    if isinstance(output_path, str):
        output_path = [output_path]

    for optpath in output_path:
        fig.savefig(optpath)

    plt.close()


def main():
    '''The main entry of the script.
    '''
    args = getargs().parse_args()
    input_file = args.input_file
    output_pref = args.output_pref
    output_fmt = args.output_fmt

    y_col_idx = args.y_col_idx
    left_x_col_idx = args.left_x_col_idx
    right_x_col_idx = args.right_x_col_idx
    bubble_size_col_idx = args.bubble_size_col_idx
    bubble_color_col_idx = args.bubble_color_col_idx

    figsize = args.figsize

    output_path = ['{}-upsetplot.{}'.format(output_pref, fmt)
                   for fmt in output_fmt]

    dtfm = load_input_file(input_file, sep='\t')
    draw_upsetplot(dtfm,
                   output_path=output_path,
                   figsize=figsize,
                   y_col=y_col_idx,
                   right_x_col=right_x_col_idx,
                   left_x_col=left_x_col_idx,
                   bcolor_col=bubble_color_col_idx,
                   bsize_col=bubble_size_col_idx,
                   bsize_func=lambda x: -np.log10(x))


if __name__ == '__main__':
    main()
