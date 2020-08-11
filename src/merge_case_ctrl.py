#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
(zcat bcg.ccr5_exp.pc/qtlmapping/per_trait/Pgated_CCR5P_CD45.tsv.gz | head -1 | cut -f1-4,6- -d',' | tr ',' '\t';
    for x in $(find ./hiv.ccr5_exp.* -name '*.5e-2.tsv'); do
        awk -F',' '{if($4<5e-8) {for(i=1; i<NF; i++) {if (i!=5) {printf $i"\t"} } print $NF}}' ${x};
    done | sort -k1,2g -t$'\t'
) > hiv.ccr5_exp.5e-8.tsv

"""

import argparse

import pandas as pd

def getargs():
    """Get command line interface arguments.
    """
    parser = argparse.ArgumentParser(description="Merge QTL mapping outputs by SNPs and traits")
    parser.add_argument('--case-fps', dest="case_fp_pl", nargs='*', required=True, help="Case file path.")
    parser.add_argument('--case-sep', dest='case_sep', type=str, default='\t', help='Splitter for case files')
    parser.add_argument('--case-name', dest='case_name', type=str, default='_case', help='The name used for the case.')

    parser.add_argument('--ctrl-fps', dest="ctrl_fp_pl", nargs='*', required=True, help="Control file path.")
    parser.add_argument('--ctrl-name', dest='ctrl_name', type=str, default='_ctrl', help='The name used for the ctrl.')
    parser.add_argument('--ctrl-sep', dest='ctrl_sep', type=str, default='\t', help='Splitter for ctrl files')

    parser.add_argument('--cols-jnon', dest='cols_jnon', default=['snps', 'gene'], nargs='*', help='Columns to join on')
    parser.add_argument('--cols-asis', dest='cols_asis', default=['SequenceName', 'Position', 'EffectAllele', 'AlternativeAllele'], nargs='*', help='The columns kept as them is.')

    parser.add_argument('--otpt-pf', dest='otpt_pf', default='./', help="Output file output prefix. Default: %(default)s")

    return parser


def read_file_bundle(csvfppl, **kwargs):
    """Read bundle files."""
    return pd.concat([pd.read_csv(csvfp, **kwargs) for csvfp in csvfppl],
                     ignore_index=True)

def partial_join(left, right, both_on, outer_asis, how='outer',
                 suffixes=('_l', '_r')):
    """Partially join two DataFrames.
    """
    if isinstance(both_on, str):
        both_on = [both_on]

    if set(both_on) & set(outer_asis):
        raise ValueError("The elements in both_on should not be in outer_asis")

    asis_cols = both_on + outer_asis
    left_asis, right_asis = left.loc[:, asis_cols], right.loc[:, asis_cols]
    asis_dtfm = pd.concat([left_asis, right_asis], ignore_index=True) \
            .drop_duplicates()

    comb_cols = [x for x in left.columns if x not in outer_asis]
    left_comb, right_comb = left.loc[:, comb_cols], right.loc[:, comb_cols]
    comb_dtfm = pd.merge(left_comb, right_comb, 'outer', both_on, suffixes=suffixes)

    joint_dtfm = pd.merge(asis_dtfm, comb_dtfm, how, both_on)

    return joint_dtfm


def main():
    """The main entry for the script.
    """
    args = getargs().parse_args()
    case_fp_pl = args.case_fp_pl
    ctrl_fp_pl = args.ctrl_fp_pl

    cols_jnon = args.cols_jnon
    ctrl_name = args.ctrl_name
    case_sep = args.case_sep

    cols_asis = args.cols_asis
    case_name = args.case_name
    ctrl_sep = args.ctrl_sep

    otpt_pf = args.otpt_pf
    otpt_path = otpt_pf + 'merged.tsv'

    casedf = read_file_bundle(case_fp_pl, sep=case_sep)
    ctrldf = read_file_bundle(ctrl_fp_pl, sep=ctrl_sep)

    mgdf = partial_join(casedf, ctrldf, cols_jnon, cols_asis, 'outer',
                        [case_name, ctrl_name])

    mgdf.sort_values(by='gene', inplace=True)
    mgdf.to_csv(otpt_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
