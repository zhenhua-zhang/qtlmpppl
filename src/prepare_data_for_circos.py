#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''A ustility script to preapre data for Circos tool.

Usage:
    prepare_data_for_circos.py [-h] -s [SNP-FILE [SNP-FILE ...]]
                                  [-S SNPS_FILE_SEP]
                                  [-c TRAIT-TYPE CHROM POS P-VAL]
                                  [-p MAX_PVAL]
                                  [--chrom-max-pval CHROM_MAX_PVAL]
                                  [-o OUTPUT_PREF] [-f {ssv,csv,tsv}]

Prepare data for circos from QTL mapping results by MatrixEQTL

Options:
  -h, --help            show this help message and exit
  -s [SNP-FILE [SNP-FILE ...]], --snps-files [SNP-FILE [SNP-FILE ...]]
                        Files including QTLs from MatrixEQTL
  -S SNPS_FILE_SEP, --snps-file-sep SNPS_FILE_SEP
                        The filed seperator in the SNPs files. Default: ,
  -c TRAIT-TYPE CHROM POS P-VAL, --kept-cols-names TRAIT-TYPE CHROM POS P-VAL
                        The columns for the Circos input file. Default:
                        ('gene', 'SequenceName', 'Position', 'pvalue')
  -p MAX_PVAL, --max-pval MAX_PVAL
                        The minimum p-value used. Default: 0.005
  --chrom-max-pval CHROM_MAX_PVAL
                        The maximal p-value for a chrom to pickup. Default:
                        5e-08
  -o OUTPUT_PREF, --output-prefix OUTPUT_PREF
                        The output file. Default: snps_density.per_cell_type
  -f {ssv,csv,tsv}, --output-fmt {ssv,csv,tsv}
                        The format of the output file. Default: ssv

Notes:
    TBA

Todo:
    TBA
'''

import argparse
import pandas as pd

from qmutils import CHROM_LEN_GRCH37


def getopt():
    '''Prepare data for Circos.
    '''
    parser = argparse.ArgumentParser(description='Prepare data for circos from QTL mapping results by MatrixEQTL')
    parser.add_argument('-s', '--snps-files', dest='snps_file_pl', metavar='SNP-FILE', required=True, nargs='*', help='Files including QTLs from MatrixEQTL')
    parser.add_argument('-S', '--snps-file-sep', dest='snps_file_sep', default=',', help='The filed seperator in the SNPs files. Default: %(default)s')
    parser.add_argument('-c', '--kept-cols-names', dest='kept_cols', metavar=('TRAIT-TYPE', 'CHROM', 'POS', 'P-VAL'), default=('gene', 'SequenceName', 'Position', 'pvalue'), nargs=4, help='The columns for the Circos input file. Default: %(default)s')
    parser.add_argument('-p', '--max-pval', dest='max_pval', default=5e-3, type=float, help='The minimum p-value used. Default: %(default)s')
    parser.add_argument('--chrom-max-pval', dest='chrom_max_pval', default=5e-8, type=float, help='The maximal p-value for a chrom to pickup. Default: %(default)s')
    parser.add_argument('-o', '--output-prefix', dest='output_pref', default='snps_density.per_cell_type', help='The output file. Default: %(default)s')
    parser.add_argument('-f', '--output-fmt', dest='output_fmt', default='ssv', choices=('ssv', 'csv', 'tsv'), help='The format of the output file. Default: %(default)s')

    return parser


def make_shift_dict(chosen_chr, chr_fmt_func=None):
    '''Make shift dict.
    '''
    if chr_fmt_func is None:
        chr_fmt_func = lambda y: ('chr{}'.format(x) for x in y)

    chr_pair = sorted(zip(chosen_chr, chr_fmt_func(chosen_chr)), key=lambda x: int(x[0]))

    pos = [CHROM_LEN_GRCH37[_chr] for _, _chr in chr_pair]
    shift_pos = [sum(pos[:i]) for i in range(len(pos))]

    return dict(zip(chosen_chr, shift_pos))


def update_record(row, shift_dict, index, trait_func=None):
    """Update coordnation by shift_dict.

    Note: The bug could be wrong order of row.
    """
    trait_type, chrom, pos, pval = row

    if trait_func:
        trait_type = trait_func(trait_type)

    pos += pos + shift_dict[chrom]
    return pd.Series([trait_type, pos, pos, pval], index=index)


def update_trait_name(trait):
    '''Update trait name.
    '''
    return trait.replace('Pgated_CCR5P_', '') \
            .replace('GM_CCR5P_', '') \
            .replace('_log10', '') \
            .replace('(', '') \
            .replace(')', '')


def main():
    '''Main entry of the script.
    '''
    output_header = ['#chrom', 'start', 'end', 'value']
    fmt_dict = {'ssv': ' ', 'csv': ',', 'tsv': '\t'}

    parser = getopt().parse_args()
    snps_file_pl = parser.snps_file_pl
    snps_file_sep = parser.snps_file_sep
    max_pval = parser.max_pval
    kept_cols = parser.kept_cols
    chrom_max_pval = parser.chrom_max_pval
    output_pref = parser.output_pref
    output_fmt = parser.output_fmt


    _trait_col, chr_col, _pos_col, pval_col = kept_cols

    qry_str = '{} < {}'.format(pval_col, max_pval)
    circos_dtfm = pd.concat([
        pd.read_csv(flpt, sep=snps_file_sep).query(qry_str).loc[:, kept_cols]
        for flpt in snps_file_pl
    ])

    qry_str = '{} < {}'.format(pval_col, chrom_max_pval)
    band_chr_pl = circos_dtfm.query(qry_str)[chr_col].drop_duplicates()
    base_shift = make_shift_dict(band_chr_pl)

    qry_str = '{} in {}'.format(chr_col, list(band_chr_pl))
    sep = fmt_dict.get(output_fmt, ' ')
    output_path = '.'.join([output_pref, output_fmt])
    circos_dtfm = circos_dtfm.query(qry_str) \
            .apply(update_record, axis=1, shift_dict=base_shift,
                   index=output_header, trait_func=update_trait_name) \
            .to_csv(output_path, sep=sep, index=False)


if __name__ == '__main__':
    main()
