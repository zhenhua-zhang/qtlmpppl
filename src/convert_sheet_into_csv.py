#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

import pandas as pd


def get_args():
    """Parse CLI options.
    """
    parser = argparse.ArgumentParser(description="A simple script to convert Excel sheet into .csv or .tsv file")
    parser.add_argument("-e", "--excel-file", dest="excel_file", required=True, help="The input excel file.")
    parser.add_argument("-o", "--ouput-pre", dest="output_pre", default="output", help="The output csv file. Default: %(default)s")
    parser.add_argument("-i", "--sheet-index", dest="sheet_index", nargs="*", help="The index of the sheet to convert.  Default: %(default)s")
    parser.add_argument("-f", "--sheet-format", dest="sheet_fmt", default="csv", choices=["csv", "tsv"], help="The format of the sheet to convert.  Default: %(default)s")

    return parser


def main():
    args = get_args().parse_args()

    excel_file = args.excel_file
    output_pre = args.output_pre
    sheet_index = args.sheet_index
    sheet_fmt = args.sheet_fmt

    if sheet_fmt == "csv":
        splitter = ","
    else:
        splitter = "\t"

    dtfm_pool = pd.read_excel(excel_file, sheet_name=sheet_index)
    if isinstance(dtfm_pool, dict):
        for key, val in dtfm_pool.items():
            sheet_name = key.replace(" ", "_")
            csv_opt_name = output_pre + "." + sheet_name + "." + sheet_fmt
            val.to_csv(csv_opt_name, sep=splitter, index=False)
    elif isinstance(dtfm_pool, pd.DataFrame):
        csv_opt_name = output_pre + "." + sheet_index + "." + sheet_fmt
        dtfm_pool.to_csv(csv_opt_name, sep=splitter, index=False)
    else:
        raise ValueError("Unknow type.")


if __name__ == "__main__":
    main()
