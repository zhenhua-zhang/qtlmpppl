#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A script to pre-process traits for QTL mapping analysis.
"""

import argparse
import json
import logging

import matplotlib.pyplot as plt
import pandas as pd

logger = logging.getLogger("matplotlib")
logger.setLevel(logging.ERROR)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

cs_handle = logging.StreamHandler()
cs_handle.setLevel(logging.DEBUG)
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}",
                        datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_handle.setFormatter(fmt)

logger.addHandler(cs_handle)


def get_args():
    """Get CLI options.
    """
    parser = argparse.ArgumentParser(description="A script to preprocessing the phenotype data.")
    parser.add_argument("-p", "--phtp-file", dest="phtp_file", required=True, metavar="PHENOTYPE_FILE", help="The file from which read the phenotypes.")
    parser.add_argument("-t", "--traits", dest="traits", nargs="*", help="Traits to be processed. If is None, all numerical variable will be processed")
    parser.add_argument("-o", "--output-pre", dest="output_pre", default="phenotype", metavar="OUTPUT_PREFIX", help="Prefix of output file after processed. Default: %(default)s")
    parser.add_argument("-x", "--exp-sd", dest="expandding_sd", default=2, type=int, metavar="TIMES", help="Times of standard deviation. Default: %(default)s")
    parser.add_argument("-c", "--conf-file", dest="conf_file", metavar="CONFIG", help="The file from which read configurations. Default: %(default)s")
    parser.add_argument("--trans-func", dest="trans_fun", default="log2", choices=["log2", "log10", "ivr"], help="The transforming function applied on the given measurements. Default: %(default)s")
    parser.add_argument("--as-covar", dest="as_covar", nargs="*", metavar="COVARIATE", help="Traits used as covariates.")

    return parser


class PreProcess:
    """A class to pre-process phenotypes.
    """
    def __init__(self, phtp_file, conf_file=None, traits_to_proc=None):
        self.phtp_file = phtp_file
        self.conf_file = conf_file
        self.traits_to_proc = traits_to_proc

        self.configs = None
        self.dataframe = None
        self.numerical_cols = None

    def _pr_load_dtfm(self, **kwargs):
        if "sep" in kwargs:
            sep = kwargs.pop("sep")
        else:
            sep = "\t"

        self.dataframe = pd.read_csv(self.phtp_file, sep=sep, **kwargs)

    def _pr_load_config(self, **kwargs):
        self.configs = json.loads(self.conf_file, **kwargs)

    def init(self):
        """Initialize the processor.
        """
        self._pr_load_dtfm(sep=",")

        if self.conf_file:
            self._pr_load_config()
        else:
            logger.info("No configuration file.")

        self.numerical_cols = self.dataframe.dtypes[self.dataframe.dtypes != "object"]
        return self

    def check_dist(self):
        """Draw figures to show the distribution of the data.
        """
        for col_name in self.numerical_cols.index:
            fig, ax = plt.subplots()
            self.dataframe.loc[:, col_name].plot(ax=ax, kind="hist")
            hist_opt_name = col_name.replace(" ", "_").replace("/", "").replace("(", "").replace(")", "") + ".png"
            fig.savefig(hist_opt_name, figsize=(10, 10))
            fig.clear()
            plt.close()
        return self

    def check_outliers(self):
        """Perform a PCA analysis and show the results.

        Note:
            1. It works for dataset with at least 2 traits.
        """
        return self

    def remove_outliers(self, n_times_sd=3):
        """Remove outliers.
        """
        return self

    def transform(self, method=None, columns=None):
        """Transform the target columns using given math function.

        Note:
            1. The `pandas.DataFrame.transform()` method could be a good choice,
            as it allow to specify the function to apply on each column by a
            axis labels->functions implementation.
        """
        return self

    def transpose(self):
        """Transpose the loaded data-frame.

        Note: The data feed to MatrixEQTL should be individuals as columns and
        traits as rows.
        """
        return self

    def save_results(self):
        """Save the processed results into disk.
        """
        return self


def main():
    """The main entry of the module.
    """
    args = get_args().parse_args()

    traits = args.traits
    as_covar = args.as_covar
    trans_fun = args.trans_fun
    phtp_file = args.phtp_file
    conf_file = args.conf_file
    output_pre = args.output_pre
    expandding_sd = args.expandding_sd

    pre_processer = PreProcess(phtp_file, conf_file=conf_file, traits_to_proc=traits)
    pre_processer.init() \
            .check_dist() \
            .transform() \
            .check_outliers() \
            .remove_outliers() \
            .transpose() \
            .save_results()


if __name__ == "__main__":
    main()
