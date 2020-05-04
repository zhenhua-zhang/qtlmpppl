#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""A script to pre-process traits for QTL mapping analysis.
"""

import argparse
import json
import logging
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pingouin as pg
import seaborn as sb
from sklearn.decomposition import PCA

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
    parser.add_argument("-c", "--conf-file", dest="conf_file", metavar="CONFIG", help="The file from which read configurations. Default: %(default)s")
    parser.add_argument("-p", "--phtp-file", dest="phtp_file", metavar="PHENOTYPE_FILE", help="The file from which read the phenotypes.")
    parser.add_argument("-t", "--phtp-list", dest="phtp_list", nargs="*", help="Traits to be processed. If is None, all numerical variable will be processed")
    parser.add_argument("-v", "--cvrt-list", dest="cvrt_list", nargs="*", metavar="COVARIATE", help="Traits used as covariates.")
    parser.add_argument("-o", "--output-pre", dest="output_pre", default="phenotype", metavar="OUTPUT_PREFIX", help="Prefix of output file after processed. Default: %(default)s")
    parser.add_argument("-s", "--stdv-times", dest="stdv_times", default=3, type=int, metavar="TIMES", help="Times of standard deviation. Default: %(default)s")

    return parser


def inverse_rank(x):
    """Inverse rank transformer.
    """
    return x


class PreProcess:
    """A class to pre-process phenotypes.
    """
    def __init__(self, conf_file=None, phtp_file=None):
        self.conf_file = conf_file
        self.phtp_file = phtp_file

        self.configs = None   # Configurations parsed from self.conf_file
        self.dataframe = None # The data frame loaded from self.phtp_file

        self.phtp_list = None
        self.phtp_dtfm = None
        self.phtp_trans_dict = None
        self.pntp_transed_name_list = None

        self.sample_idx = None
        self.sample_id_col = None

        self.cvrt_list = None
        self.cvrt_dtfm = None

        self.pcor_dtfm = None

        self.transpose = False
        self.output_prefix = None

    def _pr_load_config(self, **kwargs):
        with open(self.conf_file, "r") as config_file_handle:
            self.configs = json.load(config_file_handle, **kwargs)

    def _pr_parse_phtp_list(self):
        self.phtp_list = self.configs["preprocess"]["phtp_list"]

    def _pr_parse_cvrt_list(self):
        self.cvrt_list = self.configs["preprocess"]["cvrt_list"]

    def _pr_parse_pntp_tran_dict(self):
        phtp_trans_dict = self.configs["preprocess"]["trans_dict"]
        _kept_cols = phtp_trans_dict.pop("orig")

        col_to_func_dict = {}
        transformed_cols_vec = []
        for key, val_vec in phtp_trans_dict.items():
            if key in ["log2", "log10"]:
                col_to_func_dict.update({x: key for x in val_vec})
            elif key in ["ivrk"]:
                col_to_func_dict.update({x: inverse_rank for x in val_vec})
            else:
                logger.error("Only log2, log10, and inverse rank are supported")
                continue

            transformed_cols_vec.extend([x + "_" + key for x in val_vec])
        self.phtp_trans_dict = col_to_func_dict

    def _pr_parse_misc(self):
        if self.phtp_file is None:
            self.phtp_file = self.configs["preprocess"]["phtp_file"]

        self.transpose = self.configs["preprocess"]["transpose"]
        self.output_prefix = self.configs["preprocess"]["output_prefix"]
        self.stdv_times = self.configs["preprocess"]["stdv_times"]
        self.sample_id_col = self.configs["preprocess"]["sample_id_col"]

    def _pr_parse_sample_idx(self):
        self.sample_idx = self.configs["preprocess"]["sample_idx"]

        sample_idx = []
        for idx_range in self.sample_idx.split(","):
            if idx_range:
                if "-" in idx_range:
                    start, stop = idx_range.split("-")
                    sample_idx.extend(list(range(int(start), int(stop))))
                else:
                    sample_idx.append(int(idx_range))
            else:
                logger.debug("Empty range of sample index")

        self.sample_idx = sample_idx

    def _pr_load_dtfm(self, **kwargs):
        sep = kwargs.pop("sep") if "sep" in kwargs else "\t"

        self.dataframe = pd.read_csv(self.phtp_file, sep=sep, **kwargs)
        if self.sample_idx is not None:
            self.dataframe = self.dataframe.loc[self.sample_idx, :]

    def init(self, **kwargs):
        """Initialize the processor.
        """
        if self.conf_file:
            self._pr_load_config()
        else:
            logger.info("No configuration file.")

        self._pr_parse_misc()
        self._pr_parse_phtp_list()
        self._pr_parse_cvrt_list()
        self._pr_parse_sample_idx()
        self._pr_parse_pntp_tran_dict()

        sep = kwargs.get("sep") if "sep" in kwargs else ","
        self._pr_load_dtfm(sep=sep)
        return self

    def encode_sex(self, sex_col="Sex", mapping=None):
        """Encode gender from string into binary"""
        mapping = {"Male": 0, "Female": 1} if mapping is None else mapping
        self.dataframe[sex_col] = self.dataframe[sex_col].apply(lambda x: mapping[x])

        return self

    def transform(self, skip=False):
        """Transform the target columns using given math function.

        Note:
            1. The `pandas.DataFrame.transform()` method could be a good choice,
            as it allow to specify the function to apply on each column by a
            axis labels->functions implementation.
        """
        if skip:
            return self

        phtp_trans_dict = self.phtp_trans_dict
        pntp_to_be_transed = phtp_trans_dict.keys()
        pntp_transed_dtfm = self.dataframe.loc[:, pntp_to_be_transed].transform(phtp_trans_dict, axis=0)

        pntp_transed_name_list = []
        for col_name in pntp_transed_dtfm.columns:
            tran_func = phtp_trans_dict[col_name]
            func_name = tran_func if isinstance(tran_func, str) else "ivrk"
            pntp_transed_name_list.append(col_name + "_" + func_name)

        pntp_transed_dtfm.columns = pntp_transed_name_list
        self.dataframe[pntp_transed_name_list] = pntp_transed_dtfm
        self.pntp_transed_name_list = pntp_transed_name_list

        self.dataframe.replace((np.inf, -np.inf), np.nan, inplace=True)

        return self

    def check_outliers(self, figfmt="png", skip=False):
        """Perform a PCA analysis and show the results.

        Note:
            1. It works for dataset with at least 2 traits.
        """
        if skip:
            return self

        pca = PCA()
        transformed_values = pca.fit_transform(self.dataframe.loc[:,
                                                                  self.phtp_list])

        fig, axes = plt.subplots()
        axes.scatter(transformed_values[:, 0], transformed_values[:, 1], s=0.5)

        pca_save_name = ".".join([self.output_prefix, "pca", figfmt])
        fig.savefig(pca_save_name)

        return self

    def mask_outliers(self, skip=False):
        """Remove outliers.
        """
        if skip:
            return self

        stdv_times= self.stdv_times
        def _mask_outliers(vec: pd.Series, stdv_times):
            vec_mean = vec.mean()
            vec_stdv = vec.std()
            upper = vec_mean + vec_stdv * stdv_times
            lower = vec_mean - vec_stdv * stdv_times
            vec[((lower > vec) | (vec > upper))] = np.nan

            return vec

        phtp_list = self.phtp_list
        self.dataframe.loc[:, phtp_list] = self.dataframe.loc[:, phtp_list].transform(_mask_outliers, stdv_times=stdv_times)

        return self

    def check_dist(self, figfmt="png", skip=False):
        """Draw figures to show the distribution of the data.
        """
        if skip:
            return self

        if self.pntp_transed_name_list is None:
            pntp_check_dist_list = self.phtp_list
        else:
            pntp_check_dist_list = self.phtp_list + self.pntp_transed_name_list

        for col_name in pntp_check_dist_list:
            fig, axes = plt.subplots()
            self.dataframe.loc[:, col_name].plot(ax=axes, kind="hist")

            hist_opt_name = re.sub("[()/]", "", col_name.replace(" ", "_"))
            hist_opt_name = ".".join([self.output_prefix, hist_opt_name, figfmt])
            fig.savefig(hist_opt_name)

            plt.cla()
            plt.clf()
            plt.close()
        return self

    def correlate(self, x_list=None, y_list=None, c_list=None, method="spearman", figfmt="png", skip=False):
        """Correlate variables.
        """
        if skip:
            return self

        x_list = self.phtp_list if x_list is None else x_list
        y_list = self.phtp_list if y_list is None else y_list
        c_list = self.cvrt_list if c_list is None else c_list

        x_len, y_len = len(x_list), len(y_list)
        pcorr_mtrx_np = np.eye(x_len, y_len)
        for idx_x, pntp_x in enumerate(x_list):
            for idx_y, pntp_y in enumerate(y_list):
                if pntp_y != pntp_x:
                    _pcorr_dtfm = pg.partial_corr(self.dataframe, pntp_x, pntp_y, self.cvrt_list, method=method)
                    pcorr_mtrx_np[idx_x, idx_y] = _pcorr_dtfm['r']

        self.pcor_dtfm = pd.DataFrame(pcorr_mtrx_np, index=self.phtp_list,
                                         columns=self.phtp_list)
        pcorr_htmp_name = ".".join([self.output_prefix, "correlation_heatmap", figfmt])
        ctmp_grid = sb.clustermap(self.pcor_dtfm, col_cluster=True, row_cluster=True, cmap="Greens")
        ctmp_grid.fig.savefig(pcorr_htmp_name)

        plt.cla()
        plt.clf()
        plt.close()

        return self

    def save_results(self, dtfm_fmt="tsv", skip=False):
        """Save the processed results into disk.
        """
        if skip:
            return self

        if dtfm_fmt == "tsv":
            sep = "\t"
        elif dtfm_fmt == "csv":
            sep = ","
        elif dtfm_fmt == "ssv":
            sep = " "
        else:
            logger.info("Only ssv, tsv and csv are supported, use tsv by defult")
            dtfm_fmt, sep = "tsv", "\t"

        if self.sample_id_col:
            self.dataframe.index = self.dataframe.loc[:, self.sample_id_col]

        phtp_tosave_list = [x for x in self.phtp_list if x not in self.phtp_trans_dict] + self.pntp_transed_name_list
        self.phtp_dtfm = self.dataframe.loc[:, phtp_tosave_list]
        phtp_dtfm_name = ".".join([self.output_prefix, "proc_phtp", dtfm_fmt])

        if self.transpose:
            self.phtp_dtfm = self.phtp_dtfm.transpose()

        self.phtp_dtfm.to_csv(phtp_dtfm_name, sep=sep)

        if self.cvrt_list:
            cvrt_cols_to_save_list = self.cvrt_list
            self.cvrt_dtfm = self.dataframe.loc[:, cvrt_cols_to_save_list]

            if self.transpose:
                self.cvrt_dtfm = self.cvrt_dtfm.transpose()

            cvrt_dtfm_name = ".".join([self.output_prefix, "proc_cvrt", dtfm_fmt])
            self.cvrt_dtfm.to_csv(cvrt_dtfm_name, sep=sep)

        if self.pcor_dtfm is not None:
            pcorr_mtrx_name = ".".join([self.output_prefix, "pcorr_phtp", dtfm_fmt])
            self.pcor_dtfm.to_csv(pcorr_mtrx_name, sep=sep, index=True)

        return self


def main():
    """The main entry of the module.
    """
    args = get_args().parse_args()

    phtp_file = args.phtp_file
    conf_file = args.conf_file

    pre_processer = PreProcess(conf_file=conf_file, phtp_file=phtp_file)
    pre_processer.init() \
            .encode_sex() \
            .transform() \
            .check_outliers() \
            .mask_outliers() \
            .check_dist() \
            .correlate() \
            .save_results()


if __name__ == "__main__":
    main()
