#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""QTL mapping pipeline utils.
"""

import re
import json
import logging

import numpy as np
import pandas as pd
import seaborn as sb
import pingouin as pg
import scipy.stats as stt
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

from sklearn.decomposition import PCA

# Set up a global logger.
logger = logging.getLogger("matplotlib")
logger.setLevel(logging.ERROR)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

cs_stream = logging.StreamHandler()
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}", datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_stream.setFormatter(fmt)
cs_stream.setLevel(logging.INFO)

logger.addHandler(cs_stream)


CHROM_LEN_GRCH37 = {
    'chr1' : 249250621, 'chr2' : 243199373, 'chr3' : 198022430,
    'chr4' : 191154276, 'chr5' : 180915260, 'chr6' : 171115067,
    'chr7' : 159138663, 'chr8' : 146364022, "chr9" : 141213431,
    "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
    "chr16":  90354753, "chr17":  81195210, "chr18":  78077248,
    "chr19":  59128983, "chr20":  63025520, "chr21":  48129895,
    "chr22":  51304566, "chrx" : 155270560, "chry" :  59373566,
}


class DataSet:
    """A class to handle data set.
    """
    def __init__(self, file_path=None):
        self.file_path = file_path

    def load_data(self, file_path, fmt='tsv'):
        """Load data from given file path"""
        return self

    def dump_data(self, file_path, fmt='tsv'):
        """Dump processed data into given file path"""
        return self


class BoxPlot:
    """Draw boxplot.
    """
    def __init__(self, phtp_file, phtp_name, gntp_file, snps_indx, gntp_info_file,
                 cvrt_file=None, cvrt_name=None, output_dir="./"):
        self.phtp_file = phtp_file
        self.phtp_name = phtp_name
        self.gntp_file = gntp_file
        self.snps_indx = snps_indx
        self.gntp_info_file = gntp_info_file
        self.cvrt_file = cvrt_file
        self.cvrt_name = cvrt_name

        self.output_dir = output_dir

        self.phtp_dtfm = None
        self.cvrt_dtfm = None
        self.gntp_dtfm = None
        self.gntp_info_dtfm = None
        self.dtfm = None

    @staticmethod
    def _pr_clear():
        plt.cla()
        plt.clf()
        plt.close()

    @staticmethod
    def _pr_load_file(file_path, **kwargs):
        sep = kwargs.pop("sep") if "sep" in kwargs else "\t"
        dtfm = pd.read_csv(file_path, sep=sep, **kwargs)
        return dtfm

    def _pr_make_dtfm(self):
        # Load phenotype data
        phtp_dtfm = self._pr_load_file(self.phtp_file, index_col=0)

        if self.phtp_name is None:
            self.phtp_dtfm = phtp_dtfm
            self.phtp_name = phtp_dtfm.index
        else:
            self.phtp_name = [phtp for phtp in self.phtp_name if phtp in phtp_dtfm.index]
            self.phtp_dtfm = phtp_dtfm.loc[self.phtp_name]

        # Load genotype data
        gntp_dtfm = self._pr_load_file(self.gntp_file, sep=" ", index_col=0)
        gntp_info_dtfm = self._pr_load_file(self.gntp_info_file, sep=" ", index_col=0)

        if self.snps_indx is None:
            self.gntp_dtfm = gntp_dtfm
            self.snps_indx = gntp_dtfm.index
            self.gntp_info_dtfm = gntp_info_dtfm
        else:
            self.snps_indx = [
                snps_indx for snps_indx in self.snps_indx
                if snps_indx in gntp_dtfm.index and snps_indx in gntp_info_dtfm.index]

            self.gntp_dtfm = gntp_dtfm.loc[self.snps_indx]
            self.gntp_info_dtfm = gntp_info_dtfm.loc[self.snps_indx]


        if len(self.snps_indx) > 16:
            raise ValueError("More than 16 snps were selected. Unsupported!!!")

        # Load covariates
        if self.cvrt_file is not None:
            cvrt_dtfm = self._pr_load_file(self.cvrt_file, index_col=0)

            if self.cvrt_name is None:
                self.cvrt_dtfm = cvrt_dtfm
                self.cvrt_name = cvrt_dtfm.index
            else:
                self.cvrt_name = [cvrt for cvrt in self.cvrt_name if cvrt in cvrt_dtfm.index]
                self.cvrt_dtfm = cvrt_dtfm.loc[self.cvrt_name]

        # Merge the previous three dataframe
        self.dtfm = pd.concat([self.phtp_dtfm, self.gntp_dtfm, self.cvrt_dtfm], join="inner")

    def _pr_encode_gntp(self):
        for gntp in self.snps_indx:
            eff, alt = self.gntp_info_dtfm.loc[gntp, ["EffectAllele", "AlternativeAllele"]]
            dosage2code_dict = {0: alt + alt, 1: alt + eff, 2: eff + eff}
            self.dtfm.loc[gntp + "_code"] = self.dtfm \
                    .loc[gntp, :] \
                    .apply(lambda x: dosage2code_dict[round(x)])

    def _pr_correct_cvrt(self):
        for phtp in self.phtp_name:
            dtfm = self.dtfm.loc[[phtp]]
            smf.ols('', data=dtfm).fit()

    def _pr_draw_boxplots(self, svfmt="pdf", **kwargs):
        sb.set(style="ticks")
        width = kwargs.pop("width") if "width" in kwargs else 10
        height = kwargs.pop("height") if "height" in kwargs else 10

        for phtp in self.phtp_name:
            for gntp in self.snps_indx:
                gntp_code = gntp + "_code"
                dtfm = self.dtfm.loc[[phtp, gntp_code]]

                chrom, pos, eff, alt, *_ = self.gntp_info_dtfm.loc[gntp, ]
                allele_code = [alt + alt, alt + eff, eff + eff]

                allele_count = dtfm.loc[gntp_code].value_counts()
                xtick_labels = ["{}({})".format(aa, allele_count[aa]) for aa in allele_code if aa in allele_count]
                allele_code = [aa for aa in allele_code if aa in allele_count]

                axes = sb.boxplot(x=gntp_code, y=phtp, data=dtfm.transpose(), width=0.4, order=allele_code)
                axes = sb.swarmplot(x=gntp_code, y=phtp, data=dtfm.transpose(), color=".5", order=allele_code)

                axes.set_title("{}({})".format(phtp, gntp))
                axes.set_xlabel("{}({},{},{}>{})".format(gntp, chrom, pos, alt, eff))
                axes.set_ylabel(phtp)
                axes.set_xticklabels(xtick_labels)

                fig_name = ".".join(["boxplot", phtp, gntp_code, svfmt])
                fig_name = self.output_dir.strip("/") + "/" + fig_name

                fig = axes.get_figure()
                fig.savefig(fig_name, width=width, height=height)

                self._pr_clear()

    def init(self):
        """Init."""
        self._pr_make_dtfm()
        self._pr_encode_gntp()
        return self

    def stats(self):
        """Do statistics."""
        self._pr_correct_cvrt()
        return self

    def draw_boxplot(self):
        """Draw boxplot."""
        self._pr_draw_boxplots()
        return self


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
        self.stdv_times = None
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
                col_to_func_dict.update({x: self._inverse_rank for x in val_vec})
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

    @staticmethod
    def _inverse_rank(raw_list):
        # Inverse rank transformer.
        return stt.norm.ppf(stt.rankdata(raw_list) / (len(raw_list) + 1))

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
        if pntp_to_be_transed:
            pntp_transed_dtfm = self.dataframe \
                    .loc[:, pntp_to_be_transed] \
                    .transform(phtp_trans_dict, axis=0)

            pntp_transed_name_list = []
            for col_name in pntp_transed_dtfm.columns:
                tran_func = phtp_trans_dict[col_name]
                func_name = tran_func if isinstance(tran_func, str) else "ivrk"
                pntp_transed_name_list.append(col_name + "_" + func_name)

            pntp_transed_dtfm.columns = pntp_transed_name_list
            self.dataframe[pntp_transed_name_list] = pntp_transed_dtfm
            self.pntp_transed_name_list = pntp_transed_name_list
        else:
            self.pntp_transed_name_list = []

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

        stdv_times = self.stdv_times
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


class CollectAndReportSNP:
    """Collect and report SNPs.
    """
    def __init__(self, input_file_pl):
        self.input_file_pl = input_file_pl

    def _load_snps(self, file_path=None):
        if file_path is None:
            file_path = self.input_file_pl

        if isinstance(file_path, list):
            for flp in file_path:
                pd.read_csv(flp)
        elif isinstance(file_path, str):
            pd.read_csv(file_path)


class PreapreDataForCircos:
    """Preapre QTL mapping results for the Circos plot.
    """
    def __init__(self, fppl):
        self.fppl = fppl
        self.snps_dfpl = None
        self.snps_dtfm = None

    def _load_snps(self, flpt, **kwargs):
        pass

    @staticmethod
    def _filter(dtfm: pd.DataFrame):
        return dtfm

    @staticmethod
    def _merge(dtfm_pl: list):
        return dtfm_pl

    def prepare(self, **kwargs):
        '''Preapre dataset.'''
        self.snps_dtfm = self._merge([self._filter(self._load_snps(flpt, **kwargs)) for flpt in self.fppl])
        return self

    def write(self, opt_path, **kwargs):
        '''Write the dataset to disk.
        '''
        self.snps_dtfm.to_csv(opt_path, **kwargs)
        return self


if __name__ == '__main__':
    logger.error("Import only module.")
