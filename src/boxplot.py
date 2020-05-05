#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sbn
import statsmodels.formula.api as smf

logger = logging.getLogger("matplotlib")
logger.setLevel(logging.ERROR)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

cs_stream = logging.StreamHandler()
fmt = logging.Formatter("| {levelname: ^8} | {asctime} | {name}: {message}", datefmt="%Y-%m-%d %H:%M:%S %p", style="{")
cs_stream.setFormatter(fmt)
cs_stream.setLevel(logging.DEBUG)

logger.addHandler(cs_stream)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phtp-file", dest="phtp_file", required=True, help="The file read the phenotype level from.")
    parser.add_argument("-g", "--gntp-file", dest="gntp_file", required=True, help="The file read the genotype from.")
    parser.add_argument("-i", "--gntp-info-file", dest="gntp_info_file", required=True, help="The file read the genotype information from.")
    parser.add_argument("-c", "--cvrt-file", dest="cvrt_file", help="The file read the covriates from.")
    parser.add_argument("-P", "--phtp-name", dest="phtp_name", nargs="*", help="The name of phenotype for which the script draws the boxplot.")
    parser.add_argument("-G", "--snps-indx", dest="snps_indx", nargs="*", help="The genotype id for which the script draw boxplot.")
    parser.add_argument("-C", "--cvrt-name", dest="cvrt_name", nargs="*", help="The name of covriates by which the script correct.")
    parser.add_argument("-o", "--output-prefix", dest="output_dir", help="The prefix for the output aka boxplot")

    return parser


class BoxPlot:

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
            self.phtp_dtfm = phtp_dtfm.loc[self.phtp_name, :]

        # Load genotype data
        gntp_dtfm = self._pr_load_file(self.gntp_file, sep=" ", index_col=0)
        gntp_info_dtfm = self._pr_load_file(self.gntp_info_file, sep=" ", index_col=0)
        if self.snps_indx is None:
            self.gntp_dtfm = gntp_dtfm
            self.snps_indx = gntp_dtfm.index
            self.gntp_info_dtfm = gntp_info_dtfm
        else:
            self.gntp_dtfm = gntp_dtfm.loc[self.snps_indx, :]
            self.gntp_info_dtfm = gntp_info_dtfm.loc[self.snps_indx, :]


        if len(self.snps_indx) > 16:
            raise ValueError("More than 16 snps were selected. Unsupported!!!")

        # Load covariates
        if self.cvrt_file is not None:
            cvrt_dtfm = self._pr_load_file(self.cvrt_file, index_col=0)
            if self.cvrt_name is None:
                self.cvrt_dtfm = cvrt_dtfm
                self.cvrt_name = cvrt_dtfm.index
            else:
                self.cvrt_dtfm = cvrt_dtfm.loc[self.cvrt_name, :]

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
            model = smf.ols('', data=dtfm).fit()

    def _pr_draw_boxplots(self, svfmt="pdf", **kwargs):
        sbn.set(style="ticks")
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

                axes = sbn.boxplot(x=gntp_code, y=phtp, data=dtfm.transpose(), width=0.4, order=allele_code)
                axes = sbn.swarmplot(x=gntp_code, y=phtp, data=dtfm.transpose(), color=".5", order=allele_code)

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
        self._pr_make_dtfm()
        self._pr_encode_gntp()
        return self

    def stats(self):
        self._pr_correct_cvrt()
        return self

    def draw_boxplot(self):
        self._pr_draw_boxplots()
        return self


def main():
    args = get_args().parse_args()
    phtp_file = args.phtp_file
    phtp_name = args.phtp_name
    gntp_file = args.gntp_file
    snps_indx = args.snps_indx
    gntp_info_file = args.gntp_info_file
    cvrt_file = args.cvrt_file
    cvrt_name = args.cvrt_name

    output_dir = args.output_dir

    box_plot = BoxPlot(
        phtp_file, phtp_name, gntp_file, snps_indx, gntp_info_file, cvrt_file,
        cvrt_name, output_dir=output_dir)
    box_plot.init().draw_boxplot()


if __name__ == "__main__":
    main()
