#!/usr/bin/python

import os.path
import os
import pandas as pd
import csv


class Equalizer:
    """
    A program for standardizing GWAS outputs to have matching column names.
    This will make it easier to run all data through Winnow,
    and allows for files from multiple analysis tools to run simultaneously.

    Attributes:
        fastlmm_params : dict
            Parameters required for the Equalizer to run based off the
            output format of FaST-LMM.
        ridge_params : dict
            Parameters required for the Equalizer to run based off the
            output format of Ridge Predict.
        bayesr_params : dict
            Parameters required for the Equalizer to run based off the
            output format of BayesR.
        plink_params : dict
            Parameters required for the Equalizer to run based off the
            output format of PLINK.
        qxpak_params : dict
            Parameters required for the Equalizer to run based off the
            output format of QxPak.
        gemma_params : dict
            Parameters required for the Equalizer to run based off the
            output format of Gemma.
    """

    def __init__(self):
        self.fastlmm_params = {}
        self.ridge_params = {}
        self.bayesr_params = {}
        self.plink_params = {}
        self.qxpak_params = {}
        self.gemma_params = {}

    def equalize_gwas(self):
        pass

    def fastlmm_eq(self, lmm_output, snp="SNP", score="Pvalue", beta="", covar="", delim='\t'):
        """Creates a dictionary for FaST-LMM header-row labels.

        :param lmm_output: The filename of the GWAS output to be equalized.
        :param snp: The header value for the SNP column of the output.
        :param score: The header value for the score (p-value) column of the output.
        :param beta: The header value for the beta column of the output. Not required.
        :param covar: The header value for the covariance column of the output. Not required.
        :param delim: The delimited used in the GWAS output file. Should be comma, tab, or whitespace delimited.
        :return:
        """
        if not self.fastlmm_params:
            self.fastlmm_params = {
                "snp_cols": snp,
                "score_cols": score,
                "beta_cols": beta,
                "covar_cols": covar,
                "delimiter": delim
            }

        self.equalize(lmm_output, self.fastlmm_params)

        # Rename snp to 'SNP', score to 'PVAL', beta to 'BETA', covar to 'COVAR'

    def ridge_eq(self, ridge_output, snp, score, beta, covar, delim):
        """Creates a dictionary for Ridge Predict header-row labels.

        :param lmm_output: The filename of the GWAS output to be equalized.
        :param snp: The header value for the SNP column of the output.
        :param score: The header value for the score (p-value) column of the output.
        :param beta: The header value for the beta column of the output. Not required.
        :param covar: The header value for the covariance column of the output. Not required.
        :param delim: The delimited used in the GWAS output file. Should be comma, tab, or whitespace delimited.
        :return:
        """
        if not self.ridge_params:
            self.ridge_params = {
                "snp_cols": snp,
                "score_cols": score,
                "beta_cols": beta,
                "covar_cols": covar,
                "delimiter": delim
            }

        self.equalize(ridge_output, self.ridge_params)

    def bayesr_eq(self, bayesr_output, snp, score, beta, covar, delim):
        """Creates a dictionary for BayesR header-row labels.

        :param lmm_output: The filename of the GWAS output to be equalized.
        :param snp: The header value for the SNP column of the output.
        :param score: The header value for the score (p-value) column of the output.
        :param beta: The header value for the beta column of the output. Not required.
        :param covar: The header value for the covariance column of the output. Not required.
        :param delim: The delimited used in the GWAS output file. Should be comma, tab, or whitespace delimited.
        :return:
        """
        if not self.bayesr_params:
            self.bayesr_params = {
                "snp_cols": snp,
                "score_cols": score,
                "beta_cols": beta,
                "covar_cols": covar,
                "delimiter": delim
            }

        self.equalize(bayesr_output, self.bayesr_params)

    def plink_eq(self, plink_output, snp, score, beta, covar, delim):
        """Creates a dictionary for PLINK header-row labels.

        :param lmm_output: The filename of the GWAS output to be equalized.
        :param snp: The header value for the SNP column of the output.
        :param score: The header value for the score (p-value) column of the output.
        :param beta: The header value for the beta column of the output. Not required.
        :param covar: The header value for the covariance column of the output. Not required.
        :param delim: The delimited used in the GWAS output file. Should be comma, tab, or whitespace delimited.
        :return:
        """
        if not self.plink_params:
            self.plink_params = {
                "snp_cols": snp,
                "score_cols": score,
                "beta_cols": beta,
                "covar_cols": covar,
                "delimiter": delim
            }

        self.equalize(plink_output, self.plink_params)

    def qxpak_eq(self, qxpak_output, snp, score, beta, covar, delim):
        """Creates a dictionary for QxPak header-row labels.

        :param lmm_output: The filename of the GWAS output to be equalized.
        :param snp: The header value for the SNP column of the output.
        :param score: The header value for the score (p-value) column of the output.
        :param beta: The header value for the beta column of the output. Not required.
        :param covar: The header value for the covariance column of the output. Not required.
        :param delim: The delimited used in the GWAS output file. Should be comma, tab, or whitespace delimited.
        :return:
        """
        if not self.qxpak_params:
            self.qxpak_params = {
                "snp_cols": snp,
                "score_cols": score,
                "beta_cols": beta,
                "covar_cols": covar,
                "delimiter": delim
            }

        self.equalize(qxpak_output, self.qxpak_params)

    def gemma_eq(self, gemma_output, snp, score, beta, covar, delim):
        """Creates a dictionary for Gemma header-row labels.

        :param lmm_output: The filename of the GWAS output to be equalized.
        :param snp: The header value for the SNP column of the output.
        :param score: The header value for the score (p-value) column of the output.
        :param beta: The header value for the beta column of the output. Not required.
        :param covar: The header value for the covariance column of the output. Not required.
        :param delim: The delimited used in the GWAS output file. Should be comma, tab, or whitespace delimited.
        :return:
        """
        if not self.gemma_params:
            self.gemma_params = {
                "snp_cols": snp,
                "score_cols": score,
                "beta_cols": beta,
                "covar_cols": covar,
                "delimiter": delim
            }

        self.equalize(gemma_output, self.gemma_params)

    def equalize(self, output_files, param_dict):
        """Standardizes the column headers for GWAS output files for ease of use with Winnow.

        :param output_files:
        :param param_dict:
        :return:
        """
        # TODO work with on Discovery Environment instead of local machine
        pass

    def local_equalize(self, output_files, param_dict):
        for out_file in output_files:
            if os.path.splitext(out_file)[1] == ".txt":
                df = pd.read_table(os.path.abspath(out_file), header=0,
                                   delim_whitespace=True)
            elif os.path.splitext(out_file)[1] == ".csv":
                df = pd.read_csv(os.path.abspath(out_file), header=0, sep=",")
            else:
                df = pd.read_table(os.path.abspath(out_file), header=0,
                                   delim_whitespace=True)

            df.rename(columns={name: 'SNP' for name in param_dict['snp_cols']}, inplace=True)
            df.rename(columns={name: 'PVAL' for name in param_dict['score_cols']},
                      inplace=True)
            if param_dict['beta_cols'] is not None:
                df.rename(columns={name: 'BETA' for name in param_dict['beta_cols']},
                          inplace=True)
            if param_dict['covar_cols'] is not None:
                df.rename(columns={name: 'COVAR' for name in param_dict['covar_cols']},
                          inplace=True)

            # Next, write everything to a regular file, delimited based on command line argument
            df.to_csv(os.path.abspath(out_file), sep=param_dict['delimiter'], index=False, header=True,
                      quoting=csv.QUOTE_NONE)
