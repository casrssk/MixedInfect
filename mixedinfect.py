# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, too-many-locals, invalid-name

import os
import sys
from collections import defaultdict
from argparse import ArgumentParser
from traceback import format_exc
import ahocorasick
import pandas as pd
import numpy as np
import scipy.optimize as opt
from scripts import fastq_metrics as fm
from scripts import db_generate as dbg
from scripts import file_parser as fp

# import scripts.file_parser


def strain_identity(sample_df: pd.DataFrame, snps_uniq: pd.DataFrame, freq_matrix: pd.DataFrame, final_table: pd.DataFrame, unique: dict, final_serovars: list):
    """
    function to identify serovars from the input fastq file
    :param sample_df: sample input file path
    :param snps_uniq: strain specific signatures database path
    :param freq_matrix: path for kmer freqency database of all reference genomes
    :param final_table: path for the final snp database
    :param final_serovars: path for the list of genomes in the kmer database

    :return: list of identified serovars
    """
    def snpcov(sam: pd.DataFrame, final_table: pd.DataFrame, final_serovars: list) -> dict:

        # get kmer frequency coverage for top five strains
        snps_cov = pd.merge(sam, final_table, left_on=["POS", ], right_on=["POS", ], how='inner')
        snps_cov.loc['Total', final_serovars[3:-2]] = snps_cov.sum(axis=0)
        snps_covered = snps_cov.loc["Total", final_serovars[3:-2]].to_dict()
        snps_covered = {k: v for k, v in sorted(
            snps_covered.items(), key=lambda item: item[1], reverse=True)}
        top_freq = {k: snps_covered[k] for k in list(snps_covered)[:5]}
        strain_cov = {key: (top_freq[key]/strain_mapping_abundance.get(key, 0))
                      for key in top_freq if key in strain_mapping_abundance}
        freq = sum(list(strain_cov.values()))
        strain_sig = dict((k, v/freq) for k, v in strain_cov.items())
        return strain_sig

    def freq_calc(mydict1: dict, mydict2: dict) -> dict:
        """
        Calculate strain presence score in the sample based on kmer coverage
        :param mydict1: dictionary containing abundance of strain specific signatures in the sample
        :param mydict2: dictionary containing strains with most kmer coverage in the sample
        """
        serovar_rev = set()
        for k, v in mydict2.items():
            if mydict1[k] < 0.01 or mydict2[k] < 0.10:
                serovar_rev.add(k)
        for unwanted_key in serovar_rev:
            del mydict2[unwanted_key]
        freq = sum(list(mydict2.values()))
        strain_found = dict((k, v/freq) for k, v in mydict2.items())
        return strain_found

    temp = final_table.copy(deep=True)
    temp.loc['Total', final_serovars[3:-2]] = temp.sum(axis=0)
    # making a dictionary with total no of snp positions in each serovar
    strain_snps = (temp.loc['Total', final_serovars[3:-2]]).to_dict()
    # getting mapping abundance for each serovar: total snp positions/length of the matrix
    strain_mapping_abundance = {k: round(v / (len(final_table)), 3) for k, v in strain_snps.items()}

    samp_serovar = {}
    match_prop = pd.merge(sample_df, snps_uniq, left_on=["POS", ], right_on=["POS", ], how='inner')
    match_prop.loc['Total', final_serovars[3:-2]] = match_prop.sum(axis=0, skipna=True)
    total = (match_prop.loc["Total", final_serovars[3:-2]]).to_dict()
    snp_proportion = {key: (total[key]/unique.get(key, 0)) for key in total if key in unique}
    sam = sample_df.loc[(sample_df['freq'] >= 10)]
    sam = sam.reset_index(drop=True)
    # identifying top five strains with unique snp signatures
    topfive = sorted(snp_proportion, key=snp_proportion.get, reverse=True)[:5]
    serovar_rev = set()
    strain_high = []
    strain_low = set()
    for strain in topfive:
        value = snp_proportion.get(strain)
        if value < 0.01:
            serovar_rev.add(strain)
        else:
            samp_serovar[strain] = value
    for unwanted_key in serovar_rev:
        topfive.remove(unwanted_key)
    for k, v in samp_serovar.items():
        if samp_serovar[k] > 0.50:
            strain_high.append(k)
        else:
            strain_sig = snpcov(sam, freq_matrix, final_serovars)
            strain_select = freq_calc(snp_proportion, strain_sig)
            if k in strain_select:
                strain_low.add(k)

    subset = ["POS"] + strain_high + list(strain_low)
    return subset


def compute_matrix(sam: pd.DataFrame, freq_matrix: pd.DataFrame, subset: list):
    """
    Computes design matrix containing k-mer frequency at SNP positions in identified serovars
    :param sam: filtered sample file
    :param freq_matrix: frequency matrix for all genomes in kmer database
    :param subset: identifies serovars

    :return a design matrix and vector with observed coverage
    """

    comp_freq = freq_matrix[subset]
    kmer_freq = pd.merge(comp_freq, sam,  on=["POS"], how="left")
    kmer_freq = kmer_freq.fillna(0)
    kmer_freq = kmer_freq.drop(
        kmer_freq[(kmer_freq["kmer_seq"] == 0) & (kmer_freq["freq"] == 0)].index)
    col_names = subset[1:]
    kmer_freq = kmer_freq[(kmer_freq.loc[:, col_names[0:]] != 0.0).any(axis=1)]
    kmer_freq = kmer_freq.reset_index(drop=True)
    X_arr = []
    col_list = kmer_freq.columns.to_list()
    idx = col_list.index("POS")
    start = idx + 1

    for i in range(start, len(subset)):
        X_arr.append((kmer_freq.iloc[:, i]).values)
    X = (np.array(X_arr)).T
    y_glm = ((kmer_freq[["freq"]]).values).T
    return X,  y_glm


def relabundance(x: np.array, y: np.array, filename: str, avg_kmer_freq: int):
    """
    Calculate corrected abundances given a design matrix and observaed freqeuncy of SNPs using optimization.
    :param design_mat: [np.array (N,M)]: N is the kmer frequency at SNP positions in M strains
    :param y [np.array (N,)]: kmer frequency observed in the sample
    :param filename:name of test file
    :param avg_kmer_freq: average k-mer frequency at SNP positions

    :return: corrected relative abundances [numpy.array (M,)]: estimated abundance of each strain in the sample
    """
    # Make sure design matrix has the right orientation
    try:
        if x.shape[1] > x.shape[0]:
            x = x.transpose()
        ind_variable = range((x.shape[1]))
    except:
        print("Kmer coverage  of {0} is {1} and strains cannot be predicted".format(
            filename, avg_kmer_freq))
        return ["None"], [0.0]

    def objective(beta):
        return np.sum(np.square(np.dot(x, beta)-y))

    def build_con(k):
        # defining constraint 1: all the components of x should be > 0
        return lambda beta: beta[k]
    cons = [build_con(k) for k in ind_variable]
    # defining constraint 2: the sum of all components of x should not exceed 1
    cons.append(lambda beta: 1-np.sum(beta))
    # basic assunmption
    beta_0 = np.array([0.5 for i in range(x.shape[1])])
    # apply optimisation function
    abundances = opt.fmin_cobyla(objective, beta_0, cons, disp=0, rhoend=1e-10, maxfun=10000)
    total = 0
    for i in range(len(abundances)):
        total += abundances[i]
    rel_abundances = np.round((np.divide(abundances, total)), 2)

    return rel_abundances


results = []


def openfile(A: ahocorasick.Automaton, path: str, snps_uniq: pd.DataFrame, freq_matrix: pd.DataFrame, final_table: pd.DataFrame, unique: dict, final_serovars: list):
    """
    Identifies strain diversity and their relative abindances
    :param path: directory path for the test fastq files
    :param ref_genome: path  to the reference genome

    :return a list with sample_name, strains names, average k-mer coverage, SNPS found, relative abundances and  corelation co-efficient.
    """
    for filename in os.listdir(path):
        if filename.endswith(".fastq"):
            sample_fastq = os.path.join(path, (filename))
            sample_name = (os.path.splitext(filename)[0]).split(".")[0]
            print("executing {}".format(sample_name))
            gc_content = fm.gc_content(sample_fastq)
            if gc_content < 35:
                print("GC_content of {} is lower than the threshold".format(sample_fastq))
            phred = fm.phred_score(sample_fastq)
            if phred < 30:
                print("phred_score of {} is lower than the threshold".format(sample_fastq))
            sample_df = dbg.kmer_info(A, sample_fastq)
            avg_kmer_freq = sample_df["freq"].mean(axis=0)
            if avg_kmer_freq < 20:
                print("Averge k-mer coverage of the {} is less than the threshold of 20".format(filename))
            sam = sample_df.loc[(sample_df['freq'] >= 10)]
            sam = sam.reset_index(drop=True)
            subset = strain_identity(sample_df, snps_uniq, freq_matrix,
                                     final_table, unique, final_serovars)
            X,  y_glm = compute_matrix(sam, freq_matrix, subset)
            strain_abundance = relabundance(X, y_glm, filename, avg_kmer_freq)
            if not any(elem is "None" for elem in strain_abundance):
                for i in range(len(strain_abundance)):
                    if strain_abundance[i] > 0.0:
                        results.append(
                            (sample_name, subset[i+1], avg_kmer_freq, len(sample_df), strain_abundance[i]))
    return results


def main(argv=None):

    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    def dir_path(string):
        if os.path.isdir(string):
            return string
        else:
            raise NotADirectoryError(string)
    try:

        parser = ArgumentParser(
            description='Predict the starin diversity and realtive abundance of the identified strains')

        parser.add_argument('--path',   type=dir_path, required=True,
                            help='enter the path of the fastq test files. ')
        parser.add_argument('--ref_genome', required=True,
                            help='absolute path and name of reference genome .fasta file')
        # parser.add_argument('--genrate_db', required=False,
        #                     help='select 1 if interested in generating k-merdb dynamically or default is 0', default='0')

        args = parser.parse_args()
        snps = dbg.snp_info()
        kmers_allpos, ref_allpos = dbg.kmerseq(snps, args.ref_genome)
        final_table, snps_uniq, unique, final_serovars = dbg.kmerdb(snps, kmers_allpos, ref_allpos)
        print("initiating automaton")
        A = dbg.init_automaton(final_table)
        freq_matrix = pd.read_csv("freq_matrix.csv", sep=',')
        freq_matrix = freq_matrix.loc[:, ~freq_matrix.columns.str.contains('^Unnamed')]
        freq_matrix = freq_matrix.fillna(0)
        num = freq_matrix._get_numeric_data()
        num[num < 10.0] = 0
        print("frequency matrix ready!")
        results = openfile(A, args.path, snps_uniq, freq_matrix,
                           final_table, unique, final_serovars)
        outfile = pd.DataFrame(results, columns=["sample", "strain_id",
                                                 "Avg_kmer_freq", "Snps", "relabundance"])
        outfile.to_csv("results.csv")

    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error


if __name__ == "__main__":
    sys.exit(main())
