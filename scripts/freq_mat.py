# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, import-error, too-many-locals,c-extension-no-member
#This is a standalone script that generates a snp frequency matrix for all the refernce genomes
import os
import sys
import subprocess
import logging
from traceback import format_exc
from argparse import ArgumentParser
import ahocorasick
import pandas as pd
from scripts import simulate_ref_reads as sr
from scripts import db_generate as dbg
logging.basicConfig(level=logging.DEBUG, format='%(message)s')


def freqmat_gen(A: ahocorasick.Automaton, final_serovars: list, final_table: pd.DataFrame, fastqpath: str):
    """
    Generates k-mer frequency dataframe for all completed genomes in database
    :param A: Ahocorasick automaton with all the k-mers loaded in it
    :param final_serovars: list of genomes included in database
    :param final_table: a dataframe with SNP coverage for every genome in the database
    :param fastqpath: path for simulated fastq file all genomes
    """
    freq_mat = final_table[["POS", "positive_kmers"]]
    for file in (final_serovars[3:-2]):
        fastq = os.path.join(fastqpath, ((file + ".fastq")))
        if os.path.exists(fastq):
            logging.debug(" {} file already present".format(fastq))
        else:
            filename = file+".fasta"
            sr.ref_read_simulation(filename)
        file_i = dbg.kmer_info(A, fastq)
        freq_mat = pd.merge(freq_mat, file_i, left_on=['POS', 'positive_kmers'], right_on=[
                            'POS', 'kmer_seq'], how='left')
        freq_mat = freq_mat.rename(columns={col: '{}'.format(file)
                                            for col in ('kmer_freq', 'freq')})
    table_cols = [c for c in freq_mat.columns if c.lower()[:8] != 'kmer_seq']
    freq_matrix = freq_mat[table_cols]
    freq_matrix = freq_matrix.fillna(0)
    num = freq_matrix._get_numeric_data()
    num[num < 10.0] = 0
    freq_matrix.to_csv("freq_matrix.csv")
    cmd = "  ".join(["rm ", " -f ", fastqpath, "/*.fastq"])  # perform cleanup
    subprocess.call(cmd, shell=True)
    return freq_matrix


def main(argv=None):
    """
    Generates fastq files for all genomes in the database
    :return: A dataframe with k-mer frequency at every SNP position for all genomes in the database
    """

    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        parser = ArgumentParser(
            description='Generate kmer frequency matrix containing SNP information for all the reference genomes.')

        parser.add_argument('--ref_genome', required=True,
                            help='path and name of reference genome .fasta file')

        args = parser.parse_args()
        file_name = "freq_matrix.csv"
        cur_dir = os.getcwd()
        snps = dbg.snp_info()
        kmers_allpos, ref_allpos = dbg.kmerseq(snps, args.ref_genome)
        final_table,snps_uniq,unique, final_serovars = dbg.kmerdb(snps, kmers_allpos, ref_allpos)
        A = dbg.init_automaton(final_table)
        while True:
            file_list = os.listdir(cur_dir)
            if file_name in file_list:
                logging.debug("File Exists in: ", cur_dir)
                break
            else:
                logging.debug("preparing freq_matrix.csv")
                dir_name = "reffastq"
                fastqpath = os.path.join(cur_dir, dir_name)
                if os.path.isdir(fastqpath):
                    freqmat_gen(A, final_serovars, final_table, fastqpath)
                else:
                    os.mkdir(fastqpath)
                    freqmat_gen(A, final_serovars, final_table, fastqpath)
    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error


if __name__ == "__main__":
    sys.exit(main())
