# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, too-many-locals,E0401,C0116, invalid-name

import os
import sys
from argparse import ArgumentParser
from traceback import format_exc
import logging
import subprocess
import pandas as pd
from scripts import fastq_metrics as fm
from scripts import db_generate as dbg
from scripts import get_strains as gs

logging.basicConfig(level=logging.DEBUG, format='%(message)s')
results=[]

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
                            help='absolute path for the data folder ')
        parser.add_argument('--ref_genome', required=True,
                            help='absolute path and name of reference genome .fasta file')
        parser.add_argument('--genrate_db', required=False,
                            help='select 1 if interested in generating k-merdb dynamically or default is 0', default='0')

        args = parser.parse_args()
        snps = dbg.snp_info()
        kmers_allpos, ref_allpos = dbg.kmerseq(snps, args.ref_genome)
        final_table, snps_uniq, unique, final_serovars = dbg.kmerdb(snps, kmers_allpos, ref_allpos)
        logging.debug("initiating automaton")
        A = dbg.init_automaton(final_table)
        freq_matrix = pd.read_csv("freq_matrix.csv", sep=',')
        freq_matrix = freq_matrix.loc[:, ~freq_matrix.columns.str.contains('^Unnamed')]
        freq_matrix = freq_matrix.fillna(0)
        num = freq_matrix._get_numeric_data()
        num[num < 10.0] = 0
        logging.debug("frequency matrix ready!")
        input_path=os.path.join(args.path,"input")
        for filename in os.listdir(input_path):
            if filename.endswith(".fastq"):
                sample_fastq = os.path.join(input_path, (filename))
                sample_name = (os.path.splitext(filename)[0]).split(".")[0]
                logging.debug('executing "%s" ', sample_name)
                gc_content = fm.gc_content(sample_fastq)
                if gc_content < 35:
                    logging.debug(' GC_content of "%s" is lower than the threshold ', sample_fastq)
                phred = fm.phred_score(sample_fastq)
                if phred < 30:
                    logging.debug(' phred_score of "%s" is lower than the threshold ', sample_fastq)
                sample_df = dbg.kmer_info(A, sample_fastq)
                avg_kmer_freq = round(sample_df["freq"].mean(axis=0),2)
                if avg_kmer_freq < 20:
                    logging.debug('Averge k-mer coverage of the "%s" is less than the threshold of 20 ', filename)
                sam = sample_df.loc[(sample_df['freq'] >= 10)]
                sam = sam.reset_index(drop=True)
                subset = gs.strain_identity(sample_df, snps_uniq, freq_matrix,
                                         final_table, unique, final_serovars)
                X,  y_glm = gs.compute_matrix(sam, freq_matrix, subset)
                strain_abundance = gs.relabundance(X, y_glm, filename, avg_kmer_freq)
                if not any(elem is "None" for elem in strain_abundance):
                    for i in range(len(strain_abundance)):
                        if strain_abundance[i] > 0.0:
                            results.append(
                                (sample_name, subset[i+1], avg_kmer_freq, len(sample_df), strain_abundance[i]))
        output_dir=os.path.join(args.path,"output")
        if not os.path.isdir(output_dir):
            p=subprocess.Popen(' '.join(['mkdir ', output_dir]), shell = True)
        result_file=os.path.join(output_dir, "result.csv")
        outfile = pd.DataFrame(results, columns=["sample", "strain_id",
                                                 "Avg_kmer_freq", "Snps", "relabundance"])
        outfile.to_csv(result_file)

    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error


if __name__ == "__main__":
    sys.exit(main())
