# -*- coding: utf-8 -*-
# pylint: disable=line-too-long , missing-module-docstring
import fnmatch
import subprocess
import os
import sys
import logging
from traceback import format_exc
from argparse import ArgumentParser
from datetime import datetime
logging.basicConfig(level=logging.DEBUG, format='%(message)s')
# os.chdir("data/")


def pairwise_dist(fastafiles: str):
    """
    Calculates pairwise genome distance between candidate strains using mash
    :param fastafiles: path of the directory containing genome fasta files

    """
    if not fastafiles or not isinstance(fastafiles, str):
        raise ValueError("Fastafiles is not valid or empty!")
    mash_ref = "  ".join(["mash sketch -s 400 -k 16  -o  reference",  fastafiles])
    subprocess.call(mash_ref, shell=True)
    mash_dist = " ".join(["mash dist -t ", "reference.msh ",
                          fastafiles, " >", " distance.tsv"])
    subprocess.call(mash_dist, shell=True)
    logging.debug("SNP distance matrix saved in distance.tsv in the current directory")

# Performs core genome multi alignment and generates k-mer database


def seq_align(fastafiles: str, ref_genome: str):
    """
    Performs intra-specific  core genome alignment and generates SNP based k-mer database
    :param fastafiles: list of path to all fasta files
    :param ref_genome: path to reference genome
    """
    # ref_genome = "data/genomes/DUW3CX.fasta"
    if not fastafiles or not isinstance(fastafiles, str):
        raise ValueError("Fastafiles is not valid or empty!")
    year = datetime.now().year
    multiseq_align = " ".join(["parsnp -r  ", ref_genome, " -d ", fastafiles, " -x  -p 2"])
    subprocess.call(multiseq_align, shell=True)
    vcf_out = "output.vcf"
    out_file = "tab.out"
    cur_dir = os.getcwd()
    file_list = os.listdir(cur_dir)
    if not file_list:
        raise ValueError("Directory is empty!")
    for file in file_list:
        if fnmatch.fnmatch(file, 'P_'+str(year)+'_*'):
            filename = file+"/parsnp.ggr"
            logging.debug(filename)
    if not filename:
        raise FileNotFoundError("Expected file not found")
    # Convert the .ggr snp file into vcf and then tab file
    ggr2vcf = " ".join(["harvesttools -i", filename, "-V", vcf_out])
    subprocess.call(ggr2vcf, shell=True)
    vcf2tab = " ".join(["cat",  vcf_out, "|",  "grep", "-v", "'##'", ">", out_file])
    subprocess.call(vcf2tab,  shell=True)


def main(argv=None):
    """
    Allows to form clusters of highly related genomes, select one represenattive
    from each genome and avoid redundancy.

    """
    program_name = os.path.basename(sys.argv[0])
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        parser = ArgumentParser(
            description='Generate SNP distance matrix and k-mer database from all reference genomes')

        parser.add_argument('-d', dest='dbgenomes', nargs='+',  type=str, required=True,
                            help='Files you want to use to make the kmer database, Can be either .fa, .fna, .fasta')
        parser.add_argument('--ref_genome', required=True,
                            help='path and name of reference genome .fasta file')

        args = parser.parse_args()

        input_files = args.dbgenomes
        input_files_processed = []
        for input_f in input_files:
            if os.path.isdir(input_f):
                for f_in in os.listdir(input_f):
                    f_in = os.path.join(input_f, f_in)
                    if os.path.isfile(f_in) and (f_in.endswith(".fasta") or f_in.endswith(".fa") or f_in.endswith(".fna")):
                        input_files_processed.append(f_in)
            elif os.path.isfile(input_f) and (input_f.endswith(".fasta") or input_f.endswith(".fa") or input_f.endswith(".fna")):
                input_files_processed.append(input_f)
            else:
                logging_debug("{} is not a valid file".format(input_f))
        input_files = input_files_processed
        if len(input_files) < 2:
            logging.debug("Less than 2 input sequences provided...")
            sys.exit(1)

        name = " ".join(input_files[:])
        pairwise_dist(name)
        seq_align(name, args.ref_genome)

    except Exception as error:
        print(program_name + ": " + repr(error) + '\n' + format_exc() + '\n')
        raise error


if __name__ == "__main__":
    sys.exit(main())
