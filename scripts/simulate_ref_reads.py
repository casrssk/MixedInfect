# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, invalid-name


import os
import subprocess
import time
# path for the regference genomes
wdir_in = os.getcwd()+"/data/genomes/"
# make directory for refernec genomes that you intend to use
file_path = "reffastq"
parent_dir = os.getcwd()
wdir_out = os.path.join(parent_dir, file_path)


def ref_read_simulation(file: str) -> None:
    """
    Simulates fastq files with the specified coverage for the given genome fasta file
    :param file: name of the fasta file
    """
    filepath1 = os.path.join(wdir_in, (file))
    index_cmd = "  ".join(["samtools faidx", filepath1])
    subprocess.call(index_cmd, shell=True)
    time.sleep(10)
    sample_name = (os.path.splitext(file)[0]).split(".")[0]
    filepath2 = os.path.join(wdir_out, (sample_name+"_1.fastq"))
    filepath3 = os.path.join(wdir_out, (sample_name+"_2.fastq"))
    reads = 500000
    sim_cmd1 = "  ".join(["mason_simulator  -ir", filepath1, "-n", str(reads),
                          "--illumina-read-length  200",  " -o", filepath2, " -or", filepath3])
    subprocess.call(sim_cmd1, shell=True)
    print(sim_cmd1)
    ref_sample = os.path.join(wdir_out, (sample_name+".fastq"))
    cat_cmd1 = "  ".join(["cat", filepath2, filepath3, ">",  ref_sample])
    print(cat_cmd1)
    subprocess.call(cat_cmd1, shell=True)
    time.sleep(10)
