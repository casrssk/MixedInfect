# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring,C0103,E0401,R1732

import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
score=0

def fastq_to_dataframe(filename: str, size=100000):
    """
    Convert fastq to dataframe.
    :param size: limit to the first reads of total size
    :param filename:  filepath as string
    """

    ext = os.path.splitext(filename)[1]
    if ext == '.fastq':
        fastq_parser = SeqIO.parse(open(filename, "rt",encoding="utf-8"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r",encoding="utf-8"), "fastq")
    i = 0
    res = []
    for fastq_rec in fastq_parser:
        i += 1
        if i > size:
            break
        res.append([fastq_rec.id, str(fastq_rec.seq)])
    fq_df = pd.DataFrame(res, columns=['id', 'seq'])
    fq_df['length'] = fq_df.seq.str.len()
    return fq_df


def phred_score(file: str):
    """
    Calculates the phred score of the FASTQ file and generates warning if lower than the threshold
    :param file:  filepath as string
    """

    for record in SeqIO.parse(file, "fastq"):
        score = record.letter_annotations["phred_quality"]
    return np.mean(score)


def gc_content(file: str):
    """
    Calculates the GC content of the FASTQ file
    :param file:  filepath as string
    """

    fq_df = fastq_to_dataframe(file, size=100000)
    gc_count = fq_df.seq.apply(lambda x: GC(x))
    return np.mean(gc_count)
