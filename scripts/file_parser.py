# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring

import logging
import re
import os

allowed_bases = {'A', 'a', 'C', 'c', 'G', 'g', 'T',
                 't', 'N', 'n', 'X', 'x', }  # X for masked nucleotides

nuc = str.maketrans('acgtnxACGTNX', 'tgcanxTGCANX')


def rev_comp(seq: str) -> str:
    """
    Generates and return reverse complement of the sequnece
    :param seq: a string sequence

    :return: translated sequence
    """
    return seq.translate(nuc)[::-1]


REGEX_GZIPPED = re.compile(r'^.+\.gz$')


def parse_fasta(filepath: str):
    '''""
    Parse a .fasta/.fasta.gz file returning a generator yielding
    tuples of fasta headers to sequences.
    :param filepath: Fasta file path

    :return: generator: yields tuples of (<fasta header>, <fasta sequence>)
    '''""
    if REGEX_GZIPPED.match(filepath):
        logging.debug('Opening "%s" as gzipped file', filepath)
        with os.popen('zcat < {}'.format(filepath)) as fin:
            yield from _parse_fasta(fin, filepath)
    else:
        with open(filepath, 'r') as fin:
            yield from _parse_fasta(fin, filepath)


def _parse_fasta(fin: str, filepath: str):
    """
    Parse a .fasta file returning a generator yielding tuples
    of fasta headers to sequences
    :param fin: sequence fasta file
    :param filepath: .fasta file path as a string

    :return: generator: yields tuples of (<fasta header>, <fasta sequence>)
    """
    seqs = []
    header = ''
    line_count = 0
    for line in fin:
        if isinstance(line, bytes):
            line = line.decode()
        line = line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if header == '':
                header = line.replace('>', '')
            else:
                yield header, ''.join(seqs)
                seqs = []
                header = line.replace('>', '')
        else:
            invalid_bases = set(line) - allowed_bases
            if len(invalid_bases) > 0:
                msg = '{file}: Line {line} contains the following non-nucleotide characters: {chars}'.format(
                    file=filepath,
                    line=line_count,
                    chars=', '.join([str(x) for x in invalid_bases]))
                logging.warning(msg)
            seqs.append(line.upper())
        line_count += 1
    yield header, ''.join(seqs)

# getting all the snp position which are 16433


def read_fasta(fin: str):
    """
    Parse a .fasta file returning a generator yielding tuples of fasta headers to sequences
    :param fin: .fasta file path as a string

    """
    name, seq = None, []
    for line in fin:
        line = line.strip('\n')
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def parse_fastq(filepath: str):
    """
    Parse a .fatsq/.fastq.gz file
    :param filepath: .fastq/.fastq.gz filepath as string

    :return: returning a generator yielding
    tuples of fastq entry headers and sequences.
    """
    if REGEX_GZIPPED.match(filepath):
        logging.debug('Opening "%s" as gzipped file', filepath)
        with os.popen('zcat < {}'.format(filepath)) as fin:
            yield from _parse_fastq(fin)
    else:
        with open(filepath, 'r') as fin:
            yield from _parse_fastq(fin)


def _parse_fastq(fin: str):
    """
    A .fastq file parser which yields the header and sequence ignoring the quality scores
    :param fastq file
    """
    header = ''
    seq = ''
    skip = False
    for line in fin:
        if skip:
            skip = False
            continue
        line = line.strip()
        if line == '':
            continue
        if line[0] == '@':
            header = line.replace('@', '')
        elif line[0] == '+':
            yield header, seq
            skip = True
        else:
            seq = line.upper()
