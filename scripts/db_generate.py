# -*- coding: utf-8 -*-
# pylint: disable=line-too-long, missing-module-docstring, too-many-locals, c-extension-no-member,E0401

from collections import Counter, defaultdict
import os
import logging
import fnmatch
import ahocorasick
import pandas as pd
from scripts import file_parser as fp


# Clustering strains and selecting one representative strain from each cluster.



def p_dist(fastafile="distance.tsv.zip"):
    """
    Cluster analysis and selecting representative strains from
    each cluster to update the k-mer database
    :param fastafile: .fasta file path as a string
    :type str
    """
    # preprocess.pairwise_dist(fastafile)  # run mash
    mash_dist = pd.read_csv(fastafile, sep='\t')
    mash_dist = mash_dist.set_index("#query")
    corr = mash_dist.corr()
    identical_strains = {}
    for i in range(len(corr)):
        for j in range(len(corr)):
            if (corr.iloc[i, j] == 1.000000) and (i != j):
                identical_strains[corr.index[i]] = corr.columns[j]
    identical_col = identical_strains.keys()
    include_col = []
    for key, value in identical_strains.items():
        if not key in include_col and value not in include_col:
            include_col.append(key)
            include_col.append(value)
    selected_ser = [include_col[i] for i in range(len(include_col)) if i % 2 == 0]
    exclude_serovars = list((Counter(identical_col) - Counter(selected_ser)).elements())
    return exclude_serovars


# filter and quality check of the snps in the database

def snp_info(f_in="tab.out", path=os.getcwd()):
    """
    performing quality check for the identified snps from all reference genomes
    :param f_in: file path for file generated by parsnp
    :type str
    :param path: setting the path to current directory
    :type str
    """
    for file in os.listdir(path):
        if fnmatch.fnmatch(file, '*.out'):
            snp_output = pd.read_csv(file, sep="\t")
    snp_subset = snp_output.loc[(snp_output.FILTER == "PASS") & (snp_output.QUAL >= 30)]
    snp_subset = snp_subset.loc[snp_subset['ALT'].str.len() == 1]
    snp_subset = snp_subset.reset_index(drop=True)
    all_cols = snp_subset.columns.to_list()
    exclude_serovars = p_dist()
    # filtering the redundant strains from database
    selected_serovars = list((Counter(all_cols) - Counter(exclude_serovars)).elements())
    print("{} serovars selected for the k-mer database".format(len(selected_serovars)))
    select_cols = selected_serovars[1:5]+selected_serovars[9:]
    snp = snp_subset[select_cols]
    # strain_names = (snp.columns.str.split('.').str[0]).to_list()
    strain_names = (snp.columns.str.split('.').str[0]).to_list()
    snp = snp.set_axis(strain_names, axis=1, inplace=False)
    return snp

def kmerseq(snp:pd.DataFrame, ref_genome: str):
    """
    Combine all the snps with k-mer sequences; positive kmers with
    alternative snp and negative kmer with reference snp
    :param snp: enter the dataframe containing snp information
    :param ref_genome: path of reference genome
    """
    # preparing reference genome for extracting kmers
    with open(ref_genome,'r',encoding="utf-8") as fp_in:
        for name, seq in fp.read_fasta(fp_in):
            header, sequence = name, seq
    kmers_allpos = {}
    ref_allpos = {}
    positions = snp["POS"].tolist()
    alt_snps = snp["ALT"].tolist()  # alternate snps
    ref_snps = snp['REF'].tolist()  # ref snps for the same position

    for pos, alt, ref in zip(positions, alt_snps, ref_snps):
        kmer = sequence[pos-11:pos+10]
        alt_seq = kmer[:10]+alt+kmer[11:]
        kmers_allpos[pos] = alt_seq
        ref_seq = kmer[:10]+ref+kmer[11:]
        ref_allpos[pos] = ref_seq
    return kmers_allpos, ref_allpos


# snps = snp_info("tab.out", path=os.getcwd())


def kmerdb(snp, kmers_allpos, ref_allpos):
    """
    Creating the kmer database with snp information from all gennomes
    :param snps: kmer database with snp information
    :type dataframe

    """
    snp_cols = snp.columns.to_list()
    cols = snp_cols[0:1]+snp_cols[5:]
    snp_table = snp[cols]
    positive_kmers = snp_table['POS'].map(kmers_allpos)
    negative_kmers = snp_table['POS'].map(ref_allpos)
    snp_table.insert(loc=1, column='positive_kmers', value=positive_kmers)
    snp_table.insert(loc=2, column='negative_kmers', value=negative_kmers)
    # count the number of ref genomes
    pd.options.mode.chained_assignment = None  # default='warn'
    snp_table["count_positive"] = snp_table.loc[:, snp_cols[5:]].sum(axis=1)
    snp_table["count_negative"] = len(snp_cols[5:]) - snp_table["count_positive"]
    snp_table = snp_table.loc[:, ~snp_table.columns.duplicated()]
    col_list = snp_table.columns.to_list()

    # #creating a reference snp database to remove bias due to reference genome
    n_refs = len(col_list[3:-2])
    d_snp_count = int(n_refs * 0.75)
    d_snps = snp_table.loc[(snp_table["count_positive"] > d_snp_count) &
                           (snp_table["count_negative"] < (n_refs) - d_snp_count)]
    # closest to ref_genome with a difference of 11 snps and so not to lose the unique snps change 1 to 0
    d_snps["D-EC"].replace(1, 0, inplace=True)
    d_snps = d_snps.replace({0: 1, 1: 0})

    columns_titles = d_snps.columns.to_list()
    columns_titles[1], columns_titles[2] = columns_titles[2], columns_titles[1]
    d_snps.columns = columns_titles
    col_swap = columns_titles[0:1] + columns_titles[2:3]+columns_titles[1:2]+columns_titles[3:]
    d_snps = d_snps[col_swap]

    final_table = pd.concat([snp_table, d_snps]).drop_duplicates(
        ["POS"], keep='last').sort_values(["POS"])
    final_table = final_table.reset_index(drop=True)
    final_table["count_positive"] = final_table.loc[:, col_list[3:-2]].sum(axis=1)
    final_table["count_negative"] = len(col_list[3:-2]) - final_table["count_positive"]
    # get the snps uniques to all the serovars in snps_uniq dataframe
    snps_uniq = final_table.loc[(final_table['count_positive'] == 1)]
    snps_uniq.loc['Total', col_list[3:-2]] = snps_uniq.sum(axis=0)
    # # snps_uniq.to_csv("/Users/satwantkaur/Desktop/mixed_infect/unique.csv")
    unique = (snps_uniq.loc["Total", col_list[3:-2]]).to_dict()
    # zero is included as it accounts for the ref genome
    low_count = [k for k, v in unique.items() if (v > 0) and (v < 5)]
    final_serovars = [x for x in col_list if x not in low_count]
    final_table = final_table[final_serovars]
    final_table.to_csv("kmerdb.csv")
    snps_uniq = snps_uniq[final_serovars]
    snps_uniq.to_csv("uniquesnpdb.csv")
    unique = (snps_uniq.loc["Total", final_serovars[3:-2]]).to_dict()
    print("The number of reference genomes in K-mer database is:", len(final_serovars[3:-2]))
    return final_table,snps_uniq,unique,final_serovars


def init_automaton(final_table: pd.DataFrame):
    """
    Initialises Aho-Corasick Automaton with kmers with SNP information from reference fasta
    :param fasta file: .fasta file path as a string
    :return: Aho-Corasick Automaton with kmers loaded
    """
    A = ahocorasick.Automaton()
    kmers_list = dict(zip(final_table.POS, final_table.positive_kmers))
    # kmers_list = list(kmers_allpos.values())
    for i, (pos, seq) in enumerate(kmers_list.items()):
        A.add_word(seq, (pos, seq, i))
        A.add_word(fp.rev_comp(seq), (pos, seq, i))
    A.make_automaton()
    return A


def kmer_info(Aho: ahocorasick.Automaton, fastq: str) -> pd.DataFrame:
    """
    Finds k-mers in the input fastq files
    :param Aho: Ahocorasick automaton with all the k-mers loaded in it
    :param fastq: filepath for the input fastq file

    :return: k-mer frequency at SNP positions found in test fastq
    """
    kmer_seq_counts = defaultdict(int)
    kmer_df = pd.DataFrame(columns=['POS', 'kmer_seq', 'freq'])
    for _, sequence in fp.parse_fastq(fastq):
        for idx, (_, kmer_seq, _) in Aho.iter(sequence):
            kmer_seq_counts[kmer_seq] += 1
    res = []
    for kmer_seq, freq in kmer_seq_counts.items():
        kmername, sequence, _ = Aho.get(kmer_seq)
        res.append((kmername, kmer_seq, freq))

    def f_out(val, index): return tuple(i[index] for i in val)
    tup1 = f_out(res, 0)
    tup2 = f_out(res, 1)
    tup3 = f_out(res, 2)
    for x in range(len(res)):
        kmer_df = kmer_df.append(
            {'POS': tup1[x], 'kmer_seq': tup2[x], 'freq': tup3[x]}, ignore_index=True)
    return kmer_df
