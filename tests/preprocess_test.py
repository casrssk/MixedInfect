from unittest import TestCase
from unittest.mock import call, patch
from scripts.preprocess import pairwise_dist, subprocess, seq_align, os

class PreProcessTest(TestCase):

    VALID_FASTAFILES = "./data/genomes/Sweden2.fasta"
    EMPTY_FASTAFILES = ""
    REF_GEN = "./data/genomes/DUW3CX.fasta"

    @patch('subprocess.call')
    def pairwise_dist_valid(self, mock_sub):
        pairwise_dist(fastafiles=self.VALID_FASTAFILES)
        mock_sub.assert_has_calls(
                calls=[
                    call("mash sketch -s 400 -k 16  -o  reference ",VALID_FASTAFILES, shell=True),
                    call("mash dist -t ", "reference.msh ",
                                          VALID_FASTAFILES, " >", " distance.tsv", shell=True),
                ]
            )

    @patch('subprocess.call')
    def pairwise_dist_invalid(self, mock_sub):
        with self.assertRaises(ValueError):
            pairwise_dist(fastafiles=self.EMPTY_FASTAFILES)


    @patch('subprocess.call')
    def seq_align_valid(self, mock_sub):
        seq_align(fastafiles=self.VALID_FASTAFILES, ref_genome=self.REF_GEN)
        mock_sub.assert_has_calls(
                calls=[
                    call("parsnp -r  ", REF_GEN, " -d ", VALID_FASTAFILES, " -x  -p 2", shell=True),
                    call("harvesttools -i", filename, "-V", vcf_out", shell=True), call("cat",  vcf_out, "|",  "grep", "-v", "'##'", ">", out_file, shell=True),

                ]
            )

    @patch('subprocess.call')
    def seq_align_emptyfasta(self, mock_sub):
        with self.assertRaises(ValueError):
            seq_align(fastafiles=self.EMPTY_FASTAFILES, ref_genome=self.REF_GEN)

    @patch("os.listdir", return_value=[])
    @patch('subprocess.call')
    def seq_align_emptydir(self, mock_sub, mock_dir):
        with self.assertRaises(ValueError):
            seq_align(fastafiles=self.VALID_FASTAFILES, ref_genome=self.REF_GEN)


    @patch("os.listdir", return_value=["abc"])
    @patch('subprocess.call')
    def seq_align_nofilematch(self, mock_sub, mock_dir):
        with self.assertRaises(FileNotFoundError):
            seq_align(fastafiles=self.VALID_FASTAFILES, ref_genome=self.REF_GEN)
