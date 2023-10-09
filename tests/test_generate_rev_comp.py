"""
SBI129: Programming
DOPS 3: Write a set of unit tests to evaluate whether the output of a script or
set of functions is correct, including edge cases.
"""

import os
import sys
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from generate_rev_comp import generate_rev_comp as grc
from tests import TEST_DATA_DIR


class TestOpeningFile(unittest.TestCase):
    """
    Test errors raised correctly when giving FASTA file with no sequence
    headers or a string instead of a file
    """
    no_headers_fasta = os.path.join(TEST_DATA_DIR, "no_header.fa")

    def test_check_file_headers(self):
        self.assertRaises(
            SystemExit, grc.create_fasta_dict, self.no_headers_fasta
        )


    def test_check_file_exists(self):
        self.assertRaises(
            AssertionError, grc.create_fasta_dict, "This is a string"
        )


class TestCreateFASTADict(unittest.TestCase):
    """
    Check that the dictionary is created correctly from a FASTA file
    """
    fasta1 = os.path.join(TEST_DATA_DIR, "BRCA2.fna")
    fasta_dict1 = grc.create_fasta_dict(fasta1)

    fasta2 = os.path.join(TEST_DATA_DIR, "single_DNA.fa")
    fasta_dict2 = grc.create_fasta_dict(fasta2)

    def test_len_fasta_dict(self):
        assert len(self.fasta_dict1) == 2, (
            "Length of dictionary is incorrect with 2 sequences provided"
        )


    def test_len_fasta_dict2(self):
        assert len(self.fasta_dict2) == 1, (
            "Length of dictionary is incorrect with 1 sequence provided"
        )


    def test_fasta_dict_keys(self):
        key_names = [
            'NC_000013.11:32315508-32400268',
            'NC_060937.1:31532753-31617510'
        ]

        assert all(
            keys in self.fasta_dict1.keys() for keys in key_names
        ), "Sequence names not added as keys correctly when two sequences"


    def test_fasta_dict_keys2(self):
        key_name = 'The_sequence_ID'

        assert key_name in self.fasta_dict2.keys(), (
            "Key not correct when 1 sequence entered"
        )


class TestVerifyIsDNA(unittest.TestCase):
    """
    Test that error is raised correctly when sequence is given that is not DNA
    """
    not_DNA_sequence_fasta = os.path.join(TEST_DATA_DIR, "not_DNA.fa")

    def test_DNA_verify(self):
        self.assertRaises(
            SystemExit, grc.verify_is_DNA, self.not_DNA_sequence_fasta, 'record'
        )


class TestReverseComplement(unittest.TestCase):
    """
    Check that the reverse complement is made correctly
    """
    test_sequence1 = 'ACTGATNATCGGGACTA'
    test_sequence2 = 'actgatgtcan-gct'

    def test_reverse_complement1(self):
        assert grc.reverse_comp_sequence(self.test_sequence1) == (
            'TAGTCCCGATNATCAGT'
        ), "Reverse complement not correct with Ns"


    def test_reverse_complement2(self):
        assert grc.reverse_comp_sequence(self.test_sequence2) == (
            'agc-ntgacatcagt'
        ), "Reverse complement not correct with lower case letters"


class TestCreateReverseCompRecords(unittest.TestCase):
    """
    Check that new SeqRecord object with reverse complemenent is made correctly
    """
    test_input_dict = {
        'My_DNA_sequence': SeqRecord(
            seq=Seq('ACTGAGCCA'),
            id='My_DNA_sequence',
            name='My_DNA_sequence',
            description='My_DNA_sequence This is a DNA sequence',
            dbxrefs=[]
        )
    }

    def test_output_description(self):
        output = grc.generate_reverse_comp_records(self.test_input_dict)[0]

        assert output.seq == Seq('TGGCTCAGT'), "Revese complement not correct"

        assert output.id == 'My_DNA_sequence', "ID not kept the same"

        assert output.name == 'My_DNA_sequence', "Name not kept the same"

        assert output.description == (
            'My_DNA_sequence This is a DNA sequence reverse complement'
        ), "Reverse complement not added to sequence description"


if __name__ == '__main__':
    unittest.main()
