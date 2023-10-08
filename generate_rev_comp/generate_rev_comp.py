#!/usr/bin/env python3

"""
SBI129: Programming
DOPS 1: Write a program to create a reverse complement from a DNA sequence,
and explain the development environment, the data source, the libraries used
and the algorithm employed.
"""

import argparse
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information necessary to generate reverse complement'
    )

    parser.add_argument(
        '-i',
        '--input_fasta',
        type=str,
        required=True,
        help='Name of input FASTA containing DNA sequence(s)'
    )

    parser.add_argument(
        '-o',
        '--output_fasta',
        type=str,
        required=True,
        help='Name of output FASTA containing the reverse complement'
    )

    args = parser.parse_args()

    return args


def create_fasta_dict(fasta_file):
    """
    Read in the FASTA file to a dictionary where keys are sequence ID
    and the values are SeqRecords

    Parameters
    ----------
    fasta_file : file
        path to FASTA file to read into a dictionary

    Returns
    -------
    record_dict : dict
        dict with records for each of the sequences present in the FASTA file
    """
    assert os.path.isfile(fasta_file), "Sorry, your input is not a file :("

    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    if not record_dict:
        sys.exit(
            "Error: No sequences found in input FASTA file. Please check your"
            " input file is formatted correctly with a header line for each "
            "sequence"
        )

    return record_dict


def verify_is_DNA(sequence, record):
    """
    Verify that the input sequence is DNA and doesn't have any bad characters

    Parameters
    ----------
    sequence : str
        the sequence to have the reverse complement generated
    record: str
        ID of the sequence

    """
    if sequence:
        allowed_chars = 'ACTGN-'

        # Find set of input characters
        unique_chars = set((sequence).upper())

        # Get any characters which are not allowed
        bad_characters = unique_chars - set(allowed_chars)

        if not unique_chars.issubset(allowed_chars):
            sys.exit(
                "Error: Please check your input FASTA is DNA. The sequence "
                f"for record {record} contains the following characters which "
                f"are not valid: {bad_characters}"
            )
    else:
        print(f"Warning: empty sequence found for record: {record}")


def reverse_comp_sequence(sequence):
    """
    Create the reverse complement of a sequence

    Parameters
    ----------
    sequence : str
        the sequence of letters to generate reverse complement for

    Returns
    -------
    reverse_comp : Bio.Seq.Seq object
        the reverse complement of the sequence
    """
    reverse_comp = Seq(sequence).reverse_complement()

    return reverse_comp


def generate_reverse_comp_records(record_dict):
    """
    Create a list of SeqRecords the same as the input but with the reverse
    complement sequence

    Parameters
    ----------
    record_dict : dict
        dict with records for each of the sequences present in the FASTA file

    Returns
    -------
    reverse_records : list
        list of SeqRecord object(s) with reverse complement sequences
    """
    reverse_records = []

    for record in record_dict:
        # Get sequence as a stringm, verify it's DNA and then get reverse
        # complement
        sequence = str(record_dict[record].seq)
        verify_is_DNA(sequence, record)
        rev_comp_seq = reverse_comp_sequence(sequence)

        # Create a new SeqRecord object for the sequence and its attributes
        # adding 'reverse complement' to the sequence header
        new_record = SeqRecord(
            Seq(rev_comp_seq),
            id=record_dict[record].id,
            name=record_dict[record].name,
            description=record_dict[record].description+" reverse complement",
            annotations=record_dict[record].annotations
        )

        reverse_records.append(new_record)

    return reverse_records


def write_out_reverse_complement_fasta(reverse_records, outfile_name):
    """
    Write out the list of SeqRecord objects to a new FASTA file

    Parameters
    ----------
    reverse_records : list
        list of SeqRecord objects
    outfile_name : str
        name of the output file with reverse complement sequences
    """
    SeqIO.write(reverse_records, outfile_name, "fasta")


def main():
    """
    Main function to generate FASTA file with reverse complement sequences
    """
    args = parse_args()
    record_dict = create_fasta_dict(args.input_fasta)
    reverse_record_dict = generate_reverse_comp_records(record_dict)
    write_out_reverse_complement_fasta(reverse_record_dict, args.output_fasta)


if __name__ == '__main__':
    main()
