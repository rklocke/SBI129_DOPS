# SBI129 Programming: Direct Observation of Practical Skill (DOPS)

## DOPS 1
### Task
*Write a program to create a reverse complement from a DNA sequence, and explain the development environment, the data source, the libraries used and the algorithm employed.*

### Development environment
Code was developed, run and tested on the Ubuntu 22.04.3 LTS operating system. Code was written within Visual Studio Code v1.83.0 and the biopython package was installed in a Conda environment.

### Data source
A FASTA file containing the DNA sequence of the BRCA2 gene (*'BRCA2.fna'*) was obtained from the National Center for Biotechnology Information (NCBI) Gene database: https://www.ncbi.nlm.nih.gov/gene/675

### Libraries used (requirements.txt)
- The [argparse](https://docs.python.org/3/library/argparse.html) library v1.1 was used to specify command line arguments for the input FASTA file and the name of the output FASTA file.
- The [Biopython](https://biopython.org/) library v1.81 was used to parse the FASTA file and generate the reverse complement.
- The [pytest](https://docs.pytest.org/en/7.4.x/) library v7.3.1 was used to assess performance of the code via unit testing.

### Algorithm deployed
- Read the FASTA file (specified as a command line argument) into a dictionary using the Biopython `SeqIO.parse()` and `SeqIO.to_dict()` functions.
    - The dictionary contains the name of each sequence as keys and the values hold information about the sequence as a `SeqRecord` object; the sequence ID, name, description, annotation and the sequence itself.
- Perform checks to ensure the sequence is a DNA sequence and the file formatted correctly.
- Loop over the `SeqRecord` objects and use the Biopython `reverse_complement()` function to obtain the reverse complement of the sequence.
- Create a new `SeqRecord` object for each sequence, keeping the name, ID and annotations the same but adding `reverse complement` to the description in the header of the sequence. Append each object to a list.
- Write the list of new SeqRecord objects to a new FASTA file using `SeqIO.write`.

### Usage
The script takes a FASTA file (which can be a multi-FASTA file with multiple sequences) and outputs a FASTA file containing the reverse complement of each sequence. The tool can be run as follows:

```
usage: generate_rev_comp.py {ARGUMENTS}

required arguments:
 --input_fasta/-i [file]       Path to the input FASTA file
 --output_fasta/-o [file]      Name of the output FASTA file
```

## DOPS 3

*Write a set of unit tests to evaluate whether the output of a script or set of functions is correct, including edge cases.*

The `tests` directory contains test data and the `test_generate_rev_comp.py` script, which includes a selection of unit tests to test functions defined within `generate_rev_comp.py`. It can be run with the `pytest` command from within the parent directory.