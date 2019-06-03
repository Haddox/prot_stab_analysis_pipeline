"""
This script contains `Python` functions for analyzing deep-sequencing data

Hugh Haddox, February-19-2018
"""

# Import `Python` modules
import pandas
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import subprocess
import re
import gzip

def read_fastq(filename, five_prime_flanking_seq, three_prime_flanking_seq, seq_type='protein'):
    """
    Compute the number of DNA or protein counts for each sequence in an input FASTQ file.

    First, this function parses a DNA coding sequence from the sequencing read using input flanking
    sequences. If protein counts are being recorded, it then translates the DNA coding sequence into
    the corresponding protein sequence. It then adds a count for the given sequence to a dictionary
    Finally, it returns the name of the input file and a dictionary with the counts. The reason that
    the function returns the input file name along with the counts is to make it easier to use this
    function to parallelize tasks. I adapted this function from code from Gabe Rocklin and Franziska
    Seeger.

    Args:
        `filename`: the name of the FASTQ file
        `five_prime_flanking_seq`: a DNA sequence that flanks the coding sequence of interest on the
            5' end of the sequencing read (string). The coding sequence should begin immediately after
            the final nucleotide of this flanking sequence. This sequence and the sequence given by
            `three_prime_flanking_seq` will be used to extract the DNA coding sequence from the
            sequencing read. This DNA coding sequence or corresponding protein sequence will be counted
            in the counts dictionary.
        `three_prime_flanking_seq`: a DNA sequence that flanks the coding sequence of interest on the
            3' end of the sequencing read (string; default: 'CTCGAG'). The coding sequence should begin
            immediately before the first nucleotide of this flanking sequence.Note: this DNA sequence should
            be in the same 5'-to-3' orientation as `five_prime_flanking_seq`. The default is set to 'CTCGAG',
            which is the flanking sequence used in Rocklin et al., 2017, Science.
        `seq_type`: a string specifying whether counts should be recorded at the level of protein or DNA.
            The options are: 'protein' or 'dna'

    Returns:
        A tupple of (`filename`, `counts_dict`):
            `filename` : same as input
            `counts_dict` : a dictionary with counts for each unique protein sequence,
                of the form {protein sequence}:{counts}
    """

    # Go through each of the FASTQ file, parse sequences of interest, and add them to a dictionary
    # that tracks the counts of each
    counts_dict = {}
    with open(filename) as file:
        for (i, line) in enumerate(file):

            # Each entry is a set of four lines. The below line of code skips over the header
            # line, i.e., the first one in each set
            if i % 4 != 1: continue

            # ... it also skips over lines that are less than 12 characters
            if len(line) < 12: continue

            # Look for the pattern given by `five_prime_flanking_seq` in the line in the forward
            # direction, find it's position, and incriment that position by the length of this
            # pattern, setting this equal to the starting point of the sequence of interest.
            # Otherwise, if the pattern cannot be found, skip this entry. Note: find always
            # identifies the first position where the pattern occurs in the string
            if five_prime_flanking_seq in line:
                startpt = line.find(five_prime_flanking_seq) + len(five_prime_flanking_seq)
            else:
                continue

            # Look for the pattern given by `three_prime_flanking_seq` in the line starting from
            # the end of the line and find it's position, setting this equal to the ending point
            # of the sequence of interest. Then, if the DNA sequence is a multiple of three,
            # record a count for that sequence in the dictionary, first translating the DNA
            # sequence to a protein sequence if specified by `seq_type`
            if three_prime_flanking_seq in line:
                endpt = line.rfind(three_prime_flanking_seq)
                coding_dna = line[startpt:endpt]
                if (len(coding_dna) % 3 != 0):
                    continue
                if seq_type == 'protein':
                    coding_dna = Seq(coding_dna, generic_dna)
                    seq = str(coding_dna.translate(to_stop=True))
                elif seq_type == 'dna':
                    seq = coding_dna
                else:
                    raise ValueError("Could not parse `seq_type`: {0}".format(seq_type))
                if seq == '':
                    continue
                if seq not in counts_dict:
                    counts_dict[seq] = 0
                counts_dict[seq] += 1

    return (filename, counts_dict)


def read_fastq_with_c_script(filename, c_script_path):
    """
    Performs the same thing as the above functions for parsing FASTQ files, but using a C script written by Luki Goldschmidt

    We are interested in seeing if this C-based program is much faster than the Python based ones above.

    The `-p` flag tells the program to return counts for amino-acid sequences instead of DNA sequences

    Args:
        `filename`: same as the function `read_fastq`
        `c_script_path`: the path to the C script called `parse_fastq`

    Returns: same as the function `read_fastq`
    """

    # Put together the commandline argument to call the C script
    cmd = ' '.join([
        c_script_path,
        '-f CATATG',
        '-r CTCGAG',
        '-p',
        filename
    ])

    # Carry out the argument and get the output, converting it from byte form to a string parasable by Python
    output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    output = output.decode("utf-8")

    # Parse counts from the output, which has one newline-delimited entry per sequence, each of which has two tab-
    # delimited columns of the sequence followed by its counts. I don't consider the last newline-delimited entry
    # since it is a blank line
    counts_dict = dict([(line.split('\t')[0], int(line.split('\t')[1])) for line in output.split('\n')[:-1]])

    return (filename, counts_dict)


def ComputeCounts(fastq_files_and_output_file):
    """
    Compute DNA or protein counts from a list of FASTQ files and write the counts to a file

    This function computes DNA or protein counts from a list of FASTQ files, aggregating the results
    among all FASTQ files. I use this function to parallelize this process. I combine multiple
    variables in a single input tupple, which makes the parallelization command much simpler.

    Args:
        `fastq_files_and_outputprefix`: A tupple of the following variables:
            `fastq_files: A list of paths to input FASTQ files.
            `outputfile`: The name the output CSV file reporting sequence counts.
            `five_prime_flanking_seq`: Flanking sequence used to parse sequencing reads. See
                the `read_fastq` function for more details.
            `three_prime_flanking_seq`: Flanking sequence used to parse sequencing reads. See
                the `read_fastq` function for more details.
            `seq_type`: a string specifying whether counts should be recorded at the level of\
                protein or DNA. The options are: 'protein' or 'dna'.

    Returns:
        A CSV of sequence counts. The first line of this file is a header:

            `sequence,counts`

        while the remaining lines are comma-delimited entries of these variables:

            `{sequence},{counts}`
    """

    # Unpack the input
    (fastq_files, output_file, five_prime_flanking_seq, three_prime_flanking_seq, seq_type) = fastq_files_and_output_file

    # Compute counts for each FASTQ file
    assert len(fastq_files) > 0, "The list of FASTQ files is empty"
    counts_per_file = {}
    for (i, fastq_file) in enumerate(fastq_files):
        counts_per_file[i] = read_fastq(
            filename = fastq_file,
            five_prime_flanking_seq = five_prime_flanking_seq,
            three_prime_flanking_seq = three_prime_flanking_seq,
            seq_type = seq_type
        )[1] # return the dictionary of counts

    # Aggregate counts acorss all FASTQ files
    if len(fastq_files) == 1:
        assert len(counts_per_file.keys()) == 1
        counts_dict = counts_per_file[0]
    else:
        counts_dict = {}
        for i in counts_per_file.keys():
            for seq in counts_per_file[i]:
                if seq not in counts_dict.keys():
                    counts_dict[seq] = counts_per_file[i][seq]
                else:
                    counts_dict[seq] += counts_per_file[i][seq]

    # Write the results to an output file
    counts_df = pandas.DataFrame(list(counts_dict.items()), columns=['sequence', 'counts'])
    counts_df.to_csv(output_file, index=False)


def ParsePAREOutfile(outfile):
    """
    This function parses the output data generated by PARE when assemblying paired-end reads

    Args:
        `outfile`: the path to a file with the output data generated by PARE

    Returns:
        A tupple with the following three variables in the order they appear in the below list:
            `n_assembled_reads` : the total number of assembled reads
            `n_discarded_reads` : the total number of discarded reads
            `n_non_assembled_reads` : the total number of non_assembled_reads
    """

    # Pattern used to extract data
    n_reads_pattern = re.compile(r'\: (?P<n_reads>[\d,]+) /')

    # Use regular expressions to extract the relevant info from the file
    n_assembled_reads = n_discarded_reads = n_non_assembled_reads = n_total_reads = None
    with open(outfile) as f:
        #print(f.readlines())
        for line in f:
            if 'Assembled reads .' in line:
                if n_assembled_reads:
                    raise ValueError("Already found data for `n_assembled_reads`")
                n_assembled_match = re.search(n_reads_pattern, line)
                n_assembled_reads = int(n_assembled_match.group('n_reads').replace(',', ''))

            elif 'Discarded reads .' in line:
                if n_discarded_reads:
                    raise ValueError("Already found data for `n_discarded_reads`")
                n_discarded_match = re.search(n_reads_pattern, line)
                n_discarded_reads = int(n_discarded_match.group('n_reads').replace(',', ''))

            elif 'Not assembled reads .' in line:
                if n_non_assembled_reads:
                    raise ValueError("Already found data for `n_non_assembled_reads`")
                n_non_assembled_match = re.search(n_reads_pattern, line)
                n_non_assembled_reads = int(n_non_assembled_match.group('n_reads').replace(',', ''))

    return (n_assembled_reads, n_discarded_reads, n_non_assembled_reads)
