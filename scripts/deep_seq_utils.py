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

def original_read_fastq(fn):
    """
    Original function to compute the number of counts for each unique protein sequence in an input FASTQ file

    Written by Gabe Rocklin and Franziska Seeger, this script parses protein sequences from deep-sequencing reads
    based on defined flanking sequences and counts the number of occurrences of each unique sequence. I made a
    few small tweaks to make it generalizable (e.g., removed part about file indexing)

    Args:
        `fn` : the name of the FASTQ file

    Returns:
        A list of (`fn`, `counts_dict`):
            `fn` : same as input
            `sequences_found` : a dictionary with counts for each unique protein sequence,
                of the form {protein sequence}:{counts}
    """

    sequences_found = {}
    with open(fn) as file:
        lines = [x.strip() for x in file.readlines()]

    totalcounts = 0
    for i in range(len(lines)):
        if i % 4 != 1: continue
        if len(lines[i]) < 12: continue
        if 'CATATG' in lines[i]:
            startpt = lines[i].find('CATATG')+6
        else:
            continue
        if 'CTCGAG' in lines[i]:

            endpt = lines[i].rfind('CTCGAG')

            coding_dna = Seq(lines[i][startpt:endpt], generic_dna)
            protein_seq = str(coding_dna.translate(to_stop=True))
            if protein_seq == '': continue
            if protein_seq not in sequences_found:
                sequences_found[protein_seq] = 0
            sequences_found[protein_seq] += 1
    return [fn, sequences_found]


def read_fastq(filename):
    """
    Compute the number of counts for each unique protein sequence in an input FASTQ file

    This function is the same as `original_read_fastq`, but I changed it so that it is somewhat
    more memory efficient. Specifically, instead of reading the whole FASTQ file into memory, it
    iterates over it line by line.

    Args:
        `filename` : the name of the FASTQ file

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

            # Look for the pattern "CATATG" in the line in the forward direction, find it's
            # position, and incriment that position by six, setting this equal to the starting
            # point of the sequence of interest. Otherwise, if the pattern cannot be found, skip
            # this entry. Note: find always identifies the first position where the pattern
            # occurs in the string
            if 'CATATG' in line:
                startpt = line.find('CATATG')+6
            else:
                continue

            # Look for the pattern "CTCGAG" in the line in the reverse direction and find it's
            # position, setting this equal to the ending point of the sequence of interest, then
            # translate the DNA sequence into a protein sequence
            if 'CTCGAG' in line:
                endpt = line.rfind('CTCGAG')
                coding_dna = Seq(line[startpt:endpt], generic_dna)
                protein_seq = str(coding_dna.translate(to_stop=True))
                if protein_seq == '':
                    continue
                if protein_seq not in counts_dict:
                    counts_dict[protein_seq] = 0
                counts_dict[protein_seq] += 1

    return (filename, counts_dict)

def read_fastq_only_forward_read(filename, length_cds_to_parse=90):
    """
    Compute the number of counts for each unique protein sequence in an input
    FASTQ file with a single read that doesn't cover the whole coding sequence

    This function is the same as `read_fastq`, but I changed it so that it is
    able to analyze incomplete deep-sequencing data, specifically related to
    a scenario where you have a single forward read that doesn't cover the entire
    coding sequence.

    Args:
        `filename`: the name of the FASTQ file
        `length_cds_to_parse`: the number of nucleotides of the coding sequence
            to parse from each deep-sequencing read, where the start of the
            coding sequence is definied by the 5' flanking restriction cut site,
            the sequence of which is hard coded below. (default: 102)

    Returns:
        A tupple of (`filename`, `counts_dict`):
            `filename` : same as input
            `counts_dict` : a dictionary with counts for each unique protein sequence,
                of the form {protein sequence}:{counts}
    """

    # Go through each of the FASTQ file, parse sequences of interest, and add them to a dictionary
    # that tracks the counts of each
    counts_dict = {}
    
    with gzip.open(filename, 'rb') as file:
        for (i, line) in enumerate(file):

            # Each entry is a set of four lines. The below line of code skips over the header
            # line, i.e., the first one in each set
            if i % 4 != 1: continue

            # ... it also skips over lines that are less than 12 characters
            if len(line) < 12: continue

            # Look for the pattern "CATATG" in the line in the forward direction, find it's
            # position, and incriment that position by six, setting this equal to the starting
            # point of the sequence of interest. Otherwise, if the pattern cannot be found, skip
            # this entry. Note: find always identifies the first position where the pattern
            # occurs in the string
            if 'CATATG' in line:
                startpt = line.find('CATATG')+6
            else:
                continue

            endpt = startpt + length_cds_to_parse
            coding_dna = Seq(line[startpt:endpt], generic_dna)
            protein_seq = str(coding_dna.translate(to_stop=True))
            if protein_seq == '':
                continue
            if protein_seq not in counts_dict:
                counts_dict[protein_seq] = 0
            counts_dict[protein_seq] += 1

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
    Compute protein counts from a list of FASTQ files and write the counts to a file

    This function computes protein counts from a list of FASTQ files, aggregating the results
    among all FASTQ files. I use this function to parallelize this process. I combine multiple
    variables in a single input tupple, which makes the parallelization command much simpler.
    Currently, it uses the above function `read_fastq` to compute counts, which is more memory
    efficient than the function originally used in the lab `original_read_fastq`.

    Args:
        `fastq_files_and_outputprefix`: A tupple of the following variables:
            `fastq_files: A list of paths to input FASTQ files.
            `outputfile`: The name the output CSV file reporting sequence counts.
        `forward_read_only`: If the FASTQ data only has a single forward read
            (default: False)

    Returns:
        A CSV of sequence counts. The first line of this file is a header:

            `sequence,counts`

        while the remaining lines are comma-delimited entries of these variables:

            `{sequence},{counts}`
    """

    # Unpack the input
    (fastq_files, output_file) = fastq_files_and_output_file

    # Compute counts for each FASTQ file
    assert len(fastq_files) > 0, "The list of FASTQ files is empty"
    counts_per_file = {}
    for (i, fastq_file) in enumerate(fastq_files):
        counts_per_file[i] = read_fastq(fastq_file)[1] # return the dictionary of counts

    # Aggregate counts acorss all FASTQ files, or just return the one dictionary if there is only one input
    # FASTQ file
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

    
def compute_counts_forward_read_only(fastq_files_and_output_file):
    """
    Same as above, but for scenarios with just forward reads.
    """

    # Unpack the input
    (fastq_files, output_file) = fastq_files_and_output_file

    # Compute counts for each FASTQ file
    assert len(fastq_files) > 0, "The list of FASTQ files is empty"
    counts_per_file = {}
    for (i, fastq_file) in enumerate(fastq_files):
        counts_per_file[i] = read_fastq_only_forward_read(fastq_file)[1] # return the dictionary of counts

    # Aggregate counts acorss all FASTQ files, or just return the one dictionary if there is only one input
    # FASTQ file
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
