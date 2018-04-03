"""
Compute EC50 values from experimental data from the high-throughput protein-stability assay
"""

# Import `Python` modules
import os
import sys
import pandas
import argparse

def main():
    """Read in command-line arguments and execute the main code"""
    
    # Get a path to the directory of this script relative to the current
    # working directory of where it is being called from
    scriptsdir = os.path.dirname(__file__)

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--designed_sequences_file", help="a file with design names and protein sequences")
    parser.add_argument("--experimental_summary_file", help="a file with experimental metadata")
    parser.add_argument("--fastq_dir", help="the path to a directory with input FASTQ files")
    parser.add_argument("--pear_path", help="the path to the program `PEAR`")
    args = parser.parse_args()

    # Assign command-line arguments to variables
    designed_sequences_file = args.designed_sequences_file
    experimental_summary_file = args.experimental_summary_file
    fastq_dir = args.fastq_dir
    pear_path = args.pear_path

    # Read in the designed sequences
    designed_seqs_df = pandas.read_csv(designed_sequences_file)
    designed_seqs_df.set_index('protein_sequence', inplace=True)

    # Read in data from the experiments summary file and get a list of
    # unique proteases and unique selection indices
    summary_df = pandas.read_csv(experimental_summary_file)
    proteases = set(summary_df['protease_type'])
    selection_indices = set(summary_df['selection_strength'])
    print(proteases)
    print(selection_indices)
    
    
    
    
if __name__ == "__main__":
    main()
    

"""
# Read in data from the experiments summary file
proteases = ['trypsin', 'chymotrypsin']
summary_df = pandas.read_csv(experimental_summary_file)
# selection_indices =  # get these from the experimental summary file

# Read in the designed sequences
designed_seqs_df = pandas.read_csv(designed_sequences_file, sep='\s+')
designed_seqs_df.set_index('protein_sequence', inplace=True)

# Iterate through each sample in the experimental summary dataframe and compile
# a list of FASTQ files for each sample
FASTQ_files = {}
for (i, row) in summary_df.iterrows():

    # Get sample metadata
    protease_type = row['protease type']
    selection_index = row['selection_index']
    experiment_name = '{0}_{1}'.format(protease_type, selection_index)
    fastq_id = row['experiment/round'].replace('_', '-')

    # Find R1 and R2 files and append them to a list
    r1_files = glob.glob('{0}/{1}*_R1_*.fastq*'.format(fastq_dir, fastq_id))
    r2_files = [f.replace('_R1_', '_R2_') for f in r1_files]
    assert experiment_name not in FASTQ_files.keys(), \
        "Duplicate experiment name: {0}".format(experiment_name)
    FASTQ_files[experiment_name] = list(zip(r1_files, r2_files))

    # Make sure that each sample has the same number of R1 and R2 files, that
    # there are more than one of each, and that the patterns "_R1_" and "_R2_"
    # don't appear more than once
    assert(len(r1_files) == len(r2_files))
    if len(r1_files) == 0:
        raise ValueError(
            "Failed to find FASTQ files for the fastq_id: {0}".format(fastq_ID)
        )
    for (f1, f2) in zip(r1_files, r2_files):
        assert f1.count('_R1_') == f2.count('_R2_') == 1, \
            "The string '_R1_' or '_R2_' appear multiple times in file name"
        if not os.path.isfile(f2):
            raise ValueError(
                "Failed to find a matching R2 file for: {0}".format(f)
            )
"""