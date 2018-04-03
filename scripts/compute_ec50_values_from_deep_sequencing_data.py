"""
Compute EC50 values using experimental data from the high-throughput protein-stability assay
"""

# Import `Python` modules
import os
import argparse
import glob

import pandas


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
    parser.add_argument("--output_dir", help="a path to an output directory where all the results will be stored. This directory will be created if it does not already exist")
    args = parser.parse_args()

    # Assign command-line arguments to variables
    designed_sequences_file = args.designed_sequences_file
    experimental_summary_file = args.experimental_summary_file
    fastq_dir = args.fastq_dir
    pear_path = args.pear_path
    output_dir = args.output_dir
    
    # Initialize output directories if they do not already exist
    paired_FASTQ_files_dir = os.path.join(output_dir, 'paired_FASTQ_files')
    counts_dir = os.path.join(output_dir, 'counts')
    dirs = [output_dir, paired_FASTQ_files_dir, counts_dir]
    for dir_i in dirs:
        if not os.path.isdir(dir_i):
            print("Making the directory: {0}".format(dir_i))
            os.makedirs(dir_i)

    # Read the designed sequences into a dataframe
    designed_seqs_df = pandas.read_csv(designed_sequences_file)
    designed_seqs_df.set_index('protein_sequence', inplace=True)

    # Read the data from the experiments summary file into a dataframe
    # and get a list of unique proteases and unqiue selection indices
    summary_df = pandas.read_csv(experimental_summary_file)
    proteases = set(summary_df['protease_type'])
    selection_indices = set(summary_df['selection_strength'])

    # Iterate through each sample in the experimental summary dataframe,
    # make a list of paired-end FASTQ files for each sample, and use `PARE`
    # to assemble each file pair
    FASTQ_files = {}
    for (i, row) in summary_df.iterrows():

        # Get sample metadata
        protease_type = row['protease_type']
        selection_index = row['selection_strength']
        experiment_name = '{0}_{1}'.format(protease_type, selection_index)
        fastq_id = row['fastq_id'].replace('_', '-')

        # Find all R1 and R2 files
        r1_files = glob.glob('{0}/{1}*_R1_*.fastq*'.format(fastq_dir, fastq_id))
        r2_files = [f.replace('_R1_', '_R2_') for f in r1_files]

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
        
        # Assemble each pair of FASTQ files using `PARE`
        assert experiment_name not in FASTQ_files.keys(), \
            "Duplicate experiment name: {0}".format(experiment_name)
        FASTQ_files[experiment_name] = []
        for (pair_n, (r1_file, r2_file)) in enumerate(zip(r1_files, r2_files), 1):

            # Make a prefix for the outfiles, including a file with assembled reads
            # and a log file for the `PEAR` output
            outfile_prefix = os.path.join(paired_FASTQ_files_dir, '{0}-{1}'.format(experiment_name, pair_n))
            logfile = '{0}.log'.format(outfile_prefix)
            paired_output_file = '{0}.fastq'.format(outfile_prefix)
            
            # Don't rerun the assembly if the ouptut file already exists. Otherwise,
            # run the assembly
            if os.path.isfile(paired_output_file):
                print("The paired output FASTQ file already exists. Will not rerun `PARE`.")
                continue
            else:
                print("\nUsing PARE to assembling the files {0} and {1}".format(r1_file, r2_file))
                print("Storing the results in an assembled FASTQ file called {0}".format(paired_output_file))
                print("Logging the results of the run in a file called {0}".format(logfile))
                cmd = ' '.join([
                    pear_path,
                    '-f {0}'.format(r1_file),
                    '-r {0}'.format(r2_file),
                    '-o {0}'.format(outfile_prefix),
                    '-j 20',
                    '> {0}'.format(logfile)
                ])
                #{cmd}
            break
        
    # For each sample, compute protein-level counts from the deep-sequencing data. Write
    # a file giving all observed counts, even for proteins that don't match a starting
    # design.
        
    
    # For each protease, aggregate counts across all samples, only including proteins
    # that are within the set of input designs, throwing out mismatching proteins
    
    
    # Quantify the fraction of deep-sequencing counts for proteins that match one of
    # the starting designs or controls, and write the results to a file for computing
    # EC50 values
    
    
    
    

if __name__ == "__main__":
    main()
    