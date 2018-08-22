"""
Compute EC50 values using experimental data from the high-throughput protein-stability assay
"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
from multiprocessing import Pool
import time
import pandas

# Custom `Python` modules
import deep_seq_utils

# Get the path to this script
scriptsdir = os.path.dirname(__file__)

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""
    
    #---------------------------------------------------------------
    # Read in command-line arguments and experimental metadata
    #---------------------------------------------------------------
    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--designed_sequences_file", help="a file with design names and protein sequences")
    parser.add_argument("--experimental_summary_file", help="a file with experimental metadata")
    parser.add_argument("--fastq_dir", help="the path to a directory with input FASTQ files")
    parser.add_argument("--pear_path", help="the path to the program `PEAR`")
    parser.add_argument("--five_prime_flanking_seq", help="five-prime flanking sequence used to parse DNA coding sequences from deep-sequencing reads")
    parser.add_argument("--three_prime_flanking_seq", help="three-prime flanking sequence used to parse DNA coding sequences from deep-sequencing reads")
    parser.add_argument("--output_dir", help="a path to an output directory where all the results will be stored. This directory will be created if it does not already exist")
    args = parser.parse_args()    
    
    # Assign command-line arguments to variables
    designed_sequences_file = args.designed_sequences_file
    experimental_summary_file = args.experimental_summary_file
    fastq_dir = args.fastq_dir
    pear_path = args.pear_path
    five_prime_flanking_seq = args.five_prime_flanking_seq
    three_prime_flanking_seq = args.three_prime_flanking_seq
    output_dir = args.output_dir
    print("\nHere is a list of parsed input arguments")
    for arg in vars(args):
        print("{0}: {1}".format(arg, getattr(args, arg)))
    
    # Remove the terminal `/` in the `fastq_dir` arg if it exists
    # as code lower down assumes this `/` doesn't exist
    if fastq_dir[-1] == '/':
        fastq_dir = fastq_dir[:-1]
    
    # Initialize output directories if they do not already exist
    paired_FASTQ_files_dir = os.path.join(output_dir, 'paired_FASTQ_files')
    counts_dir = os.path.join(output_dir, 'counts')
    ec50s_dir = os.path.join(output_dir, 'ec50_values')
    stability_scores_dir = os.path.join(output_dir, 'stability_scores')
    dirs = [
        output_dir, paired_FASTQ_files_dir, counts_dir, ec50s_dir,
        stability_scores_dir
    ]
    for dir_i in dirs:
        if not os.path.isdir(dir_i):
            print("\nMaking the directory: {0}".format(dir_i))
            os.makedirs(dir_i)

    # Read the data from the designed sequences file into a dataframe
    # and make sure that there aren't any duplicate entires in either
    # the `name` or `protein_sequence` columns
    print("\nReading in input design sequences from the file: {0}".format(designed_sequences_file))
    designed_seqs_df = pandas.read_csv(designed_sequences_file)
    all_designed_seqs = designed_seqs_df['protein_sequence']
    all_designed_names = designed_seqs_df['name']
    if len(all_designed_seqs) != len(set(all_designed_seqs)):
        raise ValueError("The `designed_sequences_file` has duplicate sequences. All sequences must be unique.")
    if len(all_designed_names) != len(set(all_designed_names)):
        raise ValueError("The `designed_sequences_file` has duplicate names. All names must be unique.")
    designed_seqs_df.set_index('protein_sequence', inplace=True)
    
    # Read the data from the experiments summary file into a dataframe
    # and get a list of unique proteases and unqiue selection indices
    summary_df = pandas.read_csv(experimental_summary_file)
    proteases = set(summary_df['protease_type'])
    selection_indices = set(summary_df['selection_strength'])
    
    
    #---------------------------------------------------------------
    # Find and assemble all paired-end FASTQ files reads for each sample
    #---------------------------------------------------------------
    
    # Iterate through each sample in the experimental summary dataframe,
    # make a list of paired-end FASTQ files for each sample, and use `PARE`
    # to assemble each file pair
    FASTQ_files = {}
    for (i, row) in summary_df.iterrows():

        # Get sample metadata
        protease_type = row['protease_type']
        selection_index = row['selection_strength']
        experiment_name = '{0}_{1}'.format(protease_type, selection_index)
        fastq_id = row['fastq_id']
        print("\nAggergating paired-end FASTQ files for {0} and selection index {1} using the FASTQ ID {2}".format(
            protease_type, selection_index, fastq_id)
        )
        
        # Find all R1 and R2 files
        r1_files = glob.glob('{0}/{1}*_R1_*.fastq*'.format(fastq_dir, fastq_id))
        r2_files = [f.replace('_R1_', '_R2_') for f in r1_files]

        # Make sure that each sample has the same number of R1 and R2 files, that
        # there are more than one of each, and that the patterns "_R1_" and "_R2_"
        # don't appear more than once
        assert(len(r1_files) == len(r2_files))
        if len(r1_files) == 0:
            raise ValueError(
                "Failed to find FASTQ files for the fastq_id: {0}".format(fastq_id)
            )
        for (f1, f2) in zip(r1_files, r2_files):
            assert f1.count('_R1_') == f2.count('_R2_') == 1, \
                "The string '_R1_' or '_R2_' appear multiple times in file name"
            if not os.path.isfile(f2):
                raise ValueError(
                    "Failed to find a matching R2 file for: {0}".format(f)
                )
        print("Here is a list of R1 files: {0}".format(', '.join(r1_files)))
        print("Here is a list of R2 files: {0}".format(', '.join(r1_files)))
        
        # Assemble each pair of FASTQ files using `PARE`, storing the assembled files
        # in a dictionary called `FASTQ_files`, which has the following format:
        # {`experiment_name`} : {list of all assembled FASTQ files for a given sample},
        # where `experiment_name` is a sample-specific variable defined above.
        assert experiment_name not in FASTQ_files.keys(), \
            "Duplicate experiment name: {0}".format(experiment_name)
        FASTQ_files[experiment_name] = []
        for (pair_n, (r1_file, r2_file)) in enumerate(zip(r1_files, r2_files), 1):

            # Make a prefix for the outfiles, including a file with assembled reads
            # and a log file for the `PEAR` output
            outfile_prefix = os.path.join(paired_FASTQ_files_dir, '{0}-{1}'.format(experiment_name, pair_n))
            logfile = '{0}.log'.format(outfile_prefix)
            paired_FASTQ_output_file = '{0}.assembled.fastq'.format(outfile_prefix)
            FASTQ_files[experiment_name].append(paired_FASTQ_output_file)
            print("\nPreparing to assemble the files {0} and {1}, outputting the assembled reads to a file called {2}".format(
                r1_file, r2_file, paired_FASTQ_output_file
            ))
            
            # Don't rerun the assembly if the ouptut file already exists. Otherwise,
            # run the assembly
            if os.path.isfile(paired_FASTQ_output_file):
                print("\nThe paired output FASTQ file called {0} already exists. Will not rerun `PARE`.".format(
                    paired_FASTQ_output_file
                ))
                break
                continue
            else:
                cmd = [
                    pear_path,
                    '-f', r1_file,
                    '-r', r2_file,
                    '-o', outfile_prefix,
                    '-j', '20'
                ]
                print("Calling PARE with the command: {0}".format(' '.join(cmd)))
                print("Storing the results in an assembled FASTQ file called: {0}".format(paired_FASTQ_output_file))
                print("Logging the results of the run in a file called: {0}".format(logfile))
                log = subprocess.Popen(cmd, stdout=subprocess.PIPE)
                out, err = log.communicate()
                with open(logfile, 'wb') as f:
                    f.write(out)
           
    
    #---------------------------------------------------------------
    # Compute protein counts from the deep-sequencing data for each sample
    #---------------------------------------------------------------
    # For each sample, compute protein-level counts from the deep-sequencing data. Write
    # a file giving all observed counts, even for proteins that don't match a starting
    # design.
    print("\nComputing protein counts from the deep-sequencing data")
    print("Using the following 5' and 3' flanking sequences to parse coding sequences from sequencing reads:")
    print("5' flanking sequence: {0}".format(five_prime_flanking_seq))
    print("3' flanking sequence: {0}".format(three_prime_flanking_seq))
    input_data_for_computing_counts = []
    for experiment_name in FASTQ_files:
        fastq_files = FASTQ_files[experiment_name]
        output_counts_file = os.path.join(counts_dir, '{0}_counts.csv'.format(experiment_name))
        if os.path.isfile(output_counts_file):
            print("Already have counts for the sample: {0}".format(experiment_name))
            continue
        else:
            print("Computing counts for the sample: {0}".format(experiment_name))
            input_data_for_computing_counts.append(
                (fastq_files, output_counts_file, five_prime_flanking_seq, three_prime_flanking_seq)
            )

    # Compute the counts in parallel, reporting the progress of the computation
    if input_data_for_computing_counts:
        n_procs = len(input_data_for_computing_counts)
        myPool = Pool(processes=n_procs)
        start = time.time()
        all_out = myPool.map_async(deep_seq_utils.ComputeCounts, input_data_for_computing_counts)
        while not all_out.ready():
            print("%s tasks remaining after %s seconds." % (str(all_out._number_left), round(time.time() - start, 0)))
            time.sleep(300.0)
        print("Finished %s processes in %s seconds!" % (n_procs, round(time.time() - start, 2)))

        
    #---------------------------------------------------------------
    # Create input counts and metadata files for computing EC50 values
    #---------------------------------------------------------------
    # The input files include a single `experiments.csv` file with experimental metadata
    # as well as one file per protease that reports counts at each selection step for
    # each protein that matches one of the input designs. To make these input files, I will
    # first initiate a dictionary for keeping track of information for writing the
    # `experiments.csv` file
    counts_metadata_dict = {
        key : []
        for key in ['protease_type', 'selection_strength', 'input', 'column', 'matching_sequences']
    }
    
    # Next, for each protease, I will aggregate the counts and record metadata
    for protease in proteases:
        
        print("\nAggregating counts across selection steps for {0}".format(protease))
        aggregate_counts_outfile = os.path.join(counts_dir, '{0}.counts'.format(protease))
        aggregate_df = designed_seqs_df.copy(deep=True)
        
        for selection_index in selection_indices:
            
            # Read in the counts file
            counts_file = os.path.join(counts_dir, '{0}_{1}_counts.csv'.format(protease, selection_index))
            assert os.path.isfile(counts_file), "Could not find the file: {0}".format(counts_file)
            df = pandas.read_csv(counts_file)
            df.set_index('sequence', inplace=True)
            
            # Next, merge the counts with the dataframe of input sequences, using
            # `how="left"` to only merge on the rows in the dataframe of input sequences.
            # Some of the original sequences may not have counts if they weren't observed,
            # which will lead to values of nan. I will convert these values to counts of
            # zero.
            aggregate_df = aggregate_df.merge(
                df, left_index=True, right_index=True, how="left"
            )
            aggregate_df.fillna(value=0, inplace=True)
            aggregate_df.rename(
                index=str,
                columns={'counts': 'counts{0}'.format(selection_index)},
                inplace=True
            )
                        
            # Next, I will compute the total number of deep-sequencing counts, including
            # for sequences that do not match one of the input sequences
            total_counts = sum(df['counts'])
            n_unique_sequences = len(df['counts'].index.values)

            # Next, I will compute the fraction of all counts that actually match one of
            # the input sequences. To do so, I will first merge the counts dataframe
            # with a dataframe from above with sequences of input designs, using
            # `how=outer` to use the union of indices from both frames, so that no row
            # gets dropped.
            df = df.merge(designed_seqs_df, left_index=True, right_index=True, how="outer")
            total_counts_matching_starting_seqs = sum(
                df[
                    (df['counts'].notnull()) & (df['name'].notnull())
                ]['counts']
            )
            matching_sequences = round(total_counts_matching_starting_seqs / total_counts, 3)

            # Add counts metadata
            counts_metadata_dict['protease_type'].append(protease)
            counts_metadata_dict['selection_strength'].append(selection_index)
            counts_metadata_dict['input'].append('{0}.counts'.format(protease))
            counts_metadata_dict['column'].append('counts{0}'.format(selection_index))
            counts_metadata_dict['matching_sequences'].append(matching_sequences)
            
        # Write the aggregated counts to an output file, only including sequences with
        # at least one count across all experimental samples that were sequenced
        aggregate_df.set_index('name', inplace=True)
        columns_to_write = [
            'counts{0}'.format(selection_index)
            for selection_index in selection_indices
        ]
        aggregate_df['sum_all_counts'] = aggregate_df.apply(
            lambda row: sum([row[counts_col] for counts_col in columns_to_write]), axis=1
        )
        sequences_without_any_counts_bools = aggregate_df['sum_all_counts'] == 0
        sequences_with_counts_bools = aggregate_df['sum_all_counts'] > 0
        assert sum(sequences_without_any_counts_bools) + sum(sequences_with_counts_bools) == len(aggregate_df)
        print("Found {0} sequences that did not have any counts. These sequences will not be written to the output file".format(
            sum(sequences_without_any_counts_bools)
        ))
        print("Found {0} sequences that did have counts".format(
            sum(sequences_with_counts_bools)
        ))
        print("Writing the aggregated counts to the file {0}".format(aggregate_counts_outfile))
        aggregate_df[sequences_with_counts_bools][columns_to_write].to_csv(aggregate_counts_outfile, sep=' ')

    # Make a dataframe indexed in the same way as the `experiments.csv` file by merging
    # the counts metadata with the input file with metadata
    counts_metadata_df = pandas.DataFrame.from_dict(counts_metadata_dict)
    counts_metadata_df.set_index(['protease_type', 'selection_strength'], inplace=True)
    summary_df.dropna(how='all', inplace=True)
    summary_df.set_index(['protease_type', 'selection_strength'], inplace=True)
    assert(
        sorted(list(counts_metadata_df.index.values)) == \
            sorted(list(summary_df.index.values))
    ), "Error in processing the sample-specific metadata"
    summary_df = summary_df.merge(
        counts_metadata_df, left_index=True, right_index=True, how="outer"
    )
    summary_df.reset_index(inplace=True)
    experiments_column_order = [
        'input', 'column', 'parent', 'selection_strength', 'conc_factor', 'parent_expression',
        'fraction_collected', 'matching_sequences', 'cells_collected'
    ]
    output_experiments_file = os.path.join(ec50s_dir, 'experiments.csv')
    summary_df[experiments_column_order].to_csv(output_experiments_file, index=False)


    #---------------------------------------------------------------
    # Compute EC50 values from the deep-sequencing counts using the
    # script `fit_all_ec50_data.py`
    #---------------------------------------------------------------    
    ec50_logfile = os.path.join(ec50s_dir, 'fit_all_ec50_data.log')
    ec50_errfile = os.path.join(ec50s_dir, 'fit_all_ec50_data.err')
    if os.path.isfile(ec50_logfile):
        print("\nEC50 values already exist. To rerun the computation, remove the logfile called: {0}".format(ec50_logfile))
    else:
        cmd = [
            'python',
            '{0}/fit_all_ec50_data.py'.format(scriptsdir),
            '--counts_dir', counts_dir,
            '--experimental_summary_file', output_experiments_file,
            '--datasets', (','.join(proteases)),
            '--output_dir', ec50s_dir
        ]
        print("\nComputing EC50 values with the command: {0}".format(
            ' '.join(cmd)
        ))
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out, err = process.communicate()
        with open(ec50_logfile, 'w') as f:
            out = out.decode("utf-8")
            f.write(out)
        if err:
            with open(ec50_errfile, 'w') as f:
                err = err.decode("utf-8")
                f.write(err)
    

if __name__ == "__main__":
    main()
    