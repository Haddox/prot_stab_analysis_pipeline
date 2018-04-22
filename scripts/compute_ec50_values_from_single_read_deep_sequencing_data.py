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
    ec50s_dir = os.path.join(output_dir, 'ec50_values')
    stability_scores_dir = os.path.join(output_dir, 'stability_scores')
    dirs = [
        output_dir, paired_FASTQ_files_dir, counts_dir, ec50s_dir,
        stability_scores_dir
    ]
    for dir_i in dirs:
        if not os.path.isdir(dir_i):
            print("Making the directory: {0}".format(dir_i))
            os.makedirs(dir_i)

    # Read the data from the experiments summary file into a dataframe
    # and get a list of unique proteases and unqiue selection indices
    summary_df = pandas.read_csv(experimental_summary_file)
    proteases = set(summary_df['protease_type'])
    selection_indices = set(summary_df['selection_strength'])


    #---------------------------------------------------------------
    # Find all FASTQ files reads for each sample
    #---------------------------------------------------------------

    # Iterate through each sample in the experimental summary dataframe,
    # make a list of FASTQ files for each sample
    FASTQ_files = {}
    input_data_for_computing_counts = []
    for (i, row) in summary_df.iterrows():

        # Get sample metadata
        protease_type = row['protease_type']
        selection_index = row['selection_strength']
        experiment_name = '{0}_{1}'.format(protease_type, selection_index)
        fastq_id = row['fastq_id']

        # Find all R1 files
        print("Searching for FASTQ files for the experiment {0} using the FASTQ ID {1}".format(
            experiment_name, fastq_id
        ))
        r1_files = glob.glob('{0}/{1}*_R1_*.fastq*'.format(fastq_dir, fastq_id))
        print("Here is a list of R1 files: {0}".format(', '.join(r1_files)))

        # Assemble data for computing counts
        output_counts_file = os.path.join(counts_dir, '{0}_counts.csv'.format(experiment_name))
        if os.path.isfile(output_counts_file):
            print("Already have counts for the sample: {0}".format(experiment_name))
            continue
        else:
            print("Computing counts for the sample: {0}".format(experiment_name))
            input_data_for_computing_counts.append((r1_files, output_counts_file))


    #---------------------------------------------------------------
    # Compute protein counts from the deep-sequencing data for each sample
    #---------------------------------------------------------------
    # For each sample, compute protein-level counts from the deep-sequencing data. Write
    # a file giving all observed counts, even for proteins that don't match a starting
    # design. Compute the counts in parallel, reporting the progress of the computation
    print("\nComputing protein counts from the deep-sequencing data")
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
    # These input files include a single `experiments.csv` file with experimental metadata
    # including the number of protein counts that match a starting sequence. These files
    # also include one file per protease that reports counts at each selection step only
    # for proteins that match one of the input designs. To make these input files, I will
    # first read the designed sequences into a dataframe.
    designed_seqs_df = pandas.read_csv(designed_sequences_file)
    designed_seqs_df.set_index('protein_sequence', inplace=True)

    # Next, I will initiate a dictionary for keeping track of information for writing the
    # `experiments.csv` file
    counts_metadata_dict = {
        key : []
        for key in ['protease_type', 'selection_strength', 'input', 'column', 'matching_sequences']
    }

    # Next, for each protease, I will aggregate the counts and record metadata
    for protease in proteases:

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
        # starting counts in the naive library (`counts0`) that are greater than a
        # cutoff (`counts_cutoff`) that is currently set to zero
        counts_cutoff = 0
        aggregate_df.set_index('name', inplace=True)
        columns_to_write = [
            'counts{0}'.format(selection_index)
            for selection_index in selection_indices
        ]
        aggregate_df[
            aggregate_df['counts0']>counts_cutoff
        ][columns_to_write].to_csv(aggregate_counts_outfile, sep=' ')

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


    #---------------------------------------------------------------
    # Compute stability scores from the EC50 values using the script
    # called `compute_stability_scores_from_EC50_values.py`
    #---------------------------------------------------------------
    # Do this for each protease independently of the other
    for protease in proteases:

        # Get the concentration factor from the input experimental
        # summary file, and make sure that all factors for a given
        # protease are the same. In theory they could be different,
        # but the analysis framework doesn't currently support that
        conc_factors = set(summary_df[
            (summary_df.reset_index()['protease_type'] == protease) &\
            (summary_df['conc_factor'].notnull())
        ]['conc_factor'])
        assert len(conc_factors) == 1, "Not all concentration factors are the same for the protease {0}".format(protease)
        conc_factor = str(int(list(conc_factors)[0]))

        # Define the input file with EC50 values, the output file to
        # be created, and assemble the entire command-line argument
        output_file = os.path.join(stability_scores_dir, '{0}_stability_scores.txt'.format(protease))
        ec50_values_file = os.path.join(ec50s_dir, '{0}.fulloutput'.format(protease))
        cmd = [
            'python',
            '{0}/compute_stability_scores_from_EC50_values.py'.format(scriptsdir),
            designed_sequences_file,
            protease,
            ec50_values_file,
            conc_factor,
            output_file
        ]

        # Carry out the command
        print("\nComputing stability scores for the protease {0} with the command: {1}".format(protease, ' '.join(cmd)))
        stability_scores_logfile = os.path.join(stability_scores_dir, '{0}_stability_scores.log'.format(protease))
        stability_scores_errfile = os.path.join(stability_scores_dir, '{0}_stability_scores.err'.format(protease))
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out, err = process.communicate()
        with open(stability_scores_logfile, 'w') as f:
            out = out.decode("utf-8")
            f.write(out)
        if err:
            with open(stability_scores_errfile, 'w') as f:
                err = err.decode("utf-8")
                f.write(err)


if __name__ == "__main__":
    main()
