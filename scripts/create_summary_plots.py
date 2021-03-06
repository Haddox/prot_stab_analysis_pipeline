"""
Make plots that summarize the results of the high-throughput protein-stability assay

These plots are generated using the output from the script called: `compute_ec50_values_from_deep_sequencing_data.py`
"""

# Import `Python` modules
import os
import argparse
import glob
import re
import pandas
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("ticks")

def parse_PARE_outfile(outfile):
    """
    Parse output data generated by PARE when assemblying paired-end reads
    
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

def summarize_PARE_assembly(pare_log_dir, experiment_names, out_file):
    
    # Find all .log files from PARE in the specified directory
    log_files = glob.glob(os.path.join(pare_log_dir, '*.log'))
    
    # For each log file, look for an experimental name that matches the beginning
    # of the base name of the log file, making sure that each log file only matches
    # one experiment name. Organize log files by experiment name in `log_files_dict`
    log_files_dict = {
        experiment_name : []
        for experiment_name in experiment_names
    }
    for log_file_name in log_files:
        matching_experiments = []
        for experiment_name in log_files_dict:
            log_file_basename = os.path.basename(log_file_name)
            if log_file_basename.find('{0}-'.format(experiment_name)) == 0:
                log_files_dict[experiment_name].append(log_file_name)
                matching_experiments.append(experiment_name)
        if len(matching_experiments) == 0:
            raise ValueError("Could not find a matching experiment for the log file: {0}".format(log_file_name))
        if len(matching_experiments) > 1:
            print(matching_experiments)
            raise ValueError("Found multiple matching experiments for the log file: {0}".format(log_file_name))
        assert len(matching_experiments) == 1

    # For each experiment, read in depth and quality scores from `PEAR` .log files associated with that experiment and store the results in a pandas dataframe
    assembly_d = {
        key : []
        for key in ['experiment_name', 'n_assembled_reads', 'n_discarded_reads', 'n_non_assembled_reads']
    }
    
    for experiment_name in log_files_dict:
        assembly_d['experiment_name'].append(experiment_name)
        n_assembled_reads = n_discarded_reads = n_non_assembled_reads = 0
        for log_file_name in log_files_dict[experiment_name]:
            (n_assembled_reads_i, n_discarded_reads_i, n_non_assembled_reads_i) = \
                parse_PARE_outfile(log_file_name)
            n_assembled_reads += n_assembled_reads_i
            n_discarded_reads += n_discarded_reads_i
            n_non_assembled_reads += n_non_assembled_reads_i
        assembly_d['n_assembled_reads'].append(n_assembled_reads)
        assembly_d['n_discarded_reads'].append(n_discarded_reads)
        assembly_d['n_non_assembled_reads'].append(n_non_assembled_reads)

    assembly_df = pandas.DataFrame.from_dict(assembly_d)

    # Plot the data for each replicate as stacked bar charts

    # Get data for each bar
    assembly_df.set_index('experiment_name', inplace=True)
    labels = sorted(assembly_df.index.values)
    first_bar = assembly_df.loc[labels]['n_assembled_reads']
    second_bar = assembly_df.loc[labels]['n_non_assembled_reads']
    third_bar = assembly_df.loc[labels]['n_discarded_reads']
    assert(len(first_bar) == len(second_bar) == len(third_bar))

    # Make plot
    plot_indices = np.arange(len(first_bar))
    barwidth = len(first_bar) * [0.75]
    plt.barh(plot_indices, first_bar, barwidth, label='assembled', align='center')
    plt.barh(plot_indices, second_bar, barwidth, left=first_bar, label='non-assembled', color='red', align='center')
    plt.barh(plot_indices, third_bar, barwidth, left=second_bar, label='discarded', color='purple', align='center')
    plt.yticks(plot_indices, labels)
    plt.xlabel('Number of reads')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    def adjust_ylim(l):
        mn, mx = l
        mn += 1.25
        mx -= 0.5
        return (mn, mx)
    plt.ylim(adjust_ylim(plt.ylim()))
    plt.yticks()
    sns.despine()
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig(out_file, bbox_extra_artists=(lgd,), bbox_inches='tight')


# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""
    
    #---------------------------------------------------------------
    # Read in command-line arguments and experimental metadata
    #---------------------------------------------------------------
    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", help="the directory the contains the output from the script `compute_ec50_values_from_deep_sequencing_data.py`")
    parser.add_argument("--output_dir", help="a path to an output directory where all the summary plots will be stored. This directory will be created if it does not already exist")
    args = parser.parse_args()    
    
    # Assign command-line arguments to variables
    data_dir = args.data_dir
    output_dir = args.output_dir
    print("\nHere is a list of parsed input arguments:")
    for arg in vars(args):
        print("{0}: {1}".format(arg, getattr(args, arg)))
    if not os.path.isdir(output_dir):
        print("Making the output directory: {0}".format(output_dir))
        os.makedirs(output_dir)
    
    # Read in the `experiments.csv` file in the `ec50_values` and make
    # a column called `experiment_name` that is a string with the following
    # information: {protease}_{selection_strength}
    experimental_summary_file = os.path.join(data_dir, 'ec50_values/experiments.csv')
    experimental_summary_df = pandas.read_csv(experimental_summary_file)
    experimental_summary_df['experiment_name'] = experimental_summary_df.apply(
        lambda row: row['input'].replace('.counts', '') + '_' + str(row['selection_strength']),
        axis=1
    )
    
    #---------------------------------------------------------------
    # Make a plot summarizing the results of using `PARE` to assemble
    # the deep-sequencing reads, reporting read depth and quality for
    # each sample
    #---------------------------------------------------------------        
    experiment_names = sorted(list(experimental_summary_df['experiment_name']))
    assert len(experiment_names) == len(set(experiment_names)), "Found duplicate experiment names"
    pare_log_dir = os.path.join(data_dir, 'paired_FASTQ_files/')
    out_file = os.path.join(output_dir, 'deep_sequencing_depth_and_quality.png')
    print("\nMaking the file: {0}".format(out_file))
    summarize_PARE_assembly(pare_log_dir, experiment_names, out_file)

if __name__ == "__main__":
    main()