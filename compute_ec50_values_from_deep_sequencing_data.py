"""
Use deep-sequencing data from the protein-stability assay and compute EC50 values.

This script carries out the entire pipeline for computing EC50 values from deep-
sequencing data using a variety of `Python` scripts.

* Inputs:
    * `designed_seqs_file` : a file giving the DNA and protein sequence for each
    input design for a given round of designs

    * `experimental_summary_file` : a CSV summary file with the metadata, where
    samples are rows and columns are different pieces of metadata. The columns
    must include:
        * `round` (rename `replicate`): experimental replicate
        * `protease` : the protease associated with the sample
        * `selection_index` (rename `selection_strength`): the index of the selection step associated with the sample
        * `experiment/round` (rename `FASTQ_ID`) : a string that is both common to all FASTQ files for a given sample and is *not* found in any other sample, such that this string can be used to search for and retrieve all relevant FASTQ files
            * alternatively, it would be nice to have a list of FASTQ files associated with each experiment

        * ideally, the above CSV file would also include the following info from the `experiments.csv` file associated with Gabe's paper:
            * `parent` : the sample that is the parent to the previous sample
            * `conc_factor` : the fold-change in protease concentration between each step of the experiment
            * `parent_expression`:
            * `fraction_collected`:
            * `cells_collected`: total number of cells collected

    * `fastq_dir` : the relative path between the current working directory and the directory that has the FASTQ files



* Dependencies:
    * all dependencies in the below programs.

* Outputs:
    * Assembled FASTQ files for each sample
    * Counts files for each sample, one giving all counts of all proteins, the other only giving counts for the input designs
    * EC50 values for each sample

"""

# Import `Python` modules
import os
import sys
import pandas

# Read in command-line arguments
designed_seqs_file = sys.argv[1]
experimental_summary_file = sys.argv[2]
fastq_dir = sys.argv[3]
pear_path = sys.argv[4]

# Read in data from the experiments summary file
proteases = ['trypsin', 'chymotrypsin']
summary_df = pandas.read_csv(experimental_summary_file)
# selection_indices =  # get these from the experimental summary file

# Read in the designed sequences
designed_seqs_df = pandas.read_csv(designed_seqs_file, sep='\s+')
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
