"""
Compute stability scores from EC50 values
"""

# Import `Python` modules
import os
import sys
import argparse
import pandas
import numpy as np
import glob
import sequence_protease_susceptibility
from sequence_protease_susceptibility import CenterLimitedPSSMModel

# Define auxilliary functions
def model_from_df(df):
    """
    Get parameters for the unfolded-state model from an input dataframe
    """
    # Make a list of characters and numbers with most amino acids and some non-
    # amino-acid characters
    aa_nums = list(range(21))
    aa_nums.remove(1)
    aas = [x for x in 'ADEFGHIKLMNPQRSTVWYZ']

    # Get model parameters
    blank_model = CenterLimitedPSSMModel(flanking_window=4)
    blank_model.setup()
    blank_model.fit_coeffs_ = np.zeros((), blank_model.coeff_dtype)
    blank_model.fit_coeffs_['c0'] = df.query('variable == "c0"')['value'].astype(float).values[0]
    blank_model.fit_coeffs_['EC50_max'] = df.query('variable == "EC50_max"')['value'].astype(float).values[0]
    blank_model.fit_coeffs_['k_max'] = df.query('variable == "k_max"')['value'].astype(float).values[0]
    blank_model.fit_coeffs_['P1_PSSM'] = np.zeros((21))
    blank_model.fit_coeffs_['outer_PSSM'] = np.zeros((9,21))
    for rownum, row in enumerate(['P5','P4','P3','P2','P1',"P1'","P2'","P3'","P4'"]):
        for aa_num, aa in zip(aa_nums, aas):
            if row == 'P1':
                blank_model.fit_coeffs_['P1_PSSM'][aa_num] = df.query('variable == "P1" & aa == "%s"' % aa)['value'].astype(float).values[0]
                #print len(df.query('variable == "P1" & aa == "%s"' % aa)['value'].astype(float).values)
            else:
                blank_model.fit_coeffs_['outer_PSSM'][rownum][aa_num] = df.query('variable == "%s" & aa == "%s"' % (row, aa))['value'].astype(float).values[0]
                #print len(df.query('variable == "%s" & aa == "%s"' % (row, aa))['value'].astype(float).values)
    return blank_model

# Run the main code
def main():
    """Read in command-line arguments and run the main code of the script"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser(description="Compute stability scores from EC50 values")
    parser.add_argument("designed_sequences_file", help="a file with design names and protein sequences")
    parser.add_argument("protease", help="the relevant protease that will be used to compute predicted unfolded-state EC50 values, must be either 'trypsin' or 'chymotrypsin'")
    parser.add_argument("ec50_values_file", help="a file with ec50 values in a matrix with tab-delimited columns")
    parser.add_argument("output_file", help="a path to an output file where the results will be stored. The directories in this path must already exist.")
    args = parser.parse_args()

    # Assign command-line arguments to variables
    designed_sequences_file = args.designed_sequences_file
    protease = args.designed_sequences_file
    if protease not in ['trypsin', 'chymotrypsin']:
        raise ValueError("The `protease` argument must be either 'trypsin' or 'chymotrypsin', instead it is: {0}".format(protease))
    output_file = args.output_file

    # Set up a protease-specific unfolded-state models for the input protease,
    # reading in the specific model parameters fit in the Rocklin et al. study,
    # which are encoded in the file called `unfolded_state_model_params`
    print("Setting up the unfolded-state model for the protease: {0}".format(protease))
    params_input_file = 'unfolded_state_model_params'
    model_data = pandas.read_csv(params_input_file, delim_whitespace=True)
    unfolded_state_model = model_from_df(
        model_data.query('protease == "{0}"'.format(protease))
    )

    # Read in input amino-acid sequences and input EC50 values and find the
    # intersection using the `name` column of each dataframe
    designed_seqs_df = pandas.read_csv(designed_sequences_file)
    designed_seqs_df.set_index('name')
    ec50s_df = pandas.read_csv(ec50_values_file, sep='\t')
    ec50s_df.set_index('name')
    stability_scores_df = ec50s_df.merge(designed_seqs_df, left_index=True, right_index=True, how='inner')
    print("Found {0} sequences in the input file with design sequences called: {0}".format(
        len(designed_seqs_df.index.values),
        designed_sequences_file
        )
    )
    print("Found {0} sequences in the input file with EC50 values called: {0}".format(
        len(ec50s_df.index.values),
        ec50_values_file
        )
    )
    print("Found {0} sequences with names that are present in both input files, based on the 'name' column of each file".format(
        len(stability_scores_df.index.values)
        )
    )

    # Compute the predicted unfolded-state EC50 for each of the sequences
    print("Computing predicted unfolded-state EC50 values")
    designed_seqs_df['ec50_pred'] = unfolded_state_model.predict([
        stability_scores_df['protein_sequence']
    ])

    # Verify that the model reproduces a past unfolded-EC50 predictions for a
    # previously tested sequence HEEH_rd4_0097 in the rd 4 library
    #seq = "DVEEQIRRLEEVLKKNQPVTWNGTTYTDPNEIKKVIEELRKSMLESSGGS"
    #fullseq = "GGGSASHM" + seq + "LEGGGSEQ"
    #print ('trypsin', tryp.predict([fullseq]))
    #print ('chymotrypsin', chymo.predict([fullseq]))
    #should return 1.62 for tryp and 1.85 for chymo

if __name__ == "__main__":
    main()
