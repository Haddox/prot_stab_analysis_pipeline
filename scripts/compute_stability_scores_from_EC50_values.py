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
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Get the path to this script
scriptsdir = os.path.dirname(__file__)

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
    parser.add_argument("designed_sequences_file", help="a file with design names and sequences")
    parser.add_argument("protease", help="the relevant protease that will be used to compute predicted unfolded-state EC50 values, must be either 'trypsin' or 'chymotrypsin'")
    parser.add_argument("ec50_values_file", help="a file with ec50 values in a matrix with tab-delimited columns")
    parser.add_argument("conc_factor", help="the fold-change in protease concentration between selection steps (integer)", type=int)
    parser.add_argument("output_file", help="a path to an output file where the results will be stored. The directories in this path must already exist.")
    parser.add_argument("--protein_or_dna_level", default='protein', help="whether `designed_sequences_file` contains DNA or protein sequences. Options are: 'protein' or 'dna'. Default is: 'protein'")

    args = parser.parse_args()

    # Assign command-line arguments to variables
    designed_sequences_file = args.designed_sequences_file
    protease = args.protease
    if protease not in ['trypsin', 'chymotrypsin']:
        raise ValueError("The `protease` argument must be either 'trypsin' or 'chymotrypsin', instead it is: {0}".format(protease))
    ec50_values_file = args.ec50_values_file
    conc_factor = args.conc_factor
    output_file = args.output_file
    protein_or_dna_level = args.protein_or_dna_level
    print("\nHere is a list of parsed input arguments")
    for arg in vars(args):
        print("{0}: {1}".format(arg, getattr(args, arg)))

    # Set up a protease-specific unfolded-state models for the input protease,
    # reading in the specific model parameters fit in the Rocklin et al. study,
    # which are encoded in the file called `unfolded_state_model_params`
    params_input_file = os.path.join(scriptsdir, 'unfolded_state_model_params')
    print("\nSetting up the unfolded-state model for the protease {0} using the parameters encoded in the file: {1} ".format(protease, params_input_file))
    model_data = pandas.read_csv(params_input_file, delim_whitespace=True)
    unfolded_state_model = model_from_df(
        model_data.query('protease == "{0}"'.format(protease))
    )

    # Read in input sequences and input EC50 values and find the
    # intersection using the `name` column of each dataframe
    designed_seqs_df = pandas.read_csv(designed_sequences_file)
    designed_seqs_df.set_index('name', inplace=True)
    ec50s_df = pandas.read_csv(ec50_values_file, sep='\t')
    ec50s_df.set_index('name', inplace=True)
    stability_scores_df = ec50s_df.merge(
        designed_seqs_df,
        left_index=True,
        right_index=True,
        how='inner'
    )

    # Report the overlap in sequence names in the two input files
    print("Found {0} sequences in the input file with design sequences called: {1}".format(
        len(designed_seqs_df.index.values),
        designed_sequences_file
        )
    )
    print("Found {0} sequences in the input file with EC50 values called: {1}".format(
        len(ec50s_df.index.values),
        ec50_values_file
        )
    )
    print("Found {0} sequences with names that are present in both input files, based on the 'name' column of each file".format(
        len(stability_scores_df.index.values)
        )
    )

    # If the designed sequences are DNA sequences, add a column that gives the
    # translated sequence
    if protein_or_dna_level == 'protein':
        assert 'protein_sequence' in designed_seqs_df.columns.values, "Expected there to be a column called 'protein_sequence' in the file {0}".format(designed_sequences_file)
    elif protein_or_dna_level == 'dna':
        assert 'dna_sequence' in designed_seqs_df.columns.values, "Expected there to be a column called 'dna_sequence' in the file {0}".format(designed_sequences_file)
        stability_scores_df['protein_sequence'] = stability_scores_df['dna_sequence'].apply(
            lambda x: str(Seq(x, generic_dna).translate())
        )
    else:
        raise ValueError("Failed to parse the input variable `protein_or_dna_level`: {0}. Must be 'protein or 'dna'".format(protein_or_dna_level))

    # Add flanking sequences
    print("Adding on flanking sequences")
    stability_scores_df['full_protein_sequence'] = \
        stability_scores_df['protein_sequence'].apply(
            lambda x: "GGGSASHM" + x + "LEGGGSEQ"
        )

    # The unfolded-state model requires that all sequences be of the same
    # length. If sequences are not all the same length, the below code
    # determines the length of the longest sequence and then appends null
    # (`Z`) characters to the end of each of the sequences until they
    # match the max length.
    input_protein_sequences = list(stability_scores_df['full_protein_sequence'])
    protein_lengths = list(map(len, input_protein_sequences))
    if not set(protein_lengths) == 1:
        print("Evening out sequence lengths by adding null characters to the end of sequences shorter than the longest sequence")
        max_len = max(protein_lengths)
        elongated_input_protein_sequences = []
        for protein_sequence in input_protein_sequences:
            length_difference = max_len - len(protein_sequence)
            assert length_difference >= 0, "No sequence should be longer than the max length"
            elongated_sequence = protein_sequence + (length_difference * 'Z')
            elongated_input_protein_sequences.append(elongated_sequence)
        assert (len(elongated_input_protein_sequences) == len(input_protein_sequences)), "The list of elongated sequences is the wrong length"
        input_protein_sequences = elongated_input_protein_sequences

    # Compute the predicted unfolded-state EC50 for each of the sequences
    # by passing the unfolded-state model a list of protein sequences
    stability_scores_df['ec50_pred'] = list(unfolded_state_model.predict(
        input_protein_sequences
    ))

    # Next, for each sequence, compute the difference in the sequence's
    # observed EC50 value and the EC50 value predicted for the sequence in
    # an unfolded state calling this the `ec50_rise`, as in Rocklin et al.
    # From that, compute stability scores by transforming the difference
    # from a log3 scale to a log10 scale
    print("Computing stability scores based on differences in predicted and observed EC50s")
    print("Will use a concentration factor of {0} to compute stability scores".format(conc_factor))
    stability_scores_df['ec50_rise'] = stability_scores_df['ec50'] - \
        stability_scores_df['ec50_pred']
    stability_scores_df['stabilityscore'] = stability_scores_df['ec50_rise'].apply(
        lambda x: np.log10(np.power(conc_factor, x))
    )

    # Write the resulting dataframe with stability scores to an output file
    print("Writing the results to the file: {0}".format(output_file))
    stability_scores_df.to_csv(output_file, sep='\t')

if __name__ == "__main__":
    main()
