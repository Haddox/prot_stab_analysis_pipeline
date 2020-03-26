"""
Compute EC50 values from deep-sequencing and FACS data from the high-throughput protein-stability assay
"""

# Import `Python` modules
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs).03d %(name)s %(message)s",
    datefmt='%Y-%m-%dT%H:%M:%S'
)

import os
from os import path
import sys
import glob
import argparse

import collections
import pickle
import pandas
import numpy as np
import scipy.stats
import itertools
from itertools import product

# Import custom functions
import protease_sequencing_model

# Define auxiliary functions
def none_or_int(x):
    if x == None:
        return None
    return int(x)

def none_or_div(x1, x2):
    if (x1 == None) or (x2 == None):
        return None
    return x1 / x2

def fix_num_selected(x):
    if x == None:
        return 5e5
    return x

def dict_product(iterable_dict):
    return (dict(zip(iterable_dict, x)) for x in itertools.product(*iter(iterable_dict.values())))

def fit_model(model_input, dataset, parameters):
    model = (
        protease_sequencing_model.FractionalSelectionModel(**dict(parameters))
        .build_model(model_input[dataset])
    )

    return model.find_MAP()

def report_model_ec50(model_input, counts, dataset, model_parameters, fit_parameters, output_dir):
    model = (
        protease_sequencing_model.FractionalSelectionModel(**dict(model_parameters))
        .build_model(model_input[dataset])
    )

    counts_df = counts[dataset]
    counts_df['ec50'] = fit_parameters['sel_ec50']

    cis = np.array([
        model.estimate_ec50_cred(fit_parameters, i, cred_spans=[.95])["cred_intervals"][.95]
        for i in range(len(counts_df))
    ])

    predictions = {
        p : {
            k : ps[k](fit_parameters)
            for k in ("Frac_sel_pop", "P_cleave")
        }
        for p, ps in list(model.model_populations.items())
    }

    sum_llh=np.zeros(len(counts_df))
    sum_signed_llh=np.zeros(len(counts_df))

    for p, ps in list(predictions.items()):
        counts_df['downsamp_counts%s' % p] = model.population_data[p]["P_sel"]
        counts_df['pred_counts%s' % p] = np.round(model.population_data[p]["P_sel"].sum() * predictions[p]['P_cleave'])
        bn=scipy.stats.binom(n=model.population_data[p]["P_sel"].sum(),p=predictions[p]['P_cleave'])
        my_llh = bn.logpmf(counts_df['downsamp_counts%s' % p])
        best_llh = bn.logpmf(counts_df['pred_counts%s' % p])
        counts_df['delta_llh%s' % p] = my_llh - best_llh
        counts_df['signed_delta_llh%s' % p] = counts_df['delta_llh%s' % p] * np.sign( counts_df['downsamp_counts%s' % p] - counts_df['pred_counts%s' % p]  )
        sum_llh +=  counts_df['delta_llh%s' % p]
        sum_signed_llh += counts_df['signed_delta_llh%s' % p]
    counts_df['sum_delta_llh'] = sum_llh
    counts_df['sum_signed_delta_llh'] = sum_signed_llh

    counts_df["ec50_95ci_lbound"] = cis[:,0]
    counts_df["ec50_95ci_ubound"] = cis[:,1]
    counts_df["ec50_95ci"] = cis[:,1] - cis[:,0]
    counts_df['sel_k'] = dict(model_parameters)['sel_k'] # dict(model_parameters)['sel_k']


    counts_df=counts_df[['name'] +
                        ['counts%s' % p for p in sorted(model.population_data)] +
                        ['downsamp_counts%s' % p for p in sorted(model.model_populations)] +
                        ['pred_counts%s' % p for p in sorted(model.model_populations)] +
                        ['delta_llh%s' % p for p in sorted(model.model_populations)] +
                        ['signed_delta_llh%s' % p for p in sorted(model.model_populations)] +
                        ['sel_k','sum_delta_llh','sum_signed_delta_llh','ec50_95ci_lbound','ec50_95ci_ubound',
                         'ec50_95ci','ec50']]

    output_file = '{0}/{1}.fulloutput'.format(
        output_dir,
        dataset
    )
    print("Writing the output data to a file called: {0}".format(output_file))
    counts_df.to_csv(output_file, index=False, sep='\t')

    return counts_df

# The main code to run the analysis
def main():
    """Read in command-line arguments and execute the main code"""

    #---------------------------------------------------------------
    # Read in command-line arguments and experimental metadata,
    # including deep-sequencing counts and FACS data
    #---------------------------------------------------------------
    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts_dir", help="a path the directory with files giving counts for each variant in the library")
    parser.add_argument("--experimental_summary_file", help="a path to the input experimental summary file")
    parser.add_argument("--datasets", help="a string with comma-separated datasets (e.g., 'dataset1,dataset2,etc.'")
    parser.add_argument("--output_dir", help="a path to an output directory where all the results will be stored. This directory will be made if it does not already exist.")
    args = parser.parse_args()

    # Assign command-line arguments to variables
    counts_dir = args.counts_dir
    experimental_summary_file = args.experimental_summary_file
    datasets = args.datasets.split(',')
    output_dir = args.output_dir
    print("\nWill analyze the datasets: {0}".format(', '.join(datasets)))

    # Initialize a the output directory if it doesn't already exist
    print("\nWill store the resulting EC50 values in the directory: {0}".format(output_dir))
    if not os.path.isdir(output_dir):
        print("Making the output direcotry".format(output_dir))
        os.makedirs(output_dir)

    # Read in metadata contained in the `experimental_summary_file`
    summary = pandas.read_csv(experimental_summary_file)
    summary = summary.dropna(how='all')
    summary = summary.fillna('-')
    summary = summary.where( summary != '-', None)
    
    # Read in deep-sequencing counts from each experiment
    print("\nReading in counts from files stored in the directory: {0}".format(counts_dir))
    model_input = {}
    for i, r in summary.iterrows():
        # Name of experiment (e.g., "rd1_tryp.counts")
        name = r["input"].replace(".counts", "")
        if name not in model_input:
            model_input[name] = {}
        # Get the selection round (e.g., counts0-counts6)
        rnd = int(r["column"].replace("counts", ""))
        assert rnd not in model_input[name]

        # Get metadata for a given experiment and selection round including which
        # library served as input (`parent`), what the protease strength was used
        # (`selection_level`), how many cells collected (`cells_collected`; default
        # is 5e5 if no value is provided), the fraction of collected cells before
        # and after protease treatment (`parent_expression`, `fraction_collected`;
        # the default is `None` if no data is provided), and the factor by which the
        # protease concentration was increased between experiments (e.g., 3-fold
        # for most experiments). Store in dictionary keyed by experiment name (e.g.,
        # `rd1_tryp`) and then by selection round (0-6).
        model_input[name][rnd] = dict(
            parent            = r['parent'],
            selection_level   = r['selection_strength'],
            num_selected      = fix_num_selected(r['cells_collected']),
            Frac_sel_pop = none_or_div(r['fraction_collected'], r['parent_expression']),
            conc_factor       = r['conc_factor']
        )

        # If there is data on the fraction of the population that was selected
        # (`Frac_sel_pop`) and if this fraction is above 0.999, then replace
        # it with 0.999 since, in theory, this fraction should not be above 1.0
        if model_input[name][rnd]['Frac_sel_pop'] != None:
            model_input[name][rnd]['Frac_sel_pop'] = min(0.999, model_input[name][rnd]['Frac_sel_pop'])
            
            # Next, correct `Frac_sel_pop` with a factor computed from the
            # fraction of matching sequences post- vs. pre- selection. Only
            # do this if the correction factor isn't NaN.
            corr_factor = 'matching_sequences_corr_factor'
            if not np.isnan(r[corr_factor]):
                print("Correcting Frac_sel_pop with {0} for {1} {2}".format(
                      corr_factor, r['input'], r['selection_strength']
                ))
                model_input[name][rnd]['Frac_sel_pop'] *= r[corr_factor]

        # If there is data for `matching_sequences`, then estimate the number of
        # total cells collected that actually have designs matching one of the input
        # designs
        if r['matching_sequences'] != None:
            model_input[name][rnd]['num_selected']  *= r['matching_sequences']

    # Read in the counts files for each round of selection with each protease
    # (e.g. `rd1_tryp.counts`), which has counts for each variant (rows) at each
    # selection level (columns; 0-6). Store in dictionary keyed by experiment
    # (e.g, `rd1_tryp`)
    counts = {
        exper : pandas.read_csv(
            path.join(counts_dir, exper + ".counts"),
            delim_whitespace=True)
        for exper in model_input
    }

    # Iterate over experiments and read in counts data
    for exper in model_input:
        # Read in the counts for a given experiment, as above
        counts_df = pandas.read_csv(
                path.join(counts_dir, exper + ".counts"),
                delim_whitespace=True)

        # Iterate through columns, excluding the first column (`name`), leaving
        # seven columns corresponding to the seven levels of selection (0-6)
        for i, col in enumerate(counts_df.columns[1:]):
            # For a given experiment and selection strength, add an ordered list
            # of all variant names (e.g., `EHEE_0401.pdb`)
            model_input[exper][i]["name"] = counts_df["name"]
            # ... as well as an ordered list of all sequencing counts matching those
            # sequences
            model_input[exper][i]['seq_counts'] = counts_df[col].astype(np.int)

        # Iterate through each selection level compute `P_sel`
        for k, v in list(model_input[exper].items()):
            # If the total number of sequencing counts is greater than the total
            # number of cells collected (specifically, the number of cells expected
            # to have designs that actually match one of the input designs), then
            # compute frequencies from counts (`pfrac`) and then counts that are
            # normalized to the number of cells selected, rounding down to the
            # nearest integer.
            if v["seq_counts"].sum() > v["num_selected"]:
                pfrac = v["seq_counts"].values.astype(float) / v['seq_counts'].sum()
                v["P_sel"] = np.floor(pfrac * v["num_selected"])
            # Otherwise, if the sequencing depth is less than the number of cells
            # selected, set the number of selected cells equal to the sum of the
            # sequencing counts. Then compute `P_sel` as just the raw sequencing
            # counts.
            else:
                v["num_selected"] = v["seq_counts"].sum()
                v["P_sel"] = np.array(v["seq_counts"].astype(np.float))
            # Apparently, `min_fraction` is always set to None, at least to begin
            # with
            v['min_fraction'] = None


    #---------------------------------------------------------------
    # Compute EC50 values from the input data
    #---------------------------------------------------------------
    # Define features of the parameter space used in the analysis
    param_space = dict(
        response_fn = ("NormalSpaceErfResponse",),
        min_selection_mass = [5e-7],
        min_selection_rate = [0.0001],
        outlier_detection_opt_cycles = [3],
        sel_k = [0.8]
    )
    parameter_sets = [frozenset(list(d.items())) for d in dict_product(param_space)]

    # Iterate through each experiment, computing EC50 values for each sequence
    for d, p in product(datasets, parameter_sets):
        print ("Computing EC50s for the dataset: {0}".format(d))
        my_fit_model = (
            model_input,
            counts,
            d,
            p,
            fit_model(model_input, d, p),
            output_dir
        )
        report_model_ec50(*my_fit_model)

if __name__ == "__main__":
    main()
