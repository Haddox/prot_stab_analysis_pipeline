"""
This script complies input data for computing EC50 values

This script is heavily adapted from (Rocklin et al., 2017, Science). It compiles
multiple pieces of data used to compute EC50 values, including deep-sequencing
data and FACS data, described below in more detial. In the Rocklin et al. study,
this script was included in the same directory as these input data. In my
analysis, I have organized things a little differently. Instead, I gather data
from other directories using hard-coded paths defined below.

Using these hard-coded paths, this script reads in the starting data and stores
it in a dictionary called `model_input`, which is returned when this script
is called as a module.

Inputs (hard-coded paths):
	* a path to a directory with files of protein counts derived from deep-
    sequencing data from each step of each replicate experiment and a file
    called `experiments.csv`, which lists all steps of all replicate experiments
    and associated FACS data and data on the fraction of deep-sequencing counts
    that match one of the input protein designs or controls.

Outputs:
    * a dictionary called `model_input` that contains the input deep-sequencing
    and FACS data for each step of each replicate experiment. Upon importing
    this directory as a module, this dictionary will be imported under the name
    of the module.
"""


import pandas as pd
import numpy as np
import pickle
import sys
import glob
import os
from os import path

import collections

# When called as a module, this script will return the dictionary called
# `model_input`
__all__ = ["model_input"]

# Hard-coded input paths to a directory with counts files derived from the
# deep-sequencing data and a file called `experiments.csv`
basedir = path.dirname(__file__)
experiments_csv_dir = os.path.join(
    basedir, '../../results/ec50_values/'
)
countsdir = os.path.join(
    basedir, '../../results/counts/'
)

# Start compiling the data
def none_or_int(x):
    if x == None:
        return None
    return int(x)

def none_or_div(x1, x2):
    if x1 == None:
        return None
    return x1 / x2

def fix_num_selected(x):
    if x == None:
        return 5e5
    return x

# Read in the `experiments.csv` file
summary=pd.read_csv(path.join(experiments_csv_dir, 'experiments.csv'))
summary=summary.dropna(how='all')
summary=summary.fillna('-')
summary=summary.where( summary != '-', None)

# Read in metadata for each experiment.
model_input = {}

for i, r in summary.iterrows():
    # Name of experiment (e.g., "rd1_tryp.counts")
    name = r["input"].replace(".counts", "")
    # Expect multiple lines for the same `name`
    if name not in model_input:
        model_input[name] = {}
    # Selection round (e.g., counts0-counts6)
    rnd = int(r["column"].replace("counts", ""))
    # Make sure there aren't any duplicate entries for a given experiment and
    # selection round
    assert rnd not in model_input[name]

    # Get metadata for a given experiment and selection round including which
    # library served as input (`parent`), what the protease strength was used
    # (`selection_level`), how many cells collected (`cells_collected`; default
    # is 5e5 if no value is provided), the fraction of collected cells before
    # and after protease treatment (`parent_expression`, `fraction_collected`;
    # the default is `None` if no data is provided), and the factor by which the
    # protease concentration was increased between experiments (e.g., 3-fold
    # for most experiments). Sotre in dictionary keyed by experiment name (e.g.,
    # `rd1_tryp`) and then by selection round (0-6).
    model_input[name][rnd] = dict(
        parent            = r['parent'],
        selection_level   = r['selection_strength'],
        num_selected      = fix_num_selected(r['cells_collected']),
        Frac_sel_pop = none_or_div(r['fraction_collected'], r['parent_expression']),
        conc_factor       = r['conc_factor']
    )

    # If there is data on the fraction of the population that was selected,
    # (`Frac_sel_pop`), then estimate the fraction that actually had designs
    # that actually match one of the input designs (i.e., don't have sequence
    # errors) by multiplying `Frac_sel_pop` by `matching_sequences`.
    if model_input[name][rnd]['Frac_sel_pop'] != None:
        model_input[name][rnd]['Frac_sel_pop'] *= r['matching_sequences']

    # If there is data for `matching_sequences`, then estimate the number of
    # total cells collected that actually have designs matching one of the input
    # designs
    if r['matching_sequences'] != None:
        model_input[name][rnd]['num_selected']  *= r['matching_sequences']

# Read in the counts files for each round of selection with each protease
# (e.g. `rd1_tryp.counts`), which has counts for each protein (rows) at each
# selection level (columns; 0-6). Store in dictionary keyed by experiment
# (e.g, `rd1_tryp`)
counts = {
    exper : pd.read_csv(
        path.join(countsdir, exper + ".counts"),
        delim_whitespace=True)
    for exper in model_input
}

# Iterate over experiments and read in counts data
for exper in model_input:
    # Read in the counts for a given experiment, as above
    counts_df = pd.read_csv(
            path.join(countsdir, exper + ".counts"),
            delim_whitespace=True)

    # Iterate through columns, excluding the first column (`name`), leaving
    # seven columns corresponding to the seven levels of selection (0-6)
    for i, col in enumerate(counts_df.columns[1:]):
        # For a given experiment and selection strength, add an ordered list
        # of all protein names (e.g., `EHEE_0401.pdb`)
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
