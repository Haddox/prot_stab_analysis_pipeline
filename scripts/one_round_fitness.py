#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd
import io
import argparse
import scipy.stats
import tabulate

parser = argparse.ArgumentParser()
parser.add_argument("naive_sort_name_file", type=str)
parser.add_argument("-collected_cells", type=int, default=100000)
parser.add_argument("-conf_interval_exp", type=float, default=-20)
parser.add_argument("-native_name", type=str, default=None)
parser.add_argument("-ignore_all_zero", action="store_true")
parser.add_argument("-force_select_ratio", type=float, default=1, 
                                        help="If we're not normalizing by a native. Set fitnesses"
                                             "such that counts.sum() / naive.sum() == force_select_ratio")
# parser.add_argument("-in:file:silents", type=str, nargs="*", default=None)
# parser.add_argument("-in:file:silent", type=str, default="")


args = parser.parse_args(sys.argv[1:])

with open(args.naive_sort_name_file) as f:
    input_file = f.read()

header = 0
try:
    float( input_file[:input_file.find(" ")] )
    header = None
except:
    pass

df = pd.read_csv(io.StringIO(input_file), sep="\s+", header=header)

assert(len(df.columns) == 3)

df.columns = ["naive", "counts", "description"]

if ( args.ignore_all_zero ):
    df = df[(df['naive'] + df['counts']) != 0].copy()



naive = df['naive'].values.copy()
counts = df['counts'].values.copy()
og_counts = counts

naive_frac = naive.astype(float) / naive.sum()

if ( counts.sum() / args.collected_cells > 1 ):
    counts = ( counts * args.collected_cells / counts.sum() ).astype(int)

native_idx = None
if ( not args.native_name is None ):
    is_right_name = df['description'] == args.native_name
    if ( is_right_name.sum() == 0 ):
        sys.exit("-native_name does not appear in file")
    if ( is_right_name.sum() > 1 ):
        sys.exit("-native_name appears multiple times in file")

    native_idx = np.where(is_right_name)[0][0]
    args.force_select_ratio = None


const_min_fit = -10
const_upper_fit = 10
exp_const_min_fit = np.exp(const_min_fit)
def fitness_function(fitness):
    return np.exp(fitness) - exp_const_min_fit

def inv_fitness_function(fitness_sel):
    return np.log(fitness_sel + np.exp(const_min_fit)).clip(const_min_fit+1, const_upper_fit-1)


x_fitnesses = np.linspace(const_min_fit+1, const_upper_fit-1, 200)

fitness_function_x_fitnesses = fitness_function(x_fitnesses)

# what the fitness_sel would have to be to get to the counts from the aa_prob
start_fitness_sel = counts.astype(float) / counts.sum() / np.clip(naive_frac, 0.00000000000001, None)

if ( not native_idx is None ):
    # normalize native to fitness 0
    start_fitness_sel = start_fitness_sel * fitness_function(0) / start_fitness_sel[native_idx]

if ( not args.force_select_ratio is None ):
    # normalize the selected fraction up to force_select_ratio

    start_fitness_sel = start_fitness_sel * args.force_select_ratio #/ (naive_frac * start_fitness_sel).sum()

    # assert( np.isclose( (naive_frac * start_fitness_sel).sum(), args.force_select_ratio ) )

# add a tiny amount so there's no 0s
p_sel_early = naive_frac * start_fitness_sel + 0.0000000000001
p_sel_early_sum = np.sum(p_sel_early)

# global variables are bad lol
def get_logp_traces(idx):
    
    # apply fitness function in reverse to get fitness_guesses
    # fitness_guesses = inv_fitness_function(start_fitness_sel)

    # this statement is true
    # ( naive_frac * fitness_function(fitness_guesses) ).sum() == args.force_select_ratio
    
    # calculate p_sel for idx with the different x_fitnesses

    use_p_sel_early_sum = p_sel_early_sum - p_sel_early[idx]

    my_p_sels = naive_frac[idx] * fitness_function_x_fitnesses + 0.0000000000001

    p_sel_sum = use_p_sel_early_sum + my_p_sels

    # p_sel = np.tile(p_sel_early, (len(x_fitnesses), 1))
    # p_sel[:,idx] = naive_frac[idx] * fitness_function(x_fitnesses) + 0.0000000000001
    # p_sel_sum = p_sel.sum(axis=-1)
    
    # this is the slow step, but we don't even use it
    # p_sel *= 1/p_sel_sum[...,None] 

    # this normalizes the whole pool to 1
    my_p_sels /= p_sel_sum

    log_p = scipy.stats.binom.logpmf(counts[idx], counts.sum(), my_p_sels)
    
    # log_p = scipy.stats.multinomial.logpmf(counts, counts.sum(), p_sel )
    
    return log_p


xs = x_fitnesses

fitnesses = np.zeros(len(df), np.float)
cis = np.zeros(len(df), np.float)


for i in range(len(df)):
    if ( i % 1000 == 0 ):
        print("%4i / %4i"%(i, len(df)))
    # if ( i > 5000 ):
    #     break

    raw_trace = get_logp_traces(i)
    logps = raw_trace - raw_trace.max()

    pmf = np.exp(logps) / np.sum(np.exp(logps), axis=-1)[...,None]
    cdf = np.cumsum(pmf, axis=-1)
    rcdf = np.cumsum(pmf[::-1])[::-1] 
    # print(cdf)
    super_span = args.conf_interval_exp

    super_low = xs[np.searchsorted(np.log(cdf.clip(1e-99, None)), super_span, side="left").clip(0, len(xs)-1)]
    super_high = xs[np.searchsorted(-np.log(rcdf.clip(1e-99, None)), -super_span, side="right").clip(0, len(xs)-1)]
    
    cis[i] = super_high - super_low
    fitnesses[i] = (super_high + super_low) / 2


df['adjusted_counts'] = counts
df['frac_sel'] = fitness_function( fitnesses )
df['confidence_interval'] = cis
df['fitness'] = fitnesses

columns = list(df.columns)
columns.remove("description")
columns.append("description")

df = df[columns].copy()




with open(args.naive_sort_name_file + ".fitness", "w") as f:
    numeric = df.select_dtypes(include=[np.number])
    non_numeric = df.select_dtypes(exclude=[np.number])
    values = np.concatenate((np.around(numeric, 3), non_numeric), axis=-1)

    f.write(tabulate.tabulate(values.tolist(), list(numeric) + list(non_numeric), tablefmt="plain"))



