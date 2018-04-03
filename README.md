# Computational pipeline for analyzing data from the high-throughput protein-stability assay

This directory contains a computational pipeline for analyzing data from the high-throughput protein-stability assay from [Rocklin et al, 2017, Science](https://doi.org/10.1126/science.aan0693).

## How to run the pipeline

The pipeline is in the form of a series of `Python` scripts that are contained in the directory called `scripts/`. The entire pipeline can be run by calling the script `compute_ec50_values_from_deep_sequencing_data.py` with the appropriate command-line arguments, which are described below:

* `designed_seqs_file` : the path to a CSV file giving the name and protein sequence for each input design. See [here](data/Rocklin_2017_Science/designed_protein_sequences.csv) for an example. In this file, each row specifies a design and each column specifies information about that design. This file must have the following columns:
    * `name` : a unique name for the design
    * `protein_sequence` : the protein sequence of the design

* `experimental_summary_file` : the path to a CSV file with experimental metadata, where samples are rows and columns are different pieces of metadata. See [here](data/Rocklin_2017_Science/experimental_summary.csv) for an example. This file is based on the `experiments.csv` file from the Rocklin et al. study, but contains come additional columns, and is also missing certain columns that will be added in later by the analysis code. The columns must include:
    * `experiment_id` : a unique ID for the entire experiment (**need to decide if I want to use this somehow**)
    * `protease_type` : the protease associated with the sample (e.g., "trypsin" or "chymotrypsin")
    * `selection_strength` : the index of the selection step associated with the sample, following the indexing scheme used in the Rocklin et al. study. The standard range of these values is normally: 0-6, where "0" corresponds to the naive library that has not been exposed to protease and "6" corresponds to the library that has been challenged with the highest concentration of protease.
    * `conc_factor` : the fold-change in protease concentration between selection steps, i.e., when `selection_strength` is incrimented by a value of one. A value of 3 would indicate that the protease concentration is increased by 3 fold between selection steps.
    * `parent` : the value of `selection_strength` for the sample that serves as the input library for the given selection. A value of 0 would indicate that the 
    * `fastq_id` : a string that is both common to all FASTQ files for a given sample and is *not* found in any other sample, such that this string can be used to search for and retrieve all relevant FASTQ files in the directory specified by the input argument called `fastq_dir` (see below).
    * `parent_expression`: the fraction of cells (events) passing the selection threshold in the given library before proteolysis, according to the sorting instrument.
    * `fraction_collected`: the fraction of cells (events) passing the selection threshold in the given library after proteolysis, according to the sorting instrument
    * `cells_collected`: the total number of cells (events) collected during the given selection, according to the sorting instrument

* `fastq_dir` : a path to the directory that has the FASTQ files

* `pare_path` : a path to the program `PEAR`


## An example analysis

I provide an example of how to execute this pipeline in the `Jupyter` notebook called `analysis_code.ipynb`. In this notebook, I reproduce the entire analysis from the Rocklin et al. study, starting from the deep-sequencing data. Most of the input data for the analysis is stored in the directory called `data/`. However, the input FASTQ files are stored in a separate location on TACC.

## Summary of `Python` scripts in the pipeline

All of these scripts are in the `scripts/` directory.

### `compute_ec50_values_from_deep_sequencing_data.py`

Compute EC50 values from deep-sequencing data. This module serves as a wrapper for the modules listed below.

* Outputs:
    * Assembled FASTQ files for each sample
    * Counts files for each sample, one giving all counts of all proteins, the other only giving counts for the input designs
    * EC50 values for each sample

### `assemble_and_align_reads.py`

Assemble paired-end deep-sequencing reads in the form of FASTQ files and then align the reads to a set of input sequences

* Inputs:
    * deep-sequencing data
    * a file giving the DNA and protein sequence for each input design for a given round of designs
    * a path to an output directory where the results will be stored

* Dependencies:
    * `PEAR`
    * `deep_seq_utils.py`: a custom `Python` script

* Outputs:
    * assembled paired-end reads in the form of FASTQ files for each sample
    * counts files with the number of times each protein was observed in each sample:
        * one file with counts for *all* sequences that were observed, even if they do not match a starting design
        * one file with counts for *only* sequences that match a starting design

### `summarize_alignments.py`

Make plots summarizing the alignments, including plots showing:
    * sequencing depth
    * number of sequences thrown out during the alignment
    * etc.

### `compute_ec50_values.py`

Compute ec50 values for a given experiment

* Inputs:
    * `experiments.csv` file

* Dependencies:
    * scripts from Rocklin et al., 2017, Science
    * `Python` modules

* Output:
    * EC50 values
