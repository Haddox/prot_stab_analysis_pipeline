# Computational pipeline for analyzing data from the high-throughput protein-stability assay

This directory contains a computational pipeline for analyzing data from the high-throughput protein-stability assay from [Rocklin et al, 2017, Science](https://doi.org/10.1126/science.aan0693). Much of the code that makes up the pipeline is derived from the Rocklin et al. paper.


## Installing external dependencies

Carrying out the pipeline requires multiple external dependencies. I have created a [`Conda`](https://conda.io/docs/index.html) environment with nearly all required dependencies. This environment can be recreated using the [`environment.yml`](environment.yml) file, as described [here](https://conda.io/docs/user-guide/tasks/manage-environments.html) using the command:

    conda env create -f environment.yml

By default, the name of the environment will be `2018_prot_stab`. The pipeline must be run in an activated version of this environment. To activate the environment, use the command:

    source activate 2018_prot_stab
    
(or whatever environment name you chose instead of `2018_prot_stab`).

This pipeline also requires a program called [`PEAR`](https://sco.h-its.org/exelixis/web/software/pear/) to assemble paired-end deep-sequencing reads. Instructions for installing `PEAR` are provided on its [website](https://sco.h-its.org/exelixis/web/software/pear/).


## How to run the pipeline

The pipeline is in the form of a series of `Python` scripts that are contained in the directory called `scripts/`. The entire pipeline can be run by calling the script `compute_ec50_values_from_deep_sequencing_data.py`. If you are using a `Conda` environment (see above), you must first activate it using the command:

    source activate 2018_prot_stab
    
(or whatever environment name you chose instead of `2018_prot_stab`). Next, call the script `compute_ec50_values_from_deep_sequencing_data.py` with the appropriate command-line arguments, as described below. All arguments are required:

    python compute_ec50_values_from_deep_sequencing_data.py [--designed_sequences_file DESIGNED_SEQUENCES_FILE] [--experimental_summary_file EXPERIMENTAL_SUMMARY_FILE] [--fastq_dir FASTQ_DIR] [--pear_path PEAR_PATH] [--output_dir OUTPUT_DIR]

* `--designed_sequences_file` : the path to a CSV file giving the name and protein sequence for each input design. See [here](data/Rocklin_2017_Science/designed_protein_sequences.csv) for an example. In this file, each row specifies a design and each column specifies information about that design. This file must have the following comma-delimited columns:
    * `name` : a unique name for the design
    * `protein_sequence` : the protein sequence of the design

* `--fastq_dir` : a path to a directory that has the deep-sequencing data in the form of unassembled paired-end FASTQ files for all samples in the experiment.

* `--experimental_summary_file` : the path to a CSV file with a variety of metadata for each sample in the experiment, where samples are rows and columns are different pieces of metadata. See [here](data/Rocklin_2017_Science/experimental_summary.csv) for an example. This file is based on the `experiments.csv` file from the Rocklin et al. study, but contains come additional columns, and is also missing certain columns that will be added in later by the analysis code. The columns must be comma-delimited and must include:
    * `experiment_id` : a unique ID for the entire experiment (**need to decide if I want to use this somehow**)
    * `protease_type` : the protease associated with the sample (e.g., "trypsin" or "chymotrypsin")
    * `selection_strength` : the index of the selection step associated with the sample, following the indexing scheme used in the Rocklin et al. study. The standard range of these values is normally: 0-6, where "0" corresponds to the naive library that has not been exposed to protease and "6" corresponds to the library that has been challenged with the highest concentration of protease.
    * `conc_factor` : the fold-change in protease concentration between selection steps, i.e., when `selection_strength` is incrimented by a value of one. A value of 3 would indicate that the protease concentration is increased by 3 fold between selection steps. Note: leave this value blank for samples that have not been treated with any protease.
    * `parent` : the value of `selection_strength` for the sample that serves as the input library for the given selection. A value of 0 would indicate that the 
    * `fastq_id` : a string that is common to all FASTQ files for a given sample ***and*** is unique to those files. The code will search within the directory specified by `--fastq_dir` for all FASTQ files that contain this string in their name, aggregating the data among all matches. Thus, strings that are not unique to a particular dataset will result in "cross contamination" between datasets. Be careful when using numbers. For instance, the string "test1" would find not only FASTQ files with "test1" in their name, but also FASTQ files with "test10" in their name. Using a string like "test1\_" could be used to get around this problem.
    * `parent_expression`: the fraction of cells (events) passing the selection threshold in the given library before proteolysis, according to the sorting instrument.
    * `fraction_collected`: the fraction of cells (events) passing the selection threshold in the given library after proteolysis, according to the sorting instrument
    * `cells_collected`: the total number of cells (events) collected during the given selection, according to the sorting instrument

* `--pare_path` : a path to the program [`PEAR`](https://sco.h-its.org/exelixis/web/software/pear/)

* `--output_dir` : a path to an output directory where all the results will be stored. This directory will be created if it does not already exist


## Output of the pipeline

All results are stored in the directory specified by the input command `--output_dir`. Within this directory, there are multiple subdirectories that contain the results:

* `log.log`: a log file (currently, this file is not actually generated; I still need to implement this)

* `paired_FASTQ_files/`: this directory contains assembled FASTQ files and associated log files generated by `PARE`:
    * for each sample and for each pair of unassembled FASTQ files associated with that sample, a FASTQ of the assembled paired-end reads. The assembled FASTQ files are named according to the format: `{protease_type}_{selection_strength}-{n}.fastq`, where `protease_type` and `selection_strength` are derived from the corresponding entries in the input file specified by the command `--experimental_summary_file`, and where `n` corresponds to the index of the assembled pair of FASTQ files, indexed starting at one. This index is necessary since some samples may be associated with more than one pair of FASTQ files.

* `counts/`: this directory contains counts files, including:
    * for each sample, a file giving counts of all unique protein sequences observed in the deep-sequencing data. These files are named according to the format: `{protease_type}_{selection_strength}_counts.csv`.
    * for each protease, a file giving the counts aggregated across all samples that are associated with that protease. This file only includes counts for proteins that match one of the input designs included in the file specified by the input command `--designed_sequences_file` ***and*** have more than one counts in the naive library (selection index = 0). These files are named according to the format: `{protease_type}.counts`.
    
* `ec50_values/`: this directory contains a variety of output files, including:
    * `experiments.csv` : a file that serves as input into the script for computing EC50 values. This file is identical to the `experiments.csv` file described in the Rocklin et al. study.
    * for each protease, a file giving the aggregated counts across all samples that are associated with that protase. These files are the same as the ones in `counts/`, copied over for the purpose of inferring EC50s (the current script for inferring EC50s requires that all input is in the same directory; this should be changed in the future).
    * for each protease, a file called `{protease}.fulloutput` giving the EC50 values and other metadata in tab-delimited columns. These files are similar to the `.fulloutput` files from the Rocklin et al. study. See [here](data/original_Rocklin_EC50_values/rd4_chymo.sel_k0.8.erf.5e-7.0.0001.3cycles.fulloutput) for an example.
    
* `stability_scores/`: this directory contains a variety of output files, including:
    * for each protease, a file called `{protease}_stability_scores.txt`, which is the same as the corresponding `{protease}.fulloutput` file in the `ec50_values/` subdirectory, but with the following additional columns:
        * `ec50_pred`: a sequence's unfolded EC50, as predicted by the unfolded-state model from Rocklin et al.
        * `ec50_rise`: the difference between a sequence's observed EC50 value and the predicted EC50 value of the unfolded sequence (= observed - predicted unfolded)
        * `stabilityscore`: the stability score (**TODO: describe how this is computed**)


## An example analysis that reproduces the results from Rocklin et al. from the starting deep-sequencing and FACS data

I provide an example of how to execute this pipeline in the `Jupyter` notebook called [`example_analysis_code.ipynb`](example_analysis_code.ipynb). In this notebook, I reproduce the entire analysis from the Rocklin et al. study, starting from the deep-sequencing data. Most of the input data for the analysis is stored in the directory called `data/`. However, the input FASTQ files are stored in a separate location on TACC.

After executing the analysis, I test that the results of the pipeline match the original results from the Rocklin et al. study. To do so, I downloaded the following files from the paper's supplemental info and uploaded them in the following location in this repository:

* for each protease, I uploaded a file giving protein counts generated from the raw deep-sequencing data. These files are stored in the directory `data/original_Rocklin_counts/` and are called `rd4_chymo.counts` and `rd4_tryp.counts` for chymotrypsin and trypsin, respectively.
* for each protease, I uploaded a file giving EC50 values (and other related metadata) for each protein design. These files are stored in the directory `data/original_Rocklin_EC50_values/` and are called `rd4_chymo.sel_k0.8.erf.5e-7.0.0001.3cycles.fulloutput` and `rd4_tryp.sel_k0.8.erf.5e-7.0.0001.3cycles.fulloutput` for chymotrypsin and trypsin, respectively.


## Organization of data and computer code for analyzing new datasets

Below is a list of additional datasets analyzed in this project:

* `Inna_April_2016`: I generate an experimental summary file in the notebook: `create_experimental_summary_file_Inna_April_2016.ipynb` and then carry out the pipeline in the notebook: analysis_code_Inna_April_2016.ipynb


## Summary of `Python` scripts in the pipeline

All of the code for the pipeline is in the directory called `scripts/`. This includes scripts that perform the analyses, as well as scripts with functions that are imported as modules.

### `compute_ec50_values_from_deep_sequencing_data.py`:

This is the main script that performs the entire analysis.
* Inputs:
    * all inputs are described above in the section called "How to run the pipeline"
* Dependencies:
    * all scripts in the `scripts/` directory, described below
* Outputs:
    * all outputs are described above in the section called "How to run the pipeline"

### `fit_all_ec50_data.py`

This is a script from Rocklin et al. that is used to fit EC50 values. I call this script for this purpose in `compute_ec50_values_from_deep_sequencing_data.py`.

    python fit_all_ec50_data.py [--counts_dir COUNTS_DIR] [--experimental_summary_file EXPERIMENTAL_SUMMARY_FILE] [--datasets DATASETS] [--output_dir OUTPUT_DIR]

* Inputs:

    * `--counts_dir`: a path to the directory with the input counts files giving counts for a single protease across all selection levels (e.g., you might have one counts file for trypsin and one for chymotrypsin). Counts files must be structured such that rows are proteins, and columns (space-delimited) with the following names/info:
        * `name`: the name of the given protein
        * `counts{selection_index_0}`: the counts for a given protein at the first selection index in the experiment (e.g., this might be set to `counts0` for the counts in the naive library).
        * `counts{selection_index_1}`: the counts for a given protein at the second selection index in the experiment.
        * ... include a column for every selection index in the experiment (e.g., the Rocklin et al. study had columns for selection indices 0-6).
    * `--experimental_summary_file`: the path to an input file that follows the exact same format as the `experiments.csv` file from the Rocklin et al. study. This file has the exact same information as the `--experimental_summary_file` input file `compute_ec50_values_from_deep_sequencing_data.py`, but does not have the columns called `experiment_id` and `fastq_id`, and has additional columns called:
        * `input`: the name of the counts file with counts for a given sample, excluding the prefix provided in `--counts_dir` (e.g., "trypsin.counts" might have counts for all samples challenged with trypsin)
        * `column`: the name of the the column in the counts file that corresponds to the selection level of a given sample (e.g., "counts1" would correspond to the first selection level).
    * `--datasets`: a string of comma-delimited datasets to analyze (no-spaces between items, only commas). These datasets are defined in the `input` column of the `experiments.csv` file, formatted as `{dataset}.counts`
    * `--output_dir`: a path to an output directory where all the results will be stored. This directory will be made if it does not already exist.

* Dependencies:
    * the modules called `protease_sequencing_model.py`, `utility.py`, and `compile_counts_and_FACS_data/__init__.py`, which are described below
        
* Outputs:
    * all files in the `ec50_values/` results directory described above, except for the `experiments.csv` file

### `compute_stability_scores_from_EC50_values.py`

This is a script that computes stability scores from EC50 values using the unfolded-state model from Rocklin et al.

* Inputs:
    * need to add
* Dependencies:
    * the module called `sequence_protease_susceptibility.py`, which is described below
    * the parameters file called `unfolded_state_model_params`, which is described below
* Outputs:
    * all files in the `stability_scores/` results directory described above


### `deep_seq_utils.py`
This is a custom script with `Python` functions for analyzing deep-sequencing data, used in `compute_ec50_values_from_deep_sequencing_data.py`.

### `protease_sequencing_model.py` and `utility.py`
These are both scripts from Rocklin et al. that are imported in `fit_all_ec50_data.py` as modules for computing EC50 values.


### `sequence_protease_susceptibility.py`
This script has functionalities that are used to implement the unfolded-state model in the script `compute_stability_scores_from_EC50_values.py`.

### `unfolded_state_model_params`
A file with values used to parameterize the unfolded-state model in the script `compute_stability_scores_from_EC50_values.py`.


## To do

* Set up a system for logging progress of the script
* Figure out how to implement the unfolded-state model and compute EC50 values
* Set up a `conda` environment that is completely self contained
    * currently, there is a problem with importing `pymc3` as installed by `conda`
* Ask Gabe if he included samples with zero counts in the naive sample.