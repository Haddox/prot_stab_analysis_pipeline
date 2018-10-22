# Computational pipeline for analyzing data from the high-throughput protein-stability assay

This directory contains a computational pipeline for analyzing data from the high-throughput protein-stability assay from [Rocklin et al, 2017, Science](https://doi.org/10.1126/science.aan0693). Much of the code that makes up the pipeline is derived from the Rocklin et al. paper. This pipeline involves three main steps:

1) Computing EC50 values from FACS and deep-sequencing data using the script `scripts/compute_ec50_values_from_deep_sequencing_data.py`.
2) Computing stability scores from EC50 values using the script `scripts/compute_stability_scores_from_EC50_values.py`
3) Making summary plots using the script `scripts/create_summary_plots.py`.

Below are detailed instructions for conducting each step. Note: I've had to separate these into different steps since steps 1 and 2 require slightly different external dependencies (see the next section). In the future, my goal is to find a single set of external dependencies that works for everything, but this has been challenging.


## Installing external dependencies

Carrying out the pipeline requires external dependencies. As mentioned above, some steps require different external dependencies. For each step, I have encoded nearly all dependencies in [`Conda`](https://conda.io/docs/index.html) environments, which can be installed using the following YML files:

* `environment_compute_ec50_values.yml`: An environment for carrying out steps 1 and 3 of the pipeline.

* `environment_compute_stability_scores.yml`: An environment for carrying out step 2 of the pipeline.

as described [here](https://conda.io/docs/user-guide/tasks/manage-environments.html) using the command:

    conda env create -f {environment.yml}

where `environment.yml` is the name of the YML file. To run each step in the pipeline, the appropriate environment must first be activated. To do so, simply use the command:

    source activate {env_name}

where `env_name` is the name of the environment as encoded in the YML file (e.g., `2018_prot_stab_compute_stability_scores`).

In addition to the dependencies encoded in the `Conda` environments, this pipeline also requires a program called [`PEAR`](https://sco.h-its.org/exelixis/web/software/pear/) to assemble paired-end deep-sequencing reads. Instructions for installing `PEAR` are provided on its [website](https://sco.h-its.org/exelixis/web/software/pear/).


## How to run the pipeline

The pipeline is in the form of a series of `Python` scripts that are contained in the directory called `scripts/`. Below, I describe how to carry out each of the three steps of the pipeline summarized above.

### Step #1: Compute EC50 values from FACS and deep-sequencing data using the script `compute_ec50_values_from_deep_sequencing_data.py`

First, install the dependencies in the `Conda` environment encoded in `environment_compute_ec50_values.yml` (see above). Then activate it using the command:

    source activate 2018_prot_stab_compute_ec50_values

(`2018_prot_stab_compute_ec50_values` is the default name of the environment in the corresponding YML file).

Next, call the script `compute_ec50_values_from_deep_sequencing_data.py` with the appropriate command-line arguments, as described below. All arguments are required:

    python compute_ec50_values_from_deep_sequencing_data.py [--designed_sequences_file DESIGNED_SEQUENCES_FILE] [--experimental_summary_file EXPERIMENTAL_SUMMARY_FILE] [--fastq_dir FASTQ_DIR] [--pear_path PEAR_PATH] [--five_prime_flanking_seq FIVE_PRIME_FLANKING_SEQ] [--three_prime_flanking_seq THREE_PRIME_FLANKING_SEQ]  [--output_dir OUTPUT_DIR]

* `--designed_sequences_file` : the path to a CSV file giving the name and DNA or protein sequence for each input design. See [here](data/Rocklin_2017_Science/designed_protein_sequences.csv) for an example. In this file, each row specifies a design and each column specifies information about that design. This file must have the following comma-delimited columns:
    * `name` : a unique name for the design
    * `protein_sequence` or `dna_sequence` : the protein or DNA sequence of the design, respectively

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

* `--five_prime_flanking_seq`: a DNA sequence that flanks the coding sequence of interest on the 5' end of the sequencing read (string). The coding sequence should begin immediately after the final nucleotide of this flanking sequence. This flanking sequence and the one given by `three_prime_flanking_seq` will be used to extract the DNA coding sequence from each sequencing read and then translate it into a protein dsequence. Note: the default sequence used in Rocklin et al., 2017, Science was `CATATG`.

* `--three_prime_flanking_seq`: a DNA sequence that flanks the coding sequence of interest on the 3' end of the sequencing read (string). The coding sequence should begin immediately before the first nucleotide of this flanking sequence. Note: this DNA sequence should be in the same 5'-to-3' orientation as `five_prime_flanking_seq`. Note: the default sequence used in Rocklin et al., 2017, Science was `CTCGAG`

* `--output_dir`: a path to an output directory where all the results will be stored. This directory will be created if it does not already exist

* `--protein_or_dna_level`: whether `designed_sequences_file` contains DNA or protein sequences. Options are: 'protein' or 'dna'. Default is: 'protein'.

The output of this step is described in the below section, and includes all output except for the output data in the folder called `stability_scores`.


### Step #2: Compute stability scores from EC50 values using the script `compute_stability_scores_from_EC50_values.py`

First, install the dependencies in the `Conda` environment encoded in `environment_compute_stability_scores.yml` (see above). Then activate it using the command:

    source activate 2018_prot_stab_compute_stability_scores

(`2018_prot_stab_compute_stability_scores` is the default name of the environment in the corresponding YML file).

Next, call the script `compute_stability_scores_from_EC50_values.py` using the following command-line arguments:


### Step #3: Make summary plots using the script `create_summary_plots.py`

First, activate the same `Conda` environment used in step 1.

Next, call the script `create_summary_plots.py` using the following command-line argument:

    python scripts/create_summary_plots.py [--data_dir DATA_DIR] [--output_dir OUTPUT_DIR]

* `--data_dir`: the directory that contains all output data from steps 1 and 2 (specified above by `--output_dir`)
* `--output_dir`: the output directory where the plots will be stored


## Output of the pipeline

All results are stored in the directory specified by the input command `--output_dir`. Within this directory, there are multiple subdirectories that contain the results:

* `log.log`: a log file (currently, this file is not actually generated; I still need to implement this)

* `paired_FASTQ_files/`: this directory contains assembled FASTQ files and associated log files generated by `PARE`:
    * for each sample and for each pair of unassembled FASTQ files associated with that sample, a FASTQ of the assembled paired-end reads. The assembled FASTQ files are named according to the format: `{protease_type}_{selection_strength}-{n}.fastq`, where `protease_type` and `selection_strength` are derived from the corresponding entries in the input file specified by the command `--experimental_summary_file`, and where `n` corresponds to the index of the assembled pair of FASTQ files, indexed starting at one. This index is necessary since some samples may be associated with more than one pair of FASTQ files.

* `counts/`: this directory contains counts files, including:
    * for each sample, a file giving counts of all unique sequences observed in the deep-sequencing data. These files are named according to the format: `{protease_type}_{selection_strength}_counts.csv`. They contain either protein or DNA counts depending on the value of the `--protein_or_dna_level` input argument.
    * for each protease, a file giving the counts aggregated across all samples that are associated with that protease. This file only includes counts for sequences that match one of the input designs included in the file specified by the input command `--designed_sequences_file` ***and*** have more than one counts in the naive library (selection index = 0). These files are named according to the format: `{protease_type}.counts`. These files are similar to the `.fulloutput` files from the Rocklin et al. study. They contain either protein or DNA counts depending on the value of the `--protein_or_dna_level` input argument.

* `ec50_values/`: this directory contains a variety of output files, including:
    * `experiments.csv` : a file that serves as input into the script for computing EC50 values. This file is identical to the `experiments.csv` file described in the Rocklin et al. study.
    * for each protease, a file giving the aggregated counts across all samples that are associated with that protease. These files are the same as the ones in `counts/`, copied over for the purpose of inferring EC50s (the current script for inferring EC50s requires that all input is in the same directory; this should be changed in the future).
    * for each protease, a file called `{protease}.fulloutput` giving the EC50 values and other metadata in tab-delimited columns. These files are similar to the `.fulloutput` files from the Rocklin et al. study. See [here](data/original_Rocklin_EC50_values/rd4_chymo.sel_k0.8.erf.5e-7.0.0001.3cycles.fulloutput) for an example.

* `stability_scores/`: this directory contains a variety of output files, including:
    * for each protease, a file called `{protease}_stability_scores.txt`, which is the same as the corresponding `{protease}.fulloutput` file in the `ec50_values/` subdirectory, but with the following additional columns:
        * `ec50_pred`: a sequence's unfolded EC50, as predicted by the unfolded-state model from Rocklin et al.
        * `ec50_rise`: the difference between a sequence's observed EC50 value and the predicted EC50 value of the unfolded sequence (= observed - predicted unfolded)
        * `stabilityscore`: the stability score (**TODO: describe how this is computed**)


## A template notebook that compiles experimental data from the BIOFAB uploaded on TACC

The notebook [`template_notebook.ipynb`](template_notebook.ipynb) is a template for carrying out the above pipeline for computing stability scores starting from experimental data from the UW BIOFAB uploaded to TACC. See the instructions in the notebook for further details.

I ran the above notebook on TACC. To do so, you will need a TACC account. If you do not already have one, you can request one by following the instructions here: https://sd2e.org/accounts/request-access/

The way I ran this notebook on TACC was by logging onto the Maverick server using the following command:

    ssh [tacc_user_name]@maverick.tacc.utexas.edu

and then launching Jupyter using the arguments

    module use /work/03076/gzynda/public/apps/modulefiles
    module load singularity-sd2e
    submit_notebook

The last argument should prompt you to enter your email address. After doing so, you should receive an email. Use the URL and the password to launch Jupyter. Then change directories into `tacc_work/`, clone this repository, and run the notebook.


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
    * all scripts in the `scripts/` directory, described below, as well as `PARE`
* Outputs:
    * all outputs are described above in the section called "Output of the pipeline"

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
* Set up a `conda` environment that is completely self contained
    * currently, there is a problem with importing `pymc3` as installed by `conda`
* Ask Gabe if he included samples with zero counts in the naive sample.
