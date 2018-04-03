# Computational pipeline for analyzing data from the protein-stability assay

Below is a list of modules for computing EC50 values from deep-sequencing data

## `compute_ec50_values_from_NGS_data.py`

Compute EC50 values from deep-sequencing data. This module serves as a wrapper for the modules listed below.

* Inputs:
    * a CSV summary file with the metadata, where samples are rows and columns are different pieces of metadata. The columns must include:
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

    * the relative path between the current working directory and the directory that has the FASTQ files

* Dependencies:
    * all dependencies in the below programs.

* Outputs:
    * Assembled FASTQ files for each sample
    * Counts files for each sample, one giving all counts of all proteins, the other only giving counts for the input designs
    * EC50 values for each sample

## `assemble_and_align_reads.py`

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

## `summarize_alignments.py`

Make plots summarizing the alignments, including plots showing:
    * sequencing depth
    * number of sequences thrown out during the alignment
    * etc.

## `compute_ec50_values.py`

Compute ec50 values for a given experiment

* Inputs:
    * `experiments.csv` file

* Dependencies:
    * scripts from Rocklin et al., 2017, Science
    * `Python` modules

* Output:
    * EC50 values
