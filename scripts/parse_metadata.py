"""
Parse files from the BIOFAB with metadata from deep-sequencing and FACS experiments related to the 
"""

# Import `Python` modules
import os
import sys
import argparse
import re
import glob

import xml.etree.ElementTree as ET
sys.path.append("/home/jupyter/tacc-work/jupyter_packages/lib/python2.7/site-packages")
from FlowCytometryTools import *

import numpy as np
import pandas

def main():
    """Run the main code"""
    
    #---------------------------------------------------------------
    # Read in command-line arguments and experimental metadata
    #---------------------------------------------------------------
    # Read in command-line arguments using `argparse`    
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_dir", help="the path to a directory with input FASTQ files")
    parser.add_argument("--job_ids", help="specific Job_IDs specified by the BIOFAB")
    parser.add_argument("--facs_dirs", help="a comma-separated string of directories with input FCS files")
    parser.add_argument("--facs_channel_name", help="a string giving the name of the FACS channel to analyze (e.g, u'FITC-A'")
    parser.add_argument("--xml_files", help="a comma-separated string of the paths to XML files with additional FACS data")
    parser.add_argument("--library_name", help="a string giving the name of the library to analyze, corresponding to values in the `strain` column in the BIOFAB's metadata files")
    parser.add_argument("--ignore_aq_item_ids", help="a comma-separated string of aquarium IDs to ignore, corresponding to values in the `aq_item_id` column in the BIOFAB's metadata files. Put 'None' if you do not wish to ignore any, assuming it is not a valid `aq_item_id` for the experiment")
    parser.add_argument("--output_file", help="a path to an output file with all metadata aggregated into a single CSV")
    args = parser.parse_args()

    # Assign arguments to variables
    fastq_dir = args.fastq_dir
    job_ids = args.job_ids.split(',')
    facs_dirs = args.facs_dirs.split(',')
    facs_channel_name = args.facs_channel_name
    xml_files = args.xml_files.split(',')
    library_name = args.library_name
    ignore_aq_item_ids = args.ignore_aq_item_ids.split(',')
    output_file = args.output_file
    
    #
    print('hi')
    
    
if __name__ == "__main__":
    main()