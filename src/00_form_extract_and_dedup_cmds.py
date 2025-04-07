import pandas as pd
import numpy as np
import os
import subprocess
import logging
import argparse
import yaml

# read in arguments from command line
parser = argparse.ArgumentParser(description='Create pool reference')
parser.add_argument('--samplesheet', type=str, help='input file')
parser.add_argument('--output_base', type=str, help='output file')
# set as store true if passed
parser.add_argument('--use_slurm', action='store_true', help='use slurm')
args = parser.parse_args()

sample_sheet_path = args.samplesheet
output_base = args.output_base
use_slurm = args.use_slurm


current_dir = os.path.dirname(os.path.abspath(__file__))

# form extract cmds
extract_out_path = f"{output_base}_extract_cmds.txt"
extract_wrapper_file = os.path.join(current_dir, "form_commands_extract.py")
cmd = f"python {extract_wrapper_file} --samplesheet {sample_sheet_path} --output {extract_out_path}"
if use_slurm:
    cmd += " --use_slurm"
logging.info(f"Running command: {cmd}")
subprocess.run(cmd, shell=True, check=True)

# form dedup cmds
dedup_out_path = f"{output_base}_dedup_cmds.txt"
dedup_wrapper_file = os.path.join(current_dir, "form_commands_dedup_bulk.py")

cmd = f"python {dedup_wrapper_file} --samplesheet {sample_sheet_path} --output {dedup_out_path}"
if use_slurm:
    cmd += " --use_slurm"
logging.info(f"Running command: {cmd}")
subprocess.run(cmd, shell=True, check=True)


