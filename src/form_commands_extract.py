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
parser.add_argument('--output', type=str, help='output file')
# set as store true if passed
parser.add_argument('--use_slurm', action='store_true', help='use slurm')
args = parser.parse_args()

# IF using slurm, will write sbatch files to out_path, need to store this in a directory
use_slurm = args.use_slurm
out_path = os.path.abspath(args.output)

# if use_slurm ensure out_path is a directory
if use_slurm:

    MEM = "200G"
    partition="hpcx_macosko"
    # assert out_path is not a file
    assert not os.path.isfile(out_path), f"out_path {out_path} is a file, should be a directory"
    # create directory if it does not exist
    os.makedirs(out_path, exist_ok=True)

# get dirname and create dir
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# get current file path
current_dir = os.path.dirname(os.path.abspath(__file__))
config_dir = os.path.dirname(current_dir)

config_yaml_path = os.path.join(config_dir, "config.yaml")
config = yaml.safe_load(open(config_yaml_path, "r"))

# bulk_parse_fastq_path = os.path.join(base_ss_dir, "02_bulk", "extract_barcodes_bulk.py")
bulk_parse_fastq_path = os.path.join(current_dir, "extract_barcodes.py")
assert os.path.exists(bulk_parse_fastq_path), f"bulk_parse_fastq_path {bulk_parse_fastq_path} does not exist"

prepend_path = config.get("prepend_cmds_path", "")

# read in samplesheet
samplesheet_path_full = os.path.abspath(args.samplesheet)
sample_df = pd.read_csv(samplesheet_path_full, sep=',')

cmds_list = []

for inx, row in sample_df.iterrows():
    fastq_dir = row['fastq_dir']
    # read full path for r1 and r2
    # list files in fastq_dir and filter for r1 and r2
    r1 = [f for f in os.listdir(fastq_dir) if "R1" in f][0]
    r2 = [f for f in os.listdir(fastq_dir) if "R2" in f][0]
    r1_path = os.path.join(fastq_dir, r1)
    r2_path = os.path.join(fastq_dir, r2)

    assert os.path.exists(r1_path), f"r1 file {r1_path} does not exist"
    assert os.path.exists(r2_path), f"r2 file {r2_path} does not exist"

    out_dir_extract = row['out_dir_general']
    sample_id = row['sample_id']
    out_dir = os.path.join(out_dir_extract, sample_id)
    os.makedirs(out_dir, exist_ok=True)

    out_log_file = os.path.join(out_dir, f"{sample_id}_extract.log")    

    # create command
    cmd = f"python -u {bulk_parse_fastq_path} {r1_path} {r2_path} --output-dir {out_dir} > {out_log_file} 2>&1"
    if len(prepend_path) > 0:
        cmd = f"{prepend_path} {cmd}"
    cmd_tup = (sample_id, cmd)
    cmds_list.append(cmd_tup)

print(cmds_list)
if use_slurm:
    for sample_id, cmd in cmds_list:
        print(sample_id)
        # create sbatch file
        sbatch_path = os.path.join(out_path, f"{sample_id}.sh")
        print(sbatch_path)
        stdout_path = os.path.join(out_path, f"{sample_id}.out")
        err_path = os.path.join(out_path, f"{sample_id}.err")

        with open(sbatch_path, "w") as f:
            f.write(f"#!/bin/bash\n")
            f.write(f"#SBATCH --job-name={sample_id}\n")
            f.write(f"#SBATCH --output={stdout_path}\n")
            f.write(f"#SBATCH --error={err_path}\n")
            f.write(f"#SBATCH --time=48:00:00\n")
            f.write(f"#SBATCH --cpus-per-task=3\n")
            f.write(f"#SBATCH --mem={MEM}\n")
            f.write(f"#SBATCH --nodes=1\n")
            f.write(f"#SBATCH --ntasks=1\n")
            if partition is not None:
                f.write(f"#SBATCH --partition={partition}\n")
            f.write(f"{cmd}\n")

else:
    # write commands to file
    with open(out_path, "w") as f:
        for cmd in cmds_list:
            f.write(f"{cmd[1]}\n")
