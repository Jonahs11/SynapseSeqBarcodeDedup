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
out_path = args.output
# if use_slurm ensure out_path is a directory
if use_slurm:
    # assert out_path is not a file
    assert not os.path.isfile(out_path), f"out_path {out_path} is a file, should be a directory"
    # create directory if it does not exist
    os.makedirs(out_path, exist_ok=True)


# get current file path
current_dir = os.path.dirname(os.path.abspath(__file__))
config_dir = os.path.dirname(current_dir)


config_yaml_path = os.path.join(config_dir, "config.yaml")
config = yaml.safe_load(open(config_yaml_path, "r"))

prepend_path = config.get("prepend_cmds_path", "")

# read in samplesheet
sample_df = pd.read_csv(args.samplesheet, sep=',')

cmds_list = []
exec_path = os.path.join(current_dir, "02_graph_dedup_bulk_pool.py")

for inx, row in sample_df.iterrows():
    out_dir_extract = row['out_dir_general']
    out_dir_dedup = row['out_dir_general']
    sample_id = row['sample_id']

    out_log_file = os.path.join(out_dir_dedup, sample_id, f"{sample_id}_dedup.log")    

    cmd = f'{prepend_path}python {exec_path} -b {out_dir_extract} -r {sample_id} -o {out_dir_dedup} > {out_log_file} 2>&1'
    cmd_tup = (sample_id, cmd)
    cmds_list.append(cmd_tup)

full_out_path = os.path.abspath(out_path)
if use_slurm:
    for sample_id, cmd in cmds_list:
        # create sbatch file
        sbatch_path = os.path.join(full_out_path, f"{sample_id}.sh")

        out_path = os.path.join(full_out_path, f"{sample_id}.out")
        err_path = os.path.join(full_out_path, f"{sample_id}.err")
        MEM = "100G"
        partition = "hpcx_macosko"
        with open(sbatch_path, "w") as f:
            f.write(f"#!/bin/bash\n")
            f.write(f"#SBATCH --job-name={sample_id}\n")
            f.write(f"#SBATCH --output={out_path}\n")
            f.write(f"#SBATCH --error={err_path}\n")
            f.write(f"#SBATCH --time=24:00:00\n")
            f.write(f"#SBATCH --cpus-per-task=3\n")
            f.write(f"#SBATCH --mem={MEM}\n")
            f.write(f"#SBATCH --nodes=1\n")
            f.write(f"#SBATCH --ntasks=1\n")
            if partition is not None:
                f.write(f"#SBATCH --partition={partition}\n")
            f.write(f"{cmd}\n")
        # run sbatch file
        # subprocess.run(f"sbatch {sbatch_path}", shell=True)
else:
    # write commands to file
    with open(out_path, "w") as f:
        for cmd in cmds_list:
            f.write(f"{cmd[1]}\n")



