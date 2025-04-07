import argparse
import os
import numpy as np
import gzip
import subprocess
from tqdm import tqdm

def subsample_and_write_fastq(fastq_ref_file, fastq_write_file, prob_each_line, proportion_to_subsample):
    """
    Subsamples the fastq file and writes the subsampled reads to the output file.
    """

    if fastq_ref_file.endswith(".gz"):
        open_func = gzip.open
        mode = "rt"
    else:
        open_func = open
        mode = "r"

    with open_func(fastq_ref_file, mode) as ref_file, open(fastq_write_file, "w") as write_file:
        pbar = tqdm(unit="read", desc="Subsampling FASTQ")

        for i, line in enumerate(ref_file):
            if i % 4 == 0:
                header = line
            elif i % 4 == 1:
                seq = line
            elif i % 4 == 2:
                plus = line
            elif i % 4 == 3:
                quality = line

                if prob_each_line[i // 4] < proportion_to_subsample:
                    write_file.write(header)
                    write_file.write(seq)
                    write_file.write(plus)
                    write_file.write(quality)
                else:
                    header = None
                    seq = None
                    plus = None
                    quality = None
                    continue
                # updates only if the read is written
                pbar.update(1)
        pbar.close()


Parser = argparse.ArgumentParser(description="Subsample fastq files to a given number of reads.")
Parser.add_argument("--fastq1", type=str, help="Path to the first fastq file.")
Parser.add_argument("--fastq2", type=str, help="Path to the second fastq file.")
Parser.add_argument("--out_dir", type=str, help="Output directory.")
Parser.add_argument("--n_reads_present", type=int, default=None, help="Output directory.")
Parser.add_argument("--n_reads", type=int, help="Number of reads to subsample.")

args = Parser.parse_args()
print(args)
fastq1 = args.fastq1
fastq2 = args.fastq2
out_dir = args.out_dir
n_reads_present = args.n_reads_present
n_reads_target = args.n_reads


fastq1_abs_path = os.path.abspath(fastq1)
fastq2_abs_path = os.path.abspath(fastq2)
out_dir_abs_path = os.path.abspath(out_dir)
print(f"Fastq1 absolute path: {fastq1_abs_path}")
print(f"Fastq2 absolute path: {fastq2_abs_path}")


if not os.path.exists(out_dir_abs_path):
    os.makedirs(out_dir_abs_path)

subsampled_dir = os.path.join(out_dir_abs_path, f"subsampled_{int(n_reads_target)}")
if not os.path.exists(subsampled_dir):
    os.makedirs(subsampled_dir)

read_1_out_path = os.path.join(subsampled_dir, f"subsampled_R1.fastq")
read_2_out_path = os.path.join(subsampled_dir, f"subsampled_R2.fastq")
print(f"Output directory: {subsampled_dir}")
print(f"Output file for read 1: {read_1_out_path}")
print(f"Output file for read 2: {read_2_out_path}")

# ensure fastq files exist
if not os.path.exists(fastq1_abs_path):
    raise ValueError(f"Fastq file {fastq1_abs_path} does not exist.")

if not os.path.exists(fastq2_abs_path):
    raise ValueError(f"Fastq file {fastq2_abs_path} does not exist.")

if n_reads_present is None:
    # get the number of reads in the first fastq file
    print(f"Counting number of reads in {fastq1_abs_path}")
    if fastq1_abs_path.endswith(".gz"):
        cmd = f"zcat {fastq1_abs_path} | wc -l"
    else:
        cmd = f"wc -l {fastq1_abs_path}"
    n_reads_present = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout.split()[0]
    n_reads_present = int(n_reads_present) // 4

print(f"Number of reads present in {fastq1_abs_path}: {n_reads_present}")

proportion_to_subsample = n_reads_target / n_reads_present
print(f"Subsampling {proportion_to_subsample:.2%} of reads from {fastq1_abs_path} and {fastq2_abs_path} to {n_reads_target} reads.")

np.random.seed(0)
prob_each_line = np.random.rand(n_reads_present)

print("sampling R1")
subsample_and_write_fastq(fastq1_abs_path, read_1_out_path, prob_each_line, proportion_to_subsample)

print("sampling R2")
subsample_and_write_fastq(fastq2_abs_path, read_2_out_path, prob_each_line, proportion_to_subsample)    

print("Done.")


