


# read in args
import argparse
import os


Parser = argparse.ArgumentParser(description="Subsample fastq files to a given number of reads.")
Parser.add_argument("--fastq1", type=str, help="Path to the first fastq file.")
Parser.add_argument("--fastq2", type=str, help="Path to the second fastq file.")
Parser.add_argument("--out_file", type=str, help="Output directory.")
Parser.add_argument("--out_dir", type=str, help="Output directory.")
Parser.add_argument("--prefix", type=str, default="subsampled_test")


args = Parser.parse_args()
print(args)
fastq1 = args.fastq1
fastq2 = args.fastq2
out_file = args.out_file
out_dir = args.out_dir
prefix = args.prefix

print(f"fastq1: {fastq1}")
print(f"fastq2: {fastq2}")
print(f"out_file: {out_file}")
print(f"out_dir: {out_dir}")
print(f"prefix: {prefix}")

abs_path_fastq1 = os.path.abspath(fastq1)
abs_path_fastq2 = os.path.abspath(fastq2)
abs_path_out_file = os.path.abspath(out_file)
abs_path_out_dir = os.path.abspath(out_dir)
print(f"abs_path_fastq1: {abs_path_fastq1}")
print(f"abs_path_fastq2: {abs_path_fastq2}")
print(f"abs_path_out_file: {abs_path_out_file}")
print(f"abs_path_out_dir: {abs_path_out_dir}")

assert os.path.exists(abs_path_fastq1), f"File {abs_path_fastq1} does not exist."
assert os.path.exists(abs_path_fastq2), f"File {abs_path_fastq2} does not exist."

# check if out_dir exists, if not create it
if not os.path.exists(abs_path_out_dir):
    os.makedirs(abs_path_out_dir)




num_reads_sample = [5e7, 1e8, 5e8]

# for number of subsample reads, create folder for output, and write commands to subsample fastq files using seqtk to out_file
for num_reads in num_reads_sample:
    out_dir_num = os.path.join(out_dir, f"{prefix}_{int(num_reads)}")
    if not os.path.exists(out_dir_num):
        os.makedirs(out_dir_num)
    read1_out_path = os.path.join(out_dir_num, "subsampled_R1.fastq")
    read2_out_path = os.path.join(out_dir_num, "subsampled_R2.fastq")

    cmd1 = f"seqtk sample -s100 {abs_path_fastq1} {int(num_reads)} > {read1_out_path}"
    cmd2 = f"seqtk sample -s100 {abs_path_fastq2} {int(num_reads)} > {read2_out_path}"
    with open(abs_path_out_file, "a") as f:
        f.write(cmd1 + "\n")
        f.write(cmd2 + "\n")



