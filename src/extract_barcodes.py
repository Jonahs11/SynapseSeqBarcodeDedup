import gzip
import itertools
import logging
from collections import Counter, defaultdict
from pathlib import Path
import pickle
import click
import networkx as nx
import numpy as np
import scipy.io
import scipy.sparse
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.neighbors import radius_neighbors_graph

import yaml

import subprocess
import os

print("Begining Barcode Extraction")

log = logging.getLogger("barcode_matrix")

# load config
current_dir = os.path.dirname(os.path.abspath(__file__))
config_dir = os.path.dirname(current_dir)

config_yaml_path = os.path.join(config_dir, "config.yaml")
config = yaml.safe_load(open(config_yaml_path, "r"))


ug_exec_path = config["ugrep_path"]
vt_start_index_r2 = config["vt_start_index_r2"]
umi_start_pos_r1 = config["umi_start_pos_r1_1_index"]
umi_end_pos_r1 = config["umi_end_pos_r1_1_index"]

# assert umi start and end are integers
assert isinstance(umi_start_pos_r1, int), "umi_start_pos_r1 should be an integer"
assert isinstance(umi_end_pos_r1, int), "umi_end_pos_r1 should be an integer"
# assert umi start and end are positive
assert umi_start_pos_r1 > 0, "umi_start_pos_r1 should be positive"
assert umi_end_pos_r1 > 0, "umi_end_pos_r1 should be positive"
assert umi_start_pos_r1 < umi_end_pos_r1, "umi_start_pos_r1 should be less than umi_end_pos_r1"

# WSBWSDWSHWSVWSBWSWSHWSVWSBWSDWSH
DEGENERATE_BASE_DICT = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "M": {"A", "C"},
    "K": {"G", "T"},
    "S": {"C", "G"},
    "W": {"A", "T"},
    "H": {"A", "C", "T"},
    "B": {"C", "G", "T"},
    "V": {"A", "C", "G"},
    "D": {"A", "G", "T"},
    "N": {"A", "C", "G", "T"},
}

REVERSE_COMP_DICT = str.maketrans('ACGTN', 'TGCAN')

def match_tags(
    cell_ID, umis, barcode_codes, barcode_cutoff
):
    no_id_struct = 0
    pass_all = 0

    n_bases_bc = len(cell_ID[0])
    
    raw_umis_per_tag = defaultdict(set)
    raw_reads_per_umi = defaultdict(Counter)

    distribution_of_barcode_errors = {}
    for tag, umi in zip(cell_ID, umis):

        # Skip if Cell ID does not match structure
        
        num_matches = sum(b in s for b, s in zip(tag, barcode_codes))
        num_errors = n_bases_bc - num_matches
        current_num_w_error_count = distribution_of_barcode_errors.get(num_errors, 0)
        current_num_w_error_count+=1
        distribution_of_barcode_errors[num_errors] = current_num_w_error_count


        if num_matches < barcode_cutoff:
            no_id_struct += 1
            continue

        pass_all += 1

        raw_umis_per_tag[tag].add(umi)
        raw_reads_per_umi[tag][umi] += 1

    log.debug(f"IDs w/o correct ID structure: {no_id_struct}")
    log.debug(f"Passed all: {pass_all}")

    raw_umis_per_tag = {
        tag_seq: len(umi_seq) for tag_seq, umi_seq in raw_umis_per_tag.items()
    }

    return raw_umis_per_tag, raw_reads_per_umi, distribution_of_barcode_errors, pass_all


def write_matrix(upb, beads, tags, output_file):
    m = scipy.sparse.dok_matrix((len(beads), len(tags)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(beads)}
    t2j = {t: j for j, t in enumerate(tags)}

    for b, bd in upb.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = v

    m = m.tocsr()

    with gzip.open(output_file, "wb") as out:
        scipy.io.mmwrite(out, m)

    return m

def write_matrix_1d(upb, beads, output_file):
    m = scipy.sparse.dok_matrix((len(beads)), dtype=np.int32)
    b2i = {b: i for i, b in enumerate(beads)}
    t2j = {t: j for j, t in enumerate(tags)}

    for b, bd in upb.items():
        for t, v in bd.items():
            m[b2i[b], t2j[t]] = v

    m = m.tocsr()

    with gzip.open(output_file, "wb") as out:
        scipy.io.mmwrite(out, m)

    return m

def reverse_h_set(h_barcodes):
    base_d_rev = {0: "A", 1: "C", 2: "G", 3: "T", 4: "N"}

    return ["".join([ base_d_rev[c] for c in h_list ]) for h_list in h_barcodes]


# WSBWSDWSHWSVWSBWSDWSHWSVWSBWSDWSH
@click.command()
@click.argument("fastq-r1", type=click.Path(exists=True))
@click.argument("fastq-r2", type=click.Path(exists=True))
@click.option(
    "--output-dir",
    help="Path for output",
    required=True,
)
@click.option(
    "--tag-sequence",
    default="WSBWSDWSHWSVWSBWSWSHWSVWSBWSDWSH",
    help="Format for the tag sequence",
)
@click.option(
    "--constant-sequence-pretag",
    default="TGCCACCCGTAGATCTCTCGA",
    help="Constant sequence that should precede the tag",
)
@click.option(
    "--constant-sequence-posttag",
    default="CGATACCGAGCGCTGCACCGG",
    help="Constant sequence that should follow the tag",
)
@click.option(
    "--tag-mismatch",
    type=int,
    default=3,
    help="Number of mismatches to allow in tag",
    show_default=True,
)
@click.option(
    "--const-mismatch",
    type=int,
    default=3,
    help="Number of mismatches to allow in tag",
    show_default=True,
)
@click.option(
    "--only_summary",
    type=bool,
    default=False,
    show_default=True,
)
@click.option("--debug", is_flag=True, help="Turn on debug logging")
def main(
    fastq_r1,
    fastq_r2,
    output_dir,
    tag_sequence,
    constant_sequence_pretag,
    constant_sequence_posttag,
    tag_mismatch,
    const_mismatch,
    only_summary=False,
    debug=False,
):
    """
    This script generates some plots for mapping barcoded reads.

    Reads sequences from FASTQ_R1 and FASTQ_R2. Assumes that the first read
    contains a 15bp barcode split across two locations, along with an 8bp UMI.
    The second read is assumed to have TAG_SEQUENCE in bases 20-40.
    """
    # create_logger(debug, dryrun=False)

    print(len(tag_sequence))
    print(tag_sequence)
    
    # make output directory recursively
    os.makedirs(output_dir, exist_ok=True)

    output_dir = Path(output_dir)
    print(f"Saving output to {output_dir}")

    command = [f'zcat {fastq_r1} | head -n 2 | tail -n 1 | wc -c']
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    rd1_length = result.stdout.strip()

    rd1_length=int(rd1_length)
    a=rd1_length-1
    print(f'Read 1 is {a} nt long')

    preform_processing = not only_summary
    print(preform_processing)

    with_gzip = False
    merged_fname = "merged.fastq.gz"
    if not with_gzip:
        merged_fname = "merged.fastq"

    if (not (output_dir / merged_fname).exists()) and preform_processing:
        print("Merging sequence lines from read 1 and read 2 in one step")
        command = f'paste <(zcat {fastq_r1} | awk \'NR%4==2\') <(zcat {fastq_r2} | awk \'NR%4==2\')'
        if with_gzip:
            command += ' | gzip > merged.fastq.gz'
        else:
            command += ' > merged.fastq'
        print(command)
        # command = f'paste <(zcat {fastq_r1} | awk \'NR%4==2\') <(zcat {fastq_r2} | awk \'NR%4==2\') | gzip > merged.fastq.gz'
        subprocess.run(command, cwd=output_dir, shell=True, executable='/bin/bash')
    else:
        print("Merged file already exists. Skipping")


    filt_fname = "merged_preFilt_postFilt.fastq.gz"
    if not with_gzip:
        filt_fname = "merged_preFilt_postFilt.fastq"

    if (not (output_dir / filt_fname).exists()) &  preform_processing:
        print("Fuzzy grep-ing for polyA constant sequence")
        print(constant_sequence_pretag)
        command = f'cat {merged_fname} | {ug_exec_path} -Z3 {constant_sequence_pretag} |  {ug_exec_path} -Z3 {constant_sequence_posttag}'
        if with_gzip:
            command = 'z' + command
            command += f' | gzip > {filt_fname}'
        else:
            command += f' > {filt_fname}'
        print(command)
        result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True, text=True)

    bc_size = len(tag_sequence)
    
    bc_start = rd1_length + vt_start_index_r2
    bc_end = bc_start+bc_size-1
  
    print("Cutting out the identifiers/UMIs")

    if preform_processing:
        if with_gzip:
            command = [f'zcat {filt_fname} | cut -c {umi_start_pos_r1}-{umi_end_pos_r1} | gzip > umi_preFilt_postFilt.fastq.gz; zcat {filt_fname} | cut -c {bc_start}-{bc_end} | gzip > ident_preFilt_postFilt.fastq.gz']
        else:
            command = [f'cat {filt_fname} | cut -c {umi_start_pos_r1}-{umi_end_pos_r1} > umi_preFilt_postFilt.fastq; cat {filt_fname} | cut -c {bc_start}-{bc_end} > ident_preFilt_postFilt.fastq']
        print(command)
        result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True, text=True)

    print(f"Reading identifiers/UMIs")
    if with_gzip:
        cellID_path = output_dir / 'ident_preFilt_postFilt.fastq.gz'
        umi_path = output_dir / 'umi_preFilt_postFilt.fastq.gz'
        with gzip.open(cellID_path, "rt") as fh:
            cell_ID = [line.strip() for line in itertools.islice(fh, 0, None, 1)]
        with gzip.open(umi_path, "rt") as fh:
            umis = [line.strip() for line in itertools.islice(fh, 0, None, 1)]
    else:
        cellID_path = output_dir / 'ident_preFilt_postFilt.fastq'
        umi_path = output_dir / 'umi_preFilt_postFilt.fastq'
        with open(cellID_path, "r") as fh:
            cell_ID = [line.strip() for line in itertools.islice(fh, 0, None, 1)]
        with open(umi_path, "r") as fh:
            umis = [line.strip() for line in itertools.islice(fh, 0, None, 1)]

    assert len(cell_ID) == len(umis), "read different number of reads"
    num_reads_post_const_grep = len(cell_ID)
    
    print(f"Total of {len(cell_ID)} reads")

    barcode_codes = [DEGENERATE_BASE_DICT[b] for b in tag_sequence]
    barcode_cutoff = len(tag_sequence) - tag_mismatch

    print(f"barcode_cutoff: {barcode_cutoff}")

    raw_umis_per_tag, raw_reads_per_umi, distribution_of_barcode_errors, num_reads_passing = match_tags(
        cell_ID, umis, barcode_codes, barcode_cutoff
    )

    print("Writing output files")
    raw_tags = sorted(raw_umis_per_tag.keys())
    umis_count = [str(raw_umis_per_tag[t]) for t in raw_tags]
    raw_umis = sorted({u for t in raw_tags for u in raw_reads_per_umi[t]})

    print(f"{len(raw_tags)} unique bead IDs")

    with gzip.open(output_dir / "raw_tags.txt.gz", "wt") as out:
        print("\n".join(raw_tags), file=out)
    with gzip.open(output_dir / "umis_per_tags.txt.gz", "wt") as out:
        print("\n".join(umis_count), file=out)
    with gzip.open(output_dir / "raw_umis.txt.gz", "wt") as out:
        print("\n".join(raw_umis), file=out)

    log.debug("Writing raw umi/reads matrix")
    write_matrix(raw_reads_per_umi, raw_tags, raw_umis, output_dir / "raw_read_umi_matrix.mtx.gz")

    print("saving filtering metadata")
    filter_results = {
        "raw_nreads": n_lines_total,
        "post_const_grep_nreads": num_reads_post_const_grep,
        "post_barcode_structure_nreads": num_reads_passing,
        "distribution_matches_struc_errors": distribution_of_barcode_errors
    }

    meta_out_path = os.path.join(output_dir, "bc_extraction_meta.pkl")
    pickle.dump(filter_results, open(meta_out_path, "wb"))

    # remove intermediate files
    print("Removing intermediate files")
    command = [f'rm rd1_stripped.fastq rd2_stripped.fastq']
    result = subprocess.run(command, cwd=output_dir, shell=True, capture_output=True)
    
    print("Done!")
    
if __name__ == "__main__":
    main()

