import sys
from pathlib import Path
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import shutil


def get_options():
    """
    Parse command-line arguments.
    """
    parser = ArgumentParser()

    # Directory containing original FASTA files
    parser.add_argument(
        "-i", "--input_data",
        type=Path,
        required=True,
        help="Directory containing original FASTA files"
    )

    # Directory where copied & modified FASTA files will be written
    parser.add_argument(
        "-o", "--output_dir",
        type=Path,
        required=True,
        help="Directory for copied and modified FASTA files"
    )

    # Root directory containing per-sample cluster subfolders
    parser.add_argument(
        "-c", "--cluster_root",
        type=Path,
        required=True,
        help="Root directory containing per-sample cluster subfolders"
    )

    return parser.parse_args()


def copy_input_files(input_dir: Path, output_dir: Path):
    """
    Copy only .pep and .cds files from input directory to output directory.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    for file_path in input_dir.iterdir():
        if file_path.is_file() and file_path.suffix in {".pep", ".cds"}:
            shutil.copy2(file_path, output_dir / file_path.name)
            print(f"Copied: {file_path.name}")
            

def load_cluster_map(cluster_file: Path):
    """
    Load a cluster mapping file into a dictionary.
    """
    cluster_map = {}

    with open(cluster_file, "r") as f:
        for line in f:
            if line.strip():
                transcript_id, cluster_id = line.strip().split()
                cluster_map[transcript_id] = cluster_id

    return cluster_map


def get_sample_prefix(fasta_name: str):
    """
    Extract sample prefix from FASTA filename.

    Example:
    DE115_cephalantha_f.largest_cluster_transcripts.fa.transdecoder.pep
    → DE115_cephalantha_f
    """
    return fasta_name.split(".largest_cluster_transcripts")[0]


def main():
    args = get_options()

    # ---------------------------------------------------------------
    # STEP 1: Copy FASTA files to output directory
    # ---------------------------------------------------------------
    copy_input_files(args.input_data, args.output_dir)

    # ---------------------------------------------------------------
    # STEP 2: Process each copied FASTA file
    # ---------------------------------------------------------------
    for fasta_path in args.output_dir.iterdir():
        if not fasta_path.is_file():
            continue

        # Extract sample prefix from FASTA filename
        sample_prefix = get_sample_prefix(fasta_path.name)

        # Construct expected cluster file path:
        # cluster_root / SAMPLE / SAMPLE_salmon-clusters.txt
        cluster_file = (
            args.cluster_root
            / sample_prefix
            / f"{sample_prefix}_salmon-clusters.txt"
        )

        if not cluster_file.exists():
            print(f"ERROR: Cluster file not found for {fasta_path.name}")
            print(f"Expected: {cluster_file}")
            continue

        print(f"Using cluster file: {cluster_file}")

        # Load cluster mappings for this sample
        cluster_map = load_cluster_map(cluster_file)

        records = []

        # Read and modify FASTA records
        for record in SeqIO.parse(fasta_path, "fasta"):

            # Transcript ID may include isoform suffix (.pN)
            transcript_id = record.id
            base_transcript_id = transcript_id.split('.')[0]

            cluster_id = cluster_map.get(base_transcript_id)
            if cluster_id is None:
                print(f"WARNING: No cluster found for {base_transcript_id}")
                continue

            # Preserve entire original FASTA header
            original_header = record.description

            # New header format:
            # sample|cluster__original_header
            new_header = f"{sample_prefix}|{cluster_id}__{original_header}"

            new_record = SeqRecord(
                seq=record.seq,
                id=new_header,
                description=""
            )

            records.append(new_record)

        # Overwrite copied FASTA file
        SeqIO.write(records, fasta_path, "fasta")
        print(f"Processed: {fasta_path.name}")

    print("Done!")


if __name__ == "__main__":
    main()
