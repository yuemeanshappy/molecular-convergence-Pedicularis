#!/usr/bin/env python3

"""
Map RNA-seq reads to a Trinity transcriptome and cluster transcripts with Corset.

Recommended workflow:
    python corset_wrapper.py transcript.fasta samples.tsv threads output_dir salmon-tsv [salmon_index_dir]

TSV format (four tab-separated columns):
    species    sample_name    R1_fastq    R2_fastq

Each TSV row is quantified independently with Salmon.
All resulting equivalence-class files are passed to one Corset run.
Tissue groups are inferred from sample_name.
"""

import os
import re
import shutil
import subprocess
import sys


# ============================================================
# GENERAL HELPERS
# ============================================================

def get_transcript_name(transcript):
    """Return the FASTA filename without its extension."""
    return os.path.basename(transcript).split(".", 1)[0]


def prepare_directory(directory):
    """Create a directory and return its absolute path."""
    directory = os.path.abspath(directory)
    os.makedirs(directory, exist_ok=True)
    return directory


def run_command(command, description=None, stdout=None):
    """Print and execute a command, raising an error if it fails."""
    if description:
        print("\n" + description)

    print(" ".join(str(x) for x in command))
    subprocess.run(command, check=True, stdout=stdout)


# ============================================================
# SALMON INDEX
# ============================================================

def salmon_index(transcript, num_cores, output_dir):
    """
    Build/reuse a Salmon index in output_dir.

    Index directory:
        <transcript_basename>_salmon_index
    """
    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    index_dir = os.path.join(
        output_dir, transcript_name + "_salmon_index"
    )

    version_file = os.path.join(index_dir, "versionInfo.json")

    if os.path.isfile(version_file):
        print("Found existing Salmon index:", index_dir)
        return index_dir

    if os.path.exists(index_dir):
        print("Removing incomplete Salmon index:", index_dir)
        shutil.rmtree(index_dir)

    command = [
        "salmon", "index",
        "-t", transcript,
        "-i", index_dir,
        "--type", "quasi",
        "-p", str(num_cores)
    ]

    run_command(command, "Building Salmon index")

    if not os.path.isfile(version_file):
        raise RuntimeError(
            "Salmon index was not created successfully: " + index_dir
        )

    return index_dir


def get_or_build_salmon_index(
    transcript, num_cores, output_dir, salmon_index_dir=None
):
    """
    Reuse the exact index directory supplied by the user.
    If omitted, build/reuse the index inside output_dir.
    """

    if salmon_index_dir:
        index_dir = os.path.abspath(salmon_index_dir)

        # Check for a valid Salmon index
        version_file = os.path.join(index_dir, "versionInfo.json")
        info_file = os.path.join(index_dir, "info.json")

        if not (
            os.path.isfile(version_file)
            or os.path.isfile(info_file)
        ):
            raise FileNotFoundError(
                "A valid Salmon index was not found in: " + index_dir
            )

        print("Reusing supplied Salmon index:")
        print(index_dir)
        return index_dir

    print("No external Salmon index supplied.")
    return salmon_index(transcript, num_cores, output_dir)


# ============================================================
# SALMON QUANTIFICATION
# ============================================================

def _salmon_quantify(
    index, read1_files, read2_files, output_dir, num_cores,
    single_end=False
):
    """Internal helper for Salmon quantification."""

    eq_classes = os.path.join(
        output_dir, "aux_info", "eq_classes.txt"
    )

    if os.path.isfile(eq_classes) and os.path.getsize(eq_classes) > 0:
        print("Found existing Salmon quantification:", output_dir)
        return eq_classes

    if os.path.exists(output_dir):
        print("Removing incomplete Salmon quantification:", output_dir)
        shutil.rmtree(output_dir)

    command = [
        "salmon", "quant",
        "-i", index,
        "--dumpEq",
        "--libType", "A",
        "-p", str(num_cores)
    ]

    if single_end:
        command += ["-r"] + read1_files
    else:
        if len(read1_files) != len(read2_files):
            raise ValueError(
                "The numbers of R1 and R2 FASTQ files do not match."
            )

        command += ["-1"] + read1_files
        command += ["-2"] + read2_files

    command += ["-o", output_dir]

    run_command(command, "Running Salmon quantification")

    if not os.path.isfile(eq_classes):
        raise RuntimeError(
            "Salmon equivalence classes were not created: " + eq_classes
        )

    return eq_classes


def salmon_quantify_pe(
    transcript, index, pe_fq1, pe_fq2, num_cores, output_dir
):
    """Original paired-end, single-library Salmon workflow."""
    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    quant_dir = os.path.join(
        output_dir, transcript_name + "_salmon_quant"
    )

    return _salmon_quantify(
        index=index,
        read1_files=pe_fq1.split(","),
        read2_files=pe_fq2.split(","),
        output_dir=quant_dir,
        num_cores=num_cores,
        single_end=False
    )


def salmon_quantify_se(
    transcript, index, se_fq, num_cores, output_dir
):
    """Original single-end Salmon workflow."""
    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    quant_dir = os.path.join(
        output_dir, transcript_name + "_salmon_quant"
    )

    return _salmon_quantify(
        index=index,
        read1_files=se_fq.split(","),
        read2_files=[],
        output_dir=quant_dir,
        num_cores=num_cores,
        single_end=True
    )


def salmon_quantify_pe_tissue(
    transcript, index, pe_fq1, pe_fq2,
    num_cores, output_dir, sample_name
):
    """Quantify one paired-end sample into its own output directory."""
    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    quant_dir = os.path.join(
        output_dir,
        transcript_name + "_salmon_quant_" + sample_name
    )

    return _salmon_quantify(
        index=index,
        read1_files=pe_fq1.split(","),
        read2_files=pe_fq2.split(","),
        output_dir=quant_dir,
        num_cores=num_cores,
        single_end=False
    )


# ============================================================
# INFER TISSUE FROM SAMPLE NAME
# ============================================================

def infer_tissue(sample_name):
    """
    Infer flower/leaf tissue from exact sample-name tokens.

    Supported examples:
        DE127_T_flowers -> flower
        DE127_T_leaves  -> leaf
        DE131_F_FT_I1   -> flower
        DE131_F_LT_I1   -> leaf

    Tokens are split on underscores, hyphens, and periods.

    The single-letter token F is deliberately NOT used because it
    could mean something else in a sample identifier.

    FT and LT are recognized only as complete tokens.
    """

    tokens = [
        token.lower()
        for token in re.split(r"[_\-.]+", sample_name)
        if token
    ]

    flower_tokens = {"flower", "flowers", "ft"}
    leaf_tokens = {"leaf", "leaves", "lt"}

    has_flower = any(
        token in flower_tokens for token in tokens
    )

    has_leaf = any(
        token in leaf_tokens for token in tokens
    )

    if has_flower and has_leaf:
        raise ValueError(
            "Ambiguous tissue in sample name {!r}: it contains both "
            "flower and leaf markers. Please rename the sample or "
            "extend the infer_tissue() rules.".format(sample_name)
        )

    if has_flower:
        return "flower"

    if has_leaf:
        return "leaf"

    raise ValueError(
        "Could not infer tissue from sample name {!r}. Expected a "
        "complete token such as 'flowers', 'leaves', 'FT', or 'LT'."
        .format(sample_name)
    )


# ============================================================
# READ SAMPLE TSV
# ============================================================

def read_samples_tsv(samples_tsv):
    """
    Read a four-column TSV:

        species    sample_name    R1_fastq    R2_fastq

    A header is optional.
    Blank lines and lines beginning with # are ignored.
    Relative FASTQ paths are resolved relative to the TSV directory.
    """

    samples_tsv = os.path.abspath(samples_tsv)

    if not os.path.isfile(samples_tsv):
        raise FileNotFoundError(
            "Samples TSV file not found: " + samples_tsv
        )

    tsv_dir = os.path.dirname(samples_tsv)
    samples = []

    with open(samples_tsv, "r", encoding="utf-8") as handle:

        for line_number, line in enumerate(handle, start=1):

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) != 4:
                raise ValueError(
                    "Expected exactly 4 tab-separated columns at "
                    "line {} in {}".format(line_number, samples_tsv)
                )

            species, sample_name, r1, r2 = [
                field.strip() for field in fields
            ]

            # Detect an optional header.
            if (
                sample_name.lower() in ("sample", "sample_name")
                and r1.lower() in ("r1", "r1_fastq", "read1")
            ):
                continue

            if not species or not sample_name or not r1 or not r2:
                raise ValueError(
                    "Missing species, sample name, or FASTQ path "
                    "at line {}".format(line_number)
                )

            # Resolve relative FASTQ paths relative to the TSV.
            if not os.path.isabs(r1):
                r1 = os.path.join(tsv_dir, r1)

            if not os.path.isabs(r2):
                r2 = os.path.join(tsv_dir, r2)

            r1 = os.path.abspath(r1)
            r2 = os.path.abspath(r2)

            if not os.path.isfile(r1):
                raise FileNotFoundError("R1 FASTQ not found: " + r1)

            if not os.path.isfile(r2):
                raise FileNotFoundError("R2 FASTQ not found: " + r2)

            # Filesystem-safe sample name for output directories.
            safe_name = "".join(
                char if char.isalnum() or char in "._-" else "_"
                for char in sample_name
            )

            # Corset uses comma-separated sample names.
            if "," in sample_name:
                raise ValueError(
                    "Sample names cannot contain commas because "
                    "Corset's -n argument is comma-separated: "
                    + sample_name
                )

            tissue = infer_tissue(sample_name)

            samples.append({
                "species": species,
                "sample": safe_name,
                "sample_label": sample_name,
                "tissue": tissue,
                "r1": r1,
                "r2": r2
            })

    if not samples:
        raise ValueError(
            "No samples were found in " + samples_tsv
        )

    names = [sample["sample"] for sample in samples]

    if len(names) != len(set(names)):
        raise ValueError(
            "Sample names must be unique within the TSV."
        )

    print("\nSamples read from TSV:")

    for sample in samples:
        print(
            "{} (tissue={}): R1={}, R2={}".format(
                sample["sample_label"],
                sample["tissue"],
                sample["r1"],
                sample["r2"]
            )
        )

    return samples


# ============================================================
# CORSET: SALMON
# ============================================================

def corset_salmon_eq_classes(transcript, eq_classes, output_dir):
    """Original Corset workflow with one Salmon equivalence-class file."""

    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    prefix = os.path.join(
        output_dir, transcript_name + "_salmon"
    )

    clusters_file = prefix + "-clusters.txt"
    counts_file = prefix + "-counts.txt"

    if os.path.isfile(clusters_file) and os.path.isfile(counts_file):
        print(
            "Found existing Corset Salmon files for",
            transcript_name
        )
        return

    command = [
        "corset",
        "-i", "salmon_eq_classes",
        "-m", "5",
        "-p", prefix,
        eq_classes
    ]

    run_command(command, "Running Corset")

    if not (
        os.path.isfile(clusters_file)
        and os.path.isfile(counts_file)
    ):
        raise RuntimeError(
            "Corset did not finish successfully."
        )


def corset_salmon_multiple(
    transcript,
    eq_classes_files,
    output_dir,
    sample_names=None,
    groups=None
):
    """
    Run Corset once with one separate equivalence-class file per sample.

    sample_names and groups must follow exactly the same order as
    eq_classes_files.

    Groups are tissue labels such as flower and leaf.
    """

    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    prefix = os.path.join(
        output_dir, transcript_name + "_salmon"
    )

    clusters_file = prefix + "-clusters.txt"
    counts_file = prefix + "-counts.txt"

    if not eq_classes_files:
        raise ValueError(
            "No Salmon equivalence-class files were supplied."
        )

    for eq_file in eq_classes_files:
        if not os.path.isfile(eq_file):
            raise FileNotFoundError(
                "Salmon equivalence-class file not found: " + eq_file
            )

    if (
        sample_names is not None
        and len(sample_names) != len(eq_classes_files)
    ):
        raise ValueError(
            "Number of Corset sample names does not match input files."
        )

    if (
        groups is not None
        and len(groups) != len(eq_classes_files)
    ):
        raise ValueError(
            "Number of Corset groups does not match input files."
        )

    command = [
        "corset",
        "-i", "salmon_eq_classes",
        "-m", "5"
    ]

    if groups is not None:
        command += ["-g", ",".join(groups)]

    if sample_names is not None:

        if any("," in name for name in sample_names):
            raise ValueError(
                "Corset sample names supplied with -n "
                "cannot contain commas."
            )

        command += ["-n", ",".join(sample_names)]

    command += [
        "-p", prefix,
        "-f", "true"
    ] + eq_classes_files

    run_command(
        command,
        "Running Corset jointly across all samples"
    )

    if not (
        os.path.isfile(clusters_file)
        and os.path.isfile(counts_file)
    ):
        raise RuntimeError(
            "Corset did not finish successfully."
        )


# ============================================================
# SALMON WORKFLOWS
# ============================================================

def run_pe_salmon(
    transcript, pe_fq1, pe_fq2, num_cores, output_dir
):
    """Original paired-end, single-library Salmon workflow."""

    index = salmon_index(
        transcript, num_cores, output_dir
    )

    eq_classes = salmon_quantify_pe(
        transcript,
        index,
        pe_fq1,
        pe_fq2,
        num_cores,
        output_dir
    )

    corset_salmon_eq_classes(
        transcript, eq_classes, output_dir
    )


def run_se_salmon(
    transcript, se_fq, num_cores, output_dir
):
    """Original single-end Salmon workflow."""

    index = salmon_index(
        transcript, num_cores, output_dir
    )

    eq_classes = salmon_quantify_se(
        transcript,
        index,
        se_fq,
        num_cores,
        output_dir
    )

    corset_salmon_eq_classes(
        transcript, eq_classes, output_dir
    )


def run_pe_salmon_tissues(
    transcript,
    flower_fq1,
    flower_fq2,
    leaf_fq1,
    leaf_fq2,
    num_cores,
    output_dir,
    salmon_index_dir=None
):
    """Two-tissue paired-end workflow, retained for compatibility."""

    output_dir = prepare_directory(output_dir)

    index = get_or_build_salmon_index(
        transcript,
        num_cores,
        output_dir,
        salmon_index_dir
    )

    flower_eq = salmon_quantify_pe_tissue(
        transcript,
        index,
        flower_fq1,
        flower_fq2,
        num_cores,
        output_dir,
        "flowers"
    )

    leaf_eq = salmon_quantify_pe_tissue(
        transcript,
        index,
        leaf_fq1,
        leaf_fq2,
        num_cores,
        output_dir,
        "leaves"
    )

    corset_salmon_multiple(
        transcript,
        [flower_eq, leaf_eq],
        output_dir,
        sample_names=["flowers", "leaves"],
        groups=["flower", "leaf"]
    )


def run_salmon_from_tsv(
    transcript,
    samples_tsv,
    num_cores,
    output_dir,
    salmon_index_dir=None
):
    """
    General paired-end workflow for any number of samples in a TSV.

    Each TSV row is quantified separately.
    Corset is run once using all resulting equivalence-class files.
    Tissue groups are inferred from sample names.
    """

    output_dir = prepare_directory(output_dir)

    print("\n========================================")
    print("Starting Salmon/Corset TSV workflow")
    print("========================================")

    samples = read_samples_tsv(samples_tsv)

    transcript_name = get_transcript_name(transcript)

    # All rows in this TSV must belong to the same species.
    species_labels = sorted(
        set(sample["species"] for sample in samples)
    )

    if len(species_labels) > 1:
        raise ValueError(
            "The TSV contains more than one species label: "
            + ", ".join(species_labels)
        )

    index = get_or_build_salmon_index(
        transcript,
        num_cores,
        output_dir,
        salmon_index_dir
    )

    eq_classes_files = []

    # Run Salmon independently for every TSV row.
    for sample in samples:

        print("\n----------------------------------------")
        print("Processing sample:", sample["sample_label"])
        print("Tissue:", sample["tissue"])
        print("----------------------------------------")

        eq_file = salmon_quantify_pe_tissue(
            transcript=transcript,
            index=index,
            pe_fq1=sample["r1"],
            pe_fq2=sample["r2"],
            num_cores=num_cores,
            output_dir=output_dir,
            sample_name=sample["sample"]
        )

        eq_classes_files.append(eq_file)

    # Preserve the original sample names for Corset's count-file header.
    corset_sample_names = [
        sample["sample_label"] for sample in samples
    ]

    # Assign the same Corset group to samples from the same tissue.
    corset_groups = [
        sample["tissue"] for sample in samples
    ]

    print("\nCorset input order:")

    for i, sample in enumerate(samples):
        print(
            "{}: {} ({})".format(
                i + 1,
                sample["sample_label"],
                sample["tissue"]
            )
        )

    print(
        "\nCorset sample names:",
        ",".join(corset_sample_names)
    )

    print(
        "Corset tissue groups:",
        ",".join(corset_groups)
    )

    # One Corset run, all independent Salmon equivalence-class files.
    corset_salmon_multiple(
        transcript,
        eq_classes_files,
        output_dir,
        sample_names=corset_sample_names,
        groups=corset_groups
    )

    print("\n========================================")
    print("Salmon/Corset TSV workflow finished")
    print("Species:", species_labels[0])
    print("Number of samples:", len(samples))
    print("========================================")


# ============================================================
# BOWTIE INDEX
# ============================================================

def bowtie_index(transcript, num_cores, output_dir):

    output_dir = prepare_directory(output_dir)
    transcript = os.path.abspath(transcript)

    transcript_name = get_transcript_name(transcript)

    index_prefix = os.path.join(
        output_dir, transcript_name
    )

    index_files = [
        index_prefix + ".1.ebwt",
        index_prefix + ".2.ebwt",
        index_prefix + ".3.ebwt",
        index_prefix + ".4.ebwt",
        index_prefix + ".rev.1.ebwt",
        index_prefix + ".rev.2.ebwt"
    ]

    if all(os.path.isfile(path) for path in index_files):
        print("Found Bowtie index for", transcript_name)
        return index_prefix

    command = [
        "bowtie-build",
        "-q",
        transcript,
        index_prefix,
        "--threads",
        str(num_cores)
    ]

    run_command(command, "Building Bowtie index")

    if not all(os.path.isfile(path) for path in index_files):
        raise RuntimeError(
            "Bowtie index was not completed."
        )

    return index_prefix


# ============================================================
# BOWTIE ALIGNMENT
# ============================================================

def bowtie_align_pe(
    transcript, index, pe_fq1, pe_fq2,
    num_cores, output_dir
):

    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    sam_file = os.path.join(
        output_dir, transcript_name + "_bowtie.sam"
    )

    bam_file = os.path.join(
        output_dir, transcript_name + "_bowtie.bam"
    )

    if os.path.isfile(bam_file):
        print(
            "Found Bowtie BAM alignment for",
            transcript_name
        )
        return bam_file

    if not os.path.isfile(sam_file):

        command = [
            "bowtie",
            "--all",
            "-k", "10",
            "--chunkmbs", "2580",
            "--rf",
            "--nofw",
            "--threads", str(num_cores),
            index,
            "-1", pe_fq1,
            "-2", pe_fq2,
            "-S", sam_file
        ]

        run_command(
            command,
            "Running Bowtie paired-end alignment"
        )

    if not os.path.isfile(sam_file):
        raise RuntimeError(
            "Bowtie alignment was not completed."
        )

    with open(bam_file, "wb") as bam_out:

        command = [
            "samtools", "view",
            "-S", "-b", sam_file
        ]

        run_command(
            command,
            "Converting SAM to BAM",
            stdout=bam_out
        )

    if not os.path.isfile(bam_file):
        raise RuntimeError(
            "Samtools conversion was not completed."
        )

    os.remove(sam_file)

    return bam_file


def bowtie_align_se(
    transcript, index, se_fq, num_cores, output_dir
):

    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    sam_file = os.path.join(
        output_dir, transcript_name + "_bowtie.sam"
    )

    bam_file = os.path.join(
        output_dir, transcript_name + "_bowtie.bam"
    )

    if os.path.isfile(bam_file):
        print(
            "Found Bowtie BAM alignment for",
            transcript_name
        )
        return bam_file

    if not os.path.isfile(sam_file):

        command = [
            "bowtie",
            "--all",
            "-k", "10",
            "--threads", str(num_cores),
            index,
            se_fq,
            "-S", sam_file
        ]

        run_command(
            command,
            "Running Bowtie single-end alignment"
        )

    if not os.path.isfile(sam_file):
        raise RuntimeError(
            "Bowtie alignment was not completed."
        )

    with open(bam_file, "wb") as bam_out:

        command = [
            "samtools", "view",
            "-S", "-b", sam_file
        ]

        run_command(
            command,
            "Converting SAM to BAM",
            stdout=bam_out
        )

    if not os.path.isfile(bam_file):
        raise RuntimeError(
            "Samtools conversion was not completed."
        )

    os.remove(sam_file)

    return bam_file


# ============================================================
# CORSET: BOWTIE
# ============================================================

def corset_bowtie(transcript, bam, output_dir):

    output_dir = prepare_directory(output_dir)
    transcript_name = get_transcript_name(transcript)

    prefix = os.path.join(
        output_dir, transcript_name + "_bowtie"
    )

    clusters_file = prefix + "-clusters.txt"
    counts_file = prefix + "-counts.txt"

    if os.path.isfile(clusters_file) and os.path.isfile(counts_file):
        print(
            "Found existing Corset Bowtie files for",
            transcript_name
        )
        return

    command = [
        "corset",
        "-i", "bam",
        "-m", "5",
        "-p", prefix,
        bam
    ]

    run_command(
        command,
        "Running Corset with Bowtie BAM"
    )

    if not (
        os.path.isfile(clusters_file)
        and os.path.isfile(counts_file)
    ):
        raise RuntimeError(
            "Corset did not finish successfully."
        )


# ============================================================
# BOWTIE WORKFLOWS
# ============================================================

def run_pe_bowtie(
    transcript, pe_fq1, pe_fq2,
    num_cores, output_dir
):

    output_dir = prepare_directory(output_dir)

    index = bowtie_index(
        transcript, num_cores, output_dir
    )

    bam_file = bowtie_align_pe(
        transcript,
        index,
        pe_fq1,
        pe_fq2,
        num_cores,
        output_dir
    )

    corset_bowtie(
        transcript, bam_file, output_dir
    )


def run_se_bowtie(
    transcript, se_fq, num_cores, output_dir
):

    output_dir = prepare_directory(output_dir)

    index = bowtie_index(
        transcript, num_cores, output_dir
    )

    bam_file = bowtie_align_se(
        transcript,
        index,
        se_fq,
        num_cores,
        output_dir
    )

    corset_bowtie(
        transcript, bam_file, output_dir
    )


# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================

def print_usage():

    print("""
Usage:

1. Recommended: paired-end Salmon workflow using a TSV

   Reuse an existing index:
   python corset_wrapper.py transcript.fasta samples.tsv threads output_dir salmon-tsv [index_dir]

   Build/reuse index in output_dir:
   python corset_wrapper.py transcript.fasta samples.tsv threads output_dir salmon-tsv

2. Two-tissue paired-end Salmon workflow:
   python corset_wrapper.py transcript.fasta flower_R1 flower_R2 leaf_R1 leaf_R2 threads output_dir salmon-pe-tissues [index_dir]

3. Original single-end Salmon:
   python corset_wrapper.py transcript.fasta reads.fastq threads output_dir salmon

4. Original paired-end Salmon:
   python corset_wrapper.py transcript.fasta R1.fastq R2.fastq threads output_dir salmon

5. Original single-end Bowtie:
   python corset_wrapper.py transcript.fasta reads.fastq threads output_dir bowtie

6. Original paired-end Bowtie:
   python corset_wrapper.py transcript.fasta R1.fastq R2.fastq threads output_dir bowtie
""")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    try:

        # New general TSV workflow:
        # transcript samples.tsv threads output_dir salmon-tsv [index_dir]

        if (
            len(sys.argv) in (6, 7)
            and sys.argv[5] == "salmon-tsv"
        ):

            existing_index = (
                sys.argv[6] if len(sys.argv) == 7 else None
            )

            run_salmon_from_tsv(
                transcript=sys.argv[1],
                samples_tsv=sys.argv[2],
                num_cores=int(sys.argv[3]),
                output_dir=sys.argv[4],
                salmon_index_dir=existing_index
            )

        # Two-tissue paired-end workflow:

        elif (
            len(sys.argv) in (9, 10)
            and sys.argv[8] == "salmon-pe-tissues"
        ):

            existing_index = (
                sys.argv[9] if len(sys.argv) == 10 else None
            )

            run_pe_salmon_tissues(
                transcript=sys.argv[1],
                flower_fq1=sys.argv[2],
                flower_fq2=sys.argv[3],
                leaf_fq1=sys.argv[4],
                leaf_fq2=sys.argv[5],
                num_cores=int(sys.argv[6]),
                output_dir=sys.argv[7],
                salmon_index_dir=existing_index
            )

        # Original single-end Salmon:

        elif (
            len(sys.argv) == 6
            and sys.argv[5] == "salmon"
        ):

            run_se_salmon(
                transcript=sys.argv[1],
                se_fq=sys.argv[2],
                num_cores=int(sys.argv[3]),
                output_dir=sys.argv[4]
            )

        # Original paired-end Salmon:

        elif (
            len(sys.argv) == 7
            and sys.argv[6] == "salmon"
        ):

            run_pe_salmon(
                transcript=sys.argv[1],
                pe_fq1=sys.argv[2],
                pe_fq2=sys.argv[3],
                num_cores=int(sys.argv[4]),
                output_dir=sys.argv[5]
            )

        # Original single-end Bowtie:

        elif (
            len(sys.argv) == 6
            and sys.argv[5] == "bowtie"
        ):

            run_se_bowtie(
                transcript=sys.argv[1],
                se_fq=sys.argv[2],
                num_cores=int(sys.argv[3]),
                output_dir=sys.argv[4]
            )

        # Original paired-end Bowtie:

        elif (
            len(sys.argv) == 7
            and sys.argv[6] == "bowtie"
        ):

            run_pe_bowtie(
                transcript=sys.argv[1],
                pe_fq1=sys.argv[2],
                pe_fq2=sys.argv[3],
                num_cores=int(sys.argv[4]),
                output_dir=sys.argv[5]
            )

        else:
            print_usage()
            sys.exit(1)

    except (
        ValueError,
        FileNotFoundError,
        RuntimeError,
        subprocess.CalledProcessError
    ) as error:

        print("\nERROR:", error, file=sys.stderr)
        sys.exit(1)
