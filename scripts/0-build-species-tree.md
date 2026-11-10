# Species Tree Construction

## Objective

Construct a species-level phylogenetic tree of *Pedicularis* using 83 RNA-seq datasets (leaf and flower tissues combined) representing 51 species.

## Workflow overview

1. Download publicly available RNA-seq data from NCBI. `NCBI SRA Toolkit v2.11.3`
2. Trim adapters and remove low-quality sequences.`Trimmomatic v0.39`
3. Perform de novo transcriptome assembly.`Trinity v2.14.0`
4. Identify organelle sequences and separate nuclear transcripts.`GetOrganelle v1.7.7.1`; `fasta_force_to_gfa.py`;`separate_trinity_organelle.py`
5. Remove chimeric transcripts, cluster transcripts, and select the longest transcript from each Corset cluster.`Salmon v0.9.1`; `Corset v1.0.9`;`run_chimera_detection.py`;`detect_chimera_from_blastx_modifed.py`;`corset_wrapper.py`;`filter_corset_output.py` (these scripts are from [link](https://bitbucket.org/yanglab/phylogenomic_dataset_construction/src/master/))
6. Predict coding sequences and protein sequences, and rename sequence headers.`TransDecoder v5.7.1`;`update_rename_header.py`
7. Identify orthogroups among species.`OrthoFinder v3.1.2`
8. Perform orthogroup sequence selection, alignment, and post-processing.`twig v.0.0.5`
9. Trim alignments, infer gene trees, and construct a draft species tree. `twig v.0.0.5`;`RAxML-NG v2.0.0`;`ASTRAL-Pro3 v1.22.3.6`
10. Root individual gene trees.`twig v.0.0.5`
11. Infer final species trees using different gene tree datasets.`twig v.0.0.5`;`ASTRAL-Pro3 v1.22.3.6`

## Reference data and data availability

Selected intermediate and final phylogenomic results are available in the associated GitHub repository, including:

- Refined orthogroup sequence alignments: `./data/1.orthogroup_alignments`
- Filtered and rooted gene trees: `./data/2.orthogroup_gene_trees`
- Final species trees: `./data/3.species_tree`

---

## Step 1. Download RNA-seq data

1. **Input:** NCBI SRA accession numbers and sample metadata.
2. **Output:** Raw paired-end RNA-seq reads (`raw_fastq_dir`).
3. **Software:** NCBI SRA Toolkit v.2.11.3 

**Example code:**

```bash
# Download RNA-seq data for each SRA accession
while read accession; do

    prefetch ${accession}

    fasterq-dump ${accession} \
        --split-files \
        --threads 8 \
        --outdir ${raw_fastq_dir}

done < ${metadata_file}
```

## Step 2. Trim adapters and remove low-quality sequences

1. **Input:** Raw paired-end RNA-seq reads (`raw_fastq_dir`).
2. **Output:** Trimmed paired-end reads (`trimmed_fastq_dir`).
3. **Software:** Trimmomatic v0.39.

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Remove adapters and low-quality bases
    trimmomatic PE \
        -threads 15 \
        ${raw_fastq_dir}/${sample}_1.fq.gz \
        ${raw_fastq_dir}/${sample}_2.fq.gz \
        ${trimmed_fastq_dir}/${sample}_1_paired.fq.gz \
        ${trimmed_fastq_dir}/${sample}_1_unpaired.fq.gz \
        ${trimmed_fastq_dir}/${sample}_2_paired.fq.gz \
        ${trimmed_fastq_dir}/${sample}_2_unpaired.fq.gz \
        ILLUMINACLIP:${adapter_file}:2:30:10:2:True \
        SLIDINGWINDOW:4:5 \
        LEADING:5 \
        TRAILING:5 \
        MINLEN:25

done
```

## Step 3. Perform de novo transcriptome assembly

1. **Input:** Trimmed paired-end RNA-seq reads (`trimmed_fastq_dir`).
2. **Output:** Assembled transcript sequences (`transcriptome_assembly_dir`).
3. **Software:** Trinity v2.14.0.

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Assemble transcriptomes independently for each sample
    # Input: Paired-end trimmed FASTQ files
    # Output: Trinity transcriptome assembly
    # Software: Trinity v2.14.0

    Trinity \
        --seqType fq \
        --max_memory 80G \
        --CPU 10 \
        --left "${trimmed_fastq_dir}/${sample}_1_paired.fq.gz" \
        --right "${trimmed_fastq_dir}/${sample}_2_paired.fq.gz" \
        --output "${transcriptome_assembly_dir}/${sample}.trinity_out"

    # Remove unnecessary intermediate files only after successful assembly
    if [ $? -eq 0 ]; then

        assembly_dir="${transcriptome_assembly_dir}/${sample}.trinity_out"

        rm -f "${assembly_dir}"/*fq.gz*
        rm -f "${assembly_dir}"/jellyfish.kmers.*.fa
        rm -f "${assembly_dir}"/both.fa
        rm -f "${assembly_dir}"/inchworm.*.fa

        # Compress scaffolding alignment file if it exists
        if [ -f "${assembly_dir}/scaffolding_entries.sam" ]; then
            gzip "${assembly_dir}/scaffolding_entries.sam"
        fi

        # Remove large intermediate directories
        rm -rf "${assembly_dir}/insilico_read_normalization"
        rm -rf "${assembly_dir}/read_partitions"
        rm -rf "${assembly_dir}/chrysalis"

    else
        echo "Assembly failed for ${sample}"
    fi

done
```

## Step 4. Separate nuclear, plastid, and mitochondrial transcripts

1. **Input:** Assembled transcript sequences (`transcriptome_assembly_dir`).
2. **Output:**
   - Nuclear transcripts (`nuclear_fasta_dir`).
   - Plastid transcripts (`plastid_fasta_dir`).
   - Mitochondrial transcripts (`mitochondrial_fasta_dir`).
3. **Software:** GetOrganelle v1.7.7.1; `fasta_force_to_gfa.py`;`separate_trinity_organelle.py`

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Convert Trinity assembly FASTA to GFA format
    python fasta_force_to_gfa.py \
        "${transcriptome_assembly_dir}/${sample}.fasta"

    # Identify plastid and mitochondrial sequences
    python get_organelle_from_assembly.py \
        -F embplant_pt,embplant_mt \
        --slim-options "--perc-identity 90.0 --no-hits ex_no_hit" \
        -t 4 \
        -g "${getorganelle_output_dir}/${sample}.gfa" \
        -o "${getorganelle_output_dir}/${sample}.getorganelle"

    # Separate Trinity transcripts into nuclear and organelle sequences
    python separate_trinity_organelle.py single \
        --sample "${sample}" \
        --trinity "${transcriptome_assembly_dir}/${sample}.fasta" \
        --name-map "${getorganelle_output_dir}/${sample}.gfa.name_map" \
        --pt "${getorganelle_output_dir}/${sample}.getorganelle/embplant_pt.scaffolds.graph1.1.path_sequence.fasta" \
        --mt "${getorganelle_output_dir}/${sample}.getorganelle/embplant_mt.scaffolds.graph1.1.path_sequence.fasta" \
        --out-dir "${separated_organelle_output_dir}/${sample}"

done
```

## Step 5. Remove chimeric transcripts, cluster transcripts, and select the longest transcript

### 5.1. Remove chimeric transcripts

1. **Input:** Nuclear transcript sequences (`separated_organelle_output_dir`); reference protein sequences from *P. cranolopha*.
2. **Output:** Chimera-filtered nuclear transcript sequences (`filtered_nuclear_fasta_dir`).
3. **Software:** `run_chimera_detection.py`; `detect_chimera_from_blastx_modifed.py`

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Identify and remove potentially chimeric transcripts
    python run_chimera_detection.py \
        ${separated_organelle_output_dir}/${sample}.nucl.fasta \
        ${reference_protein_file} \
        50 \
        ${filtered_nuclear_fasta_dir}/${sample}

done
```

### 5.2. Cluster transcripts using read-mapping information

1. **Input:** Chimera-filtered nuclear transcripts (`filtered_nuclear_fasta_dir`) and corresponding paired-end RNA-seq reads (`trimmed_fastq_dir`).
2. **Output:** Transcript-to-cluster assignment files (`transcript_cluster_dir`).
3. **Software:** Salmon v0.9.1; Corset v1.0.9; `corset_wrapper.py`.

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Step 5.2: Quantify transcript abundance using Salmon and cluster
    # transcripts using Corset
    #
    # Difference between Trinity and Corset clustering:
    #
    # Trinity groups assembled transcripts based on sequence assembly
    # relationships. Its transcript IDs follow the format c1_g1_i1:
    #   c1 = Trinity component/cluster ID
    #   g1 = Trinity gene ID
    #   i1 = Trinity isoform ID
    #
    # Corset re-evaluates transcript relationships by mapping sequencing
    # reads to transcripts and grouping transcripts according to their
    # shared-read proportions and expression patterns. Therefore, Corset
    # clusters do not necessarily correspond to Trinity clusters or genes.
    # Transcripts assigned to the same Trinity gene may be separated into
    # different Corset clusters if they have insufficient shared reads.
    #
    # Difference between Corset clusters and subclusters:
    #
    #   Cluster 19.0 vs. Cluster 20.0:
    #       These are two different Corset clusters. They have no shared
    #       cluster assignment and are treated as separate groups,
    #       potentially representing different genes.
    #
    #   Cluster 19.0 vs. Cluster 19.1:
    #       These are two subclusters within the same broader Corset
    #       cluster family (19). They have been separated by Corset
    #       because of differences in read-sharing or expression patterns.
    #       They may represent different genes or distinct transcript
    #       groups, such as isoforms or paralogs.
    #
    # Examples from DE115_salmon-clusters.txt:
    #   DN152075_c0_g1_i2 -> Corset cluster 65118.0
    #   DN152075_c0_g1_i1 -> Corset cluster 68524.0
    #       Two Trinity isoforms assigned to different Corset clusters.
    #
    #   DN4432_c0_g1_i2 -> Corset cluster 19.0
    #   DN4432_c0_g1_i9 -> Corset cluster 19.1
    #       Two Trinity isoforms assigned to different Corset subclusters.
    #
    # Output interpretation:
    #   salmon-counts.txt contains the number of reads mapped to each
    #   Corset cluster, rather than counts for individual transcripts.
    #
    #   The Corset cluster assignments are subsequently used in Step 5.3
    #   to select the longest transcript from each cluster/subcluster.
    python corset_wrapper.py \
        ${filtered_nuclear_fasta_dir}/${sample}.filtered_transcripts.fa \
        ${trimmed_fastq_dir}/${sample}_1.fq.gz \
        ${trimmed_fastq_dir}/${sample}_2.fq.gz \
        50 \
        ${transcript_cluster_dir}/${sample} \
        salmon

done
```

### 5.3. Select the longest transcript from each cluster

1. **Input:** Chimera-filtered nuclear transcripts (`filtered_nuclear_fasta_dir`) and transcript-to-cluster assignment files (`transcript_cluster_dir`).
2. **Output:** Representative nuclear transcripts (`representative_transcript_dir`).
3. **Software:** `filter_corset_output.py`.

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Extract the longest transcript from each Corset cluster
    # The script combines the Corset cluster assignments with transcript
    # lengths calculated from the nuclear transcript FASTA file.
    #
    # For each Corset cluster/subcluster, transcripts are sorted by
    # sequence length in descending order, and only the longest
    # transcript is retained as the representative sequence.
    #
    # All other transcripts within the same Corset cluster/subcluster
    # are classified as redundant and written to a separate FASTA file.
    #
    # Cluster labels are treated independently (e.g., 19.0 and 19.1
    # are considered two separate clusters).

    python filter_corset_output.py \
        ${filtered_nuclear_fasta_dir}/${sample}.filtered_transcripts.fa \
        ${transcript_cluster_dir}/${sample}_salmon-clusters.txt \
        ${representative_transcript_dir}/${sample}.largest_cluster_transcripts.fa

done
```

## Step 6. Predict ORFs and obtain CDS and amino acid sequences

### 6.1. Predict ORFs

1. **Input:** Representative nuclear transcript sequences (`representative_transcript_dir`).
2. **Output:** Predicted CDS and amino acid sequences (`cds_protein_dir`).
3. **Software:** TransDecoder v5.7.1.

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Predict complete ORFs and retain the best ORF per transcript
    TransDecoder \
        -t ${representative_transcript_dir}/${sample}.largest_cluster_transcripts.fa \
        --single_best_only \
        --complete-ORF \
        -O ${cds_protein_dir}/${sample}

done
```

### 6.2. Rename sequence identifiers

1. **Input:** Predicted CDS and amino acid sequences (`cds_protein_dir`).
2. **Output:** Renamed CDS and amino acid sequences (`renamed_cds_dir`, `renamed_aa_dir`).
3. **Software:** `update_rename_header.py`

**Example code:**

```bash
for sample in ${samples[@]}; do

    # Add species name and Corset cluster ID to sequence identifiers
    python update_rename_header.py \
        --input_cds ${cds_protein_dir}/${sample}.cds.fa \
        --input_aa ${cds_protein_dir}/${sample}.pep.fa \
        --output_cds ${renamed_cds_dir}/${sample}.cds.fa \
        --output_aa ${renamed_aa_dir}/${sample}.aa.fa \
        --sample ${sample} \
        --cluster_file ${transcript_cluster_dir}/${sample}_salmon-clusters.txt

done
```

## Step 7. Identify orthogroups

1. **Input:** Renamed CDS sequences from all samples (`renamed_cds_dir`).
2. **Output:** Orthogroup sequence files (`orthogroup_fasta_dir`).
3. **Software:** OrthoFinder v3.1.2.

**Example code:**

```bash
# Identify orthogroups across all samples
orthofinder \
    -t 32 \
    -a 2 \
    -d \
    -X \
    -f ${renamed_cds_dir} \
    -o ${orthofinder_output_dir}
```

## Step 8. Filter orthogroups, remove redundant isoforms, perform sequence alignment, and trim alignments

1. **Input:** Orthogroup sequence files (`orthogroup_fasta_dir`).
2. **Output:**
   - Filtered orthogroup sequences (`filtered_orthogroup_dir`).
   - Multiple sequence alignments (`orthogroup_alignment_dir`).
   - Refined and trimmed alignments (`trimmed_alignment_dir`).
3. **Software:** twig v.0.0.5

**Example code:**

```bash
for orthogroup in ${orthogroups[@]}; do

    # Filter orthogroups and remove redundant isoforms
    twig align-pre \
        -i ${orthogroup_fasta_dir}/${orthogroup}.fa \
        -o ${filtered_orthogroup_dir}/${orthogroup}.trim.nt.fa \
        -X ${isoform_rules_file} \
        -d '_' \
        -di 0 \
        -hf 0.3 \
        -hi 0.6 \
        -ml 200 \
        -mc 25 \
        -hc 15 \
        -ti 33 \
        -te 33 \
        -f

    # Perform codon-aware CDS sequence alignment
    twig align-cds \
        -i ${filtered_orthogroup_dir}/${orthogroup}.trim.nt.fa \
        -o ${orthogroup_alignment_dir}/${orthogroup}.msa_raw.nt.fa

    # Trim poorly aligned regions
    twig align-post \
        -i ${orthogroup_alignment_dir}/${orthogroup}.msa_raw.nt.fa \
        -o ${trimmed_alignment_dir}/${orthogroup}.msa.nt.fa \
        -mo 100 \
        -ms 25 \
        -ac 0.5

done
```

## Step 9. Infer ML gene trees, filter gene trees, and construct a draft species tree

1. **Input:** Refined CDS alignments (`trimmed_alignment_dir`).
2. **Output:**
   - ML gene trees with bootstrap support (`gene_tree_dir`).
   - Filtered gene trees (`filtered_gene_tree_dir`).
   - Draft species tree (`draft_species_tree_dir`).
3. **Software:** RAxML-NG v2.0.0; twig v.0.0.5; ASTRAL-Pro3 v1.22.3.6

**Example code:**

```bash
for orthogroup in ${orthogroups[@]}; do

    # Partition CDS alignments by codon position
    twig partition-cds \
        -i ${trimmed_alignment_dir}/${orthogroup}.msa.nt.fa \
        -n 3 \
        -m "GTR+G" \
        > ${partition_dir}/${orthogroup}.partitions

    # Infer partitioned ML gene trees with 100 bootstrap replicates
    raxml-ng \
        --all \
        --msa ${trimmed_alignment_dir}/${orthogroup}.msa.nt.fa \
        --model ${partition_dir}/${orthogroup}.partitions \
        --brlen scaled \
        --seed 42 \
        --bs-trees 100 \
        --prefix ${gene_tree_dir}/${orthogroup}

    # Collapse branches with bootstrap support below 50%
    # Retain gene trees containing at least one designated outgroup
    twig tree-filter \
        -i ${gene_tree_dir}/${orthogroup}.raxml.support \
        -o ${filtered_gene_tree_dir}/${orthogroup}.collapse50.glabel \
        -d "|" \
        -di 0 \
        -s 50 \
        -eo 10 \
        -m 25 \
        -O ${outgroup_mapping_file} \
        --require

done

# Combine filtered gene trees into a single tree set
find ${filtered_gene_tree_dir} \
    -type f -name "*.glabel" -size +1b -print0 \
    | xargs -0 cat \
    > ${draft_species_tree_dir}/tree-set-glabel-filtered.nwk

# Standardize taxon labels
twig tree-filter \
    -i ${draft_species_tree_dir}/tree-set-glabel-filtered.nwk \
    -d "|" \
    -di 0 \
    -rd \
    > ${draft_species_tree_dir}/tree-set-alabel-filtered.nwk

# Infer the draft species tree
astral-pro3 \
    -i ${draft_species_tree_dir}/tree-set-alabel-filtered.nwk \
    -o ${draft_species_tree_dir}/sptree-alabel-filtered.nwk \
    -t 10 \
    --root ${outgroup_accession} \
    --seed 42

# Relabel the draft species tree using species-accession mapping
twig tree-filter \
    -i ${draft_species_tree_dir}/sptree-alabel-filtered.nwk \
    -o ${draft_species_tree_dir}/sptree-alabel-filtered-relabeled.nwk \
    -d "|" \
    -di 0 \
    -ri \
    -I ${species_mapping_file}
```

## Step 10. Root gene trees based on the draft species tree

1. **Input:** Filtered gene trees (`filtered_gene_tree_dir`) and draft species tree (`draft_species_tree_dir`).
2. **Output:** Rooted gene trees (`rooted_gene_tree_dir`).
3. **Software:** twig v.0.0.5

**Example code:**

```bash
for orthogroup in ${filtered_orthogroups[@]}; do

    # First, attempt to root gene trees using designated outgroups
    twig tree-rooter \
        -i ${filtered_gene_tree_dir}/${orthogroup}.collapse50.glabel \
        -d "|" \
        -di 0 \
        -I ${outgroup_mapping_file} \
        -s ${draft_species_tree_dir}/sptree-alabel-filtered.nwk \
        1> ${rooted_gene_tree_dir}/${orthogroup}.collapse50.glabel.rooted \
        2> ${rooted_gene_tree_dir}/${orthogroup}.log

    # If outgroup rooting fails, use P. longicaulis (JJ131)
    if [ ! -s ${rooted_gene_tree_dir}/${orthogroup}.collapse50.glabel.rooted ]; then

        twig tree-rooter \
            -i ${filtered_gene_tree_dir}/${orthogroup}.collapse50.glabel \
            -d "|" \
            -di 0 \
            -r JJ131 \
            -s ${draft_species_tree_dir}/sptree-alabel-filtered.nwk \
            1> ${rooted_gene_tree_dir}/${orthogroup}.collapse50.glabel.rooted \
            2> ${rooted_gene_tree_dir}/${orthogroup}.log

    fi

done
```

## Step 11. Infer the final species tree using ASTRAL-Pro3

1. **Input:** Rooted gene trees (`rooted_gene_tree_dir`) and species-accession mapping information.
2. **Output:**
   - Species tree inferred from multicopy gene trees.
   - Species tree inferred from single-copy gene trees.
   - Species tree files with accession-level and species-level taxon labels.
3. **Software:** twig v.0.0.5 ; ASTRAL-Pro3 v1.22.3.6

**Example code:**

```bash
# Combine all successfully rooted gene trees
find ${rooted_gene_tree_dir} \
    -type f -name "*.glabel.rooted" -size +1b -print0 \
    | xargs -0 cat \
    > ${species_tree_dir}/tree-set-glabel-filtered-rooted.nwk

# Generate multicopy gene trees with species-level labels
twig tree-filter \
    -i ${species_tree_dir}/tree-set-glabel-filtered-rooted.nwk \
    -d "|" \
    -di 0 \
    -I ${species_mapping_file} \
    -ri \
    > ${species_tree_dir}/tree-set-species-filtered-rooted-multicopy.nwk

# Generate single-copy gene trees with species-level labels
twig tree-filter \
    -i ${species_tree_dir}/tree-set-glabel-filtered-rooted.nwk \
    -d "|" \
    -di 0 \
    -c 1 \
    -I ${species_mapping_file} \
    -ri \
    > ${species_tree_dir}/tree-set-species-filtered-rooted-singlecopy.nwk

# Generate single-copy gene trees with accession-level labels
twig tree-filter \
    -i ${species_tree_dir}/tree-set-glabel-filtered-rooted.nwk \
    -d "|" \
    -di 0 \
    -rd \
    -c 1 \
    > ${species_tree_dir}/tree-set-accession-filtered-rooted-singlecopy.nwk

# Infer species tree from multicopy gene trees
astral-pro3 \
    -i ${species_tree_dir}/tree-set-species-filtered-rooted-multicopy.nwk \
    -o ${species_tree_dir}/sptree-species-filtered-rooted-multicopy.nwk \
    -t 10 \
    --root ${outgroup_species_accession} \
    --seed 42

# Infer species tree from single-copy gene trees (species-level labels)
astral-pro3 \
    -i ${species_tree_dir}/tree-set-species-filtered-rooted-singlecopy.nwk \
    -o ${species_tree_dir}/sptree-species-filtered-rooted-singlecopy.nwk \
    -t 10 \
    --root ${outgroup_species_accession} \
    --seed 42

# Infer species tree from single-copy gene trees (accession-level labels)
astral-pro3 \
    -i ${species_tree_dir}/tree-set-accession-filtered-rooted-singlecopy.nwk \
    -o ${species_tree_dir}/sptree-accession-filtered-rooted-singlecopy.nwk \
    -t 10 \
    --root ${outgroup_accession} \
    --seed 42
```


