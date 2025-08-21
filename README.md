# BioPrinCRISPR: A Biological Principle-Informed Framework for Discovery of Novel CRISPR-Cas Systems

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the source code and analysis pipeline for the manuscript: **"A deep learning and co-conservation framework enable discovery of non-canonical Cas proteins"**.

BioPrinCRISPR is a class-agnostic computational framework that integrates gene co-conservation, protein domain co-occurrence, and semantic embedding similarity to systematically identify and characterize novel CRISPR-Cas systems from prokaryotic genomes.

## Repository Structure

```
bioprincrispr/
├── LICENSE
├── README.md
├── config.yaml
├── data
│   └── reference
│       └── download_pfam.sh
├── environment.yml
├── notebooks
│   └── exploratory_analysis.ipynb
├── scripts
│   ├── 01_array.py
│   ├── 02_ap_pairs_v2.py
│   ├── 03_pairs_to_fasta.py
│   ├── 04_co-conservation_v4.py
│   ├── 05_domain_co-occurrence.py
│   ├── 05_domain_co-occurrence_3.py
│   ├── 05_domain_co-occurrence_undirected.py
│   ├── 06_rna_ss_vis.py
│   ├── 10_ArrowerSVG.py
│   ├── 10_ArrowerSVG_to_pdf.py
│   ├── 15_rank_combinations.py
│   ├── NetworkAnalysis.ipynb
│   ├── analysis
│   │   ├── plot_figure1.py
│   │   └── plot_figure5.py
│   └── embeddings
│       ├── README
│       ├── clustering_K-Means_OneStep.py
│       ├── draw_UMAP.py
│       └── get_embedding.py

```

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/bioprincrispr.git
    cd bioprincrispr
    ```

2.  **Create the Conda environment:**
     Create and activate the environment using the provided file.
    ```bash
    conda env create -f environment.yml
    conda activate crispr
    ```

3.  **Download Reference Databases:**
    Download the necessary reference databases (e.g., Pfam-A) by running the script in the `data/reference/` directory.
    ```bash
    bash data/reference/download_pfam.sh
    ```

## Quick Start / Usage

### Step 1: Identification of CRISPR Arrays using `01_array.py` 

**Description:**  
This script utilizes the `minced` tool to identify CRISPR repeat arrays within prokaryotic genome sequences provided as FASTA files. It supports multithreading to improve processing efficiency for large datasets. The input is a tab-separated file listing absolute paths to FASTA files. The user can customize detection criteria such as minimum repeats, repeat length, spacer length, and number of threads.

**Key Parameters:**

* `-i` / `--input`: Tab-separated file listing paths to FASTA files. *(Required)*
* `-rN` / `--repeat_number`: Minimum number of repeats. *(Required)*
* `-rMin`: Minimum repeat length (default 11).
* `-rMax`: Maximum repeat length (default 80).
* `-sMin`: Minimum spacer length (default 11).
* `-sMax`: Maximum spacer length (default 80).
* `-n` / `--num_threads`: Number of threads to use (default: all available).

**Example usage:**

```bash
python scripts/01_array.py \
  -i /data/genomes/file_list.tsv \
  -rN 3 \
  -rMin 11 \
  -rMax 80 \
  -sMin 11 \
  -sMax 80 \
  -n 16 \
  -o /data/output/
  
  example usage:
    python 01_array.py -i fna_path.txt -rN 3 -rMin 11 -rMax 80 -sMin 11 -sMax 80 -n 16 -o /out
```

**Output:**
* `outdir/repat_info` stores repeat info, and 
* `minced_out2fna.tsv` stores path info:
    
```bash
   outdir/
     ├── minced_out2fna.tsv
     ├── repeat_info/
     │   ├── fna_01.fa_minced_out.txt
     │   ├── fna_02.fa_minced_out.txt
     │   └── fna_03.fa_minced_out.txt

    minced_out2fna.tsv:
    fna_path       file_path
    /dada/user/David/fna_01.fa  outdir/repeat_info/fna_01.fa_minced_out.txt
    /dada/user/David/fna_02.fa  outdir/repeat_info/fna_02.fa_minced_out.txt
```

### Step 2: Generate Array-Protein Pairs using `02_ap_pairs_v2.py`

**Description:**  
Processes Step 1 output to associate CRISPR repeat arrays with nearby protein-coding sequences using Prodigal. Extracts relevant proteins with distances relative to arrays.

**Key Parameters:**

* `-i` / `--input`: TSV file linking FASTA files to minced output (`minced_out2fna.tsv`). *(Required)*
* `-w` / `--window`: Genomic window size around CRISPR array in base pairs. *(Required)*
* `-n` / `--num_threads`: Number of parallel threads (default: all available).
* `-o` / `--outdir`: Output directory. *(Required)*

**Example usage:**

```bash
python scripts/02_ap_pairs_v2.py \
  -i /data/output/minced_out2fna.tsv \
  -w 10000 \
  -n 16 \
  -o /data/output/
```
**Output:**
* `prodigal_faa/`
* `single_csv_and_faa/`
* `array2prot_pairs.csv`
```bash
  outdir/
 ├── array2prot_pairs.csv
 ├── prodigal_faa/
 ├── single_csv_and_faa
```

### Step 3: Extract and Deduplicate Array and Protein FASTA Sequences using `03_pairs_to_fasta.py`

**Description:**  
Converts array-protein pairs into FASTA format and performs deduplication. Supports filtering by repeat count, repeat length, and protein-array distance modes (`bp` or `orf`).

**Key Parameters:**

* `-i` / `--input`: Input CSV (`array2prot_pairs.csv`) from Step 2. *(Required)*
* `-mode` : Filter mode — `bp` (base pair distance) or `orf` (ORF count). *(Required)*
* `-rN` / `--repeat_number`: Minimum repeat copies (default 3).
* `-rL` / `--repeat_length`: Max repeat length (default 20).
* `-bp_num` / `--base_pair_num`: Max base pair distance for `bp` mode (default 8000).
* `-orf_num` / `--orf_num`: Max ORF distance for `orf` mode (default 5).
* `-o` / `--outdir`: Output directory. *(Required)*


**Example usage:**

```bash
python scripts/03_pairs_to_fasta.py \
  -i /data/output/array2prot_pairs.csv \
  -mode bp \
  -rN 3 \
  -rL 20 \
  -bp_num 8000 \
  -orf_num 5 \
  -o /data/output/
```
**Output**: 
* Filtered and deduplicated (on array) array and protein fasta:
```bash
    outdir/
   ├── array_fasta_all.csv
   ├── prot_fasta_all.csv
   ├── array_fasta_unique.csv
   ├── array_fasta_filter_by_5_orf_unique.csv  # we need this file
   ├── prot_fasta_filter_by_5_orf.csv   # we need this file
```
We need `array_fasta_filter_by_5_orf.csv` and `prot_fasta_filter_by_5_orf.csv` as input for the next step.

### Step 4: Co-Conservation Analysis Using `04_co-conservation_v4.py`

**Description:**
Performs co-conservation network analysis by integrating array and protein cluster graphs, filtering based on coverage and minimum cluster sizes, and identifying conserved protein sets.

First, run `mmseqs` to get prot and array clusters. 

```bash
    mmseqs for protein, using --min-seq-id = 0.3
    mmseqs easy-cluster --min-seq-id 0.3 --seq-id-mode 1 --cluster-mode 1 --cov-mode 1 --similarity-type 2 --single-step-clustering --dbtype 1 prot_fasta_filter_by_5_orf.csv 0.3_prot_out tmp
   
    mmseqs for array, using --min-seq-id = 0.35
    mmseqs easy-cluster --min-seq-id 0.35 --seq-id-mode 1 --cluster-mode 1 --cov-mode 1 --similarity-type 2 --single-step-clustering --dbtype 2 array_fasta_unique.csv 0.35_array_out tmp
                
         mmseq_cluster/
         ├── 0.3_prot_out_cluster.tsv  # we need this file
         ├── 0.35_array_out_cluster.tsv  # we need this file
         
```

Then, run `04_co-conservation_v4.py`
**Key Parameters:**

* `-a` / `--array_cluster`: TSV file for array clusters (edges). *(Required)*
* `-p` / `--prot_cluster`: TSV file for protein clusters (edges). *(Required)*
* `-c` / `--coverage`: Coverage threshold (0–1). *(Required)*
* `-m` / `--mininode_num`: Minimum cluster size (default 5).
* `-n` / `--num_threads`: Number of CPU threads (default: all available).
* `-o` / `--outdir`: Output directory. *(Required)*

**Example Usage:**

```bash
python scripts/04_co-conservation_v4.py \
  -a /data/output/0.55_array_out_cluster.tsv \
  -p /data/output/0.35_prot_out_cluster.tsv \
  -c 0.3 \
  -m 5 \
  -n 16 \
  -o /data/output/
```
**Output**: 
```bash
    output/
 ├── v1_array_kept.csv
 ├── v1_CLUSTER_INFO.txt
 ├── v1_prot_cluster.graphml
 ├── v1_prot_kept.csv  # This file contais proteins that are biologically assocoiated with array
 ├── v1_prot_representative_nodes.txt
```

### Step 5: Protein Domain Co-Occurrence Analysis using `05_domain_co-occurrence_3.py`

**Description:**
Constructs a directed domain-domain network from Pfam annotations, merging overlapping domains and filtering edges by protein counts and coverage.

**Key Parameters:**

* `-i` / `--input`: Domain annotation file (e.g., `hmmsearch -domtblout`). *(Required)*
* `-f` / `--prot_fasta`: Protein FASTA file. *(Required)*
* Optional:

  * `-r` / `--reference`: Reference domain annotation for Cas proteins.
  * `-a` / `--prot_id_focused`: File with protein IDs to restrict analysis.
  * `-p` / `--percent`: Minimum percentage threshold for edges (default 0).
  * `-m` / `--mini_edge_num`: Minimum edges count (default 0).
  * `-o` / `--outdir`: Output directory.

**Example Usage:**

```bash
python scripts/05_domain_co-occurrence_3.py \
  -i ./hmmsearch_domtblout.txt \
  -f ./proteins.fasta \
  -r ./reference_hmmout.txt \
  -a ./focused_protein_ids.txt \
  -p 0.02 \
  -m 20 \
  -o ./output_dir/
  
  python 05_domain_co-occurrence_3.py -i v1_prot_kept_hmmout.txt -f v1_prot_kept_fasta.fasta -r caspdb_and_caspedia_hmmout.txt -p 0 -m 5 -o ./1_domain_co-occurrence -a ./v1_prot_kept.csv

```

### Step 5 (Alternate): Undirected Protein Domain Co-Occurrence Analysis using `05_domain_co-occurrence_undirected.py`

**Description:**
Builds an undirected domain co-occurrence graph capturing all pairwise domain associations within proteins.

**Key Parameters:**

* `-i` / `--input`: Domain annotation file (e.g., `hmmsearch -domtblout`). *(Required)*
* `-f` / `--prot_fasta`: Protein FASTA file. *(Required)*
* Optional:

  * `-r` / `--reference`: Reference domain annotation.
  * `-a` / `--prot_id_focused`: Protein ID list for focus.
  * `-p` / `--percent`: Minimum percentage threshold (default 0).
  * `-m` / `--mini_edge_num`: Minimal proteins per domain pair (default 0).
  * `-o` / `--outdir`: Output directory.

**Example Usage:**

```bash
python scripts/05_domain_co-occurrence_undirected.py \
  -i ./hmmsearch_domtblout.txt \
  -f ./proteins.fasta \
  -r ./reference_hmmout.txt \
  -a ./focused_protein_ids.txt \
  -p 0.02 \
  -m 20 \
  -o ./output_dir/
```

### Step 6: RNA Secondary Structure Visualization using `06_rna_ss_vis.py`

**Description:**  
Visualizes RNA secondary structures using `matplotlib` and `forgi`, annotating structural elements with clash-free layout.

**Key Parameters:**  
This script does not accept command-line parameters. Input RNA sequence, secondary structure, output filename, and directory are hard-coded in the script variables and need to be modified directly in the code if customization is required.

**Example Usage:**  

```bash
python scripts/06_rna_ss_vis.py
```

### Step 6 (Alternate): RNA Secondary Structure Visualization using `06_rna_ss_vis_args.py`

**Description:**
Visualizes RNA secondary structures using `matplotlib` and `forgi`, annotating structural elements with clash-free layout. Inputs are RNA sequence and dot-bracket secondary structure files.

**Key Parameters:**

* `-s` / `--structure_file`: RNA secondary structure file (`.fx` format). *(Required)*
* `-q` / `--sequence`: RNA nucleotide sequence. *(Required)*
* `-i` / `--id`: RNA identifier for output naming. *(Required)*
* `-o` / `--outdir`: Output directory for SVG files. *(Required)*

**Example Usage:**

```bash
python scripts/06_rna_ss_vis_args.py \
  -s example_rna.fx \
  -q "GGCAACGGCGGCGGCAACGG..." \
  -i example_rna \
  -o ./results/
```

### Step 7: Gene Cluster Visualization using `10_ArrowerSVG.py`

**Description:**
Generates SVG visualizations of gene clusters from GenBank files, drawing arrows representing genes with strand orientation and coloring by Pfam domains. Annotates labels, repeat regions, and supports connecting conserved genes.

**Key Parameters:**
This script does not accept command-line parameters. Input directory, output paths, and drawing parameters are hard-coded inside the script and require manual modification prior to running.

**Example Usage:**

```bash
python scripts/10_ArrowerSVG.py
```

### Step 8: Ranking Protein Domain Combinations using `15_rank_combinations.py`

**Description:**
Extracts domain-domain combinations from a GraphML protein domain network, annotates Pfam pairs with descriptions, ranks them by protein count, and outputs a CSV summary.

**Key Parameters:**

* `-i` / `--input`: Input GraphML file of domain co-occurrence network. *(Required)*
* `-d` / `--desc`: Pfam accession-to-description mapping file (tab-delimited). *(Required)*
* `-o` / `--output`: Output CSV file for ranked domain pairs. *(Required)*

**Example Usage:**

```bash
python scripts/15_rank_combinations.py \
  -i ./domain_cooccurrence.graphml \
  -d ./pfam_acc2des.txt \
  -o ./ranked_domain_combinations.csv
```

## Citation

If you use BioPrinCRISPR or any findings from our study, please cite our manuscript:

```bibtex
@article{YourLastName2025BioPrinCRISPR,
  title={A deep learning and co-conservation framework enable discovery of non-canonical Cas proteins},
  author={B.B.H., C.Q.,Y.Y.F.,D.L.,F.W.,D.W.,Z.Y.,Y.Z.,H.X.L.,Y.Z.,Y.X.L. 
  volume={XX},
  pages={YY-ZZ}
}
```


