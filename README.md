# BioPrinCRISPR: A Biological Principle-Informed Framework for Discovery of Novel CRISPR-Cas Systems

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the source code and analysis pipeline for the manuscript: **"Systematic discovery across one million genomes reveals vast architectural and functional diversity of CRISPR-Cas systems"**.

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
└── src
    └── bioprincrispr
        ├── __init__.py
        ├── conservation.py
        ├── discovery.py
        └── utils.py
```

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/bioprincrispr.git
    cd bioprincrispr
    ```

2.  **Create the Conda environment:**
    All dependencies are managed using Conda. Create and activate the environment using the provided file.
    ```bash
    conda env create -f environment.yml
    conda activate bioprincrispr
    ```

3.  **Download Reference Databases:**
    Download the necessary reference databases (e.g., Pfam-A) by running the script in the `data/reference/` directory.
    ```bash
    bash data/reference/download_pfam.sh
    ```

## Quick Start / Usage

The entire pipeline can be run using the scripts in the `scripts/` directory. The order is indicated by the numerical prefix.

1.  **Configure Paths:** Ensure all paths in `config.yaml` are set correctly for your system.
2.  **Run the pipeline:** Execute the main workflow script (you can create a master script `run_pipeline.sh` that calls the numbered scripts in order).

    For example, to run the first step on the example data:
    ```bash
    python scripts/01_run_minced.py --config config.yaml
    ```

## Citation

If you use BioPrinCRISPR or any findings from our study, please cite our manuscript:

```bibtex
@article{YourLastName2025BioPrinCRISPR,
  title={Systematic discovery across one million genomes reveals vast architectural and functional diversity of CRISPR-Cas systems},
  author={Your Name, et al.},
  journal={Journal Name},
  year={2025},
  volume={XX},
  pages={YY-ZZ}
}
```
