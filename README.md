
# AviTag-seq

AviTag-seq is a unified bioinformatics workflow for regulatory-grade safety evaluation of gene-editing therapies, 
enabling simultaneous, high-sensitivity profiling of nuclease/editor-induced off-target events and AAV vector integration sites from the same sequencing dataset.

---

## Workflow Overview

<img width="1326" height="477" alt="image" src="https://github.com/user-attachments/assets/544a31e6-ea65-4897-8625-aebfac7924b7" />

---

## Environment

Please install dependencies manually (Conda recommended).

### System Requirements
- OS: Linux (recommended)
- Python: >= 3.7.12 
- Conda: Miniconda / Mambaforge (recommended)

### Python Dependencies
- numpy == 1.21.6
- pyfaidx == 0.8.1.1
- svgwrite == 1.4.3
- swalign == 0.3.7
- pyyaml == 6.0

---

### Recommended Manual Installation (example)

Create a clean environment:
```bash
conda create -n ENV_NAME python=3.10 -y
conda activate ENV_NAME
```
---

## Inputs

### Required files

**Reference**
- `data/ref/reference.fa` (FASTA)  
  Supported builds: **GRCh38/hg38** (human) or **GRCm38/mm10** (mouse)
- `data/ref/annotation.gff3` or `data/ref/annotation.gtf` (optional)  
  Used for genomic feature annotation of integration sites (e.g., promoter/exon/intron/UTR/intergenic).

**Samples**
- Paired-end FASTQ files generated from **Illumina NovaSeq PE150** sequencing:  
  `data/reads/<sample>_R1.fastq.gz` and `data/reads/<sample>_R2.fastq.gz`

### Input assumptions

- Sequencing mode: **Illumina paired-end (PE150)**
- Reference genome: **GRCh38/hg38** or **GRCm38/mm10**
- Read preprocessing (recommended): adapter trimming and quality filtering (e.g., `fastp`)
- Normalization (recommended): optional depth normalization via downsampling (e.g., to a fixed number of read pairs per sample) to reduce sequencing-depth–dependent bias across libraries
- Alignment: high-quality mapping is required for robust site calling (e.g., MAPQ-based filtering)

---

# optional analysis settings 
min_mapq: 50
merge_window_bp: 10
max_sgRNA_mismatches: 6


---

## Outputs

All outputs are written to `results/` (or the directory specified by `outdir`).

**Typical outputs**
- `results/qc/` — read QC summaries (e.g., FastQC/MultiQC) and preprocessing logs  
- `results/alignment/` — alignment BAM/BAI files and mapping statistics  
- `results/offtarget/` — de novo off-target site calls, sgRNA-homology annotations, and per-site summary tables  
- `results/integration/` — AAV integration breakpoint tables and genomic feature annotations (e.g., promoter/exon/intron/UTR/intergenic)  
- `results/plots/` — publication-ready plots (PDF/PNG)  
- `results/report/` — final summary tables and analysis reports  

**Paper-facing deliverables (edit as needed)**
- Figures: `results/plots/*.pdf`  
- Tables: `results/report/*.tsv`  


---
## How to Use

This repository contains two analysis branches:

- **Off-target profiling (AviTag-seq)**  
  The `AviTag-seq/` directory provides the bioinformatics workflow for **nuclease/editor-induced off-target detection** from UMI-tagged capture sequencing data.

- **AAV integration site profiling (AAV-Integrate-Seq)**  
  The `AAV-Integrate-Seq-human/` and `AAV-Integrate-Seq-mouse/` directories provide workflows for **AAV integration breakpoint calling and genomic feature annotation**, using **hg38** (human) and **mm10** (mouse) references, respectively.

### Run AviTag-seq (off-target calling)

```bash
python AviTag-seq/AviTag-seq.py all -m xxx.yaml
```
### Run AAV-Integrate-Seq
```bash
python AAV-Integrate-Seq-human/AAV-Integrate-Seq.py all -m xxx.yaml
```

```
```
