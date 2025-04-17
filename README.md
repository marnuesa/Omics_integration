# 🧬 Omics Integration in *Cucumis melo*

This repository contains all the scripts and resources used for a multi-omics analysis of stress response in melon plants (*Cucumis melo*, cultivar Piel de Sapo). The pipeline integrates small and long RNA transcriptomics, whole-genome bisulfite sequencing (WGBS) methylome data, and miRNA–mRNA–methylation associations using multiple linear regression.

---

## 📁 Project Structure

The pipeline includes several analysis stages:

### 1. Sample Collection

- Samples were collected at 3 time points (2, 4, and 11 days post-treatment).
- Five conditions: control, one biotic stress (*Monosporascus cannonballus*), and three abiotic stresses (cold, drought, short day).

### 2. Small RNA Transcriptome (miRNAs)

- **Preprocessing**: `fastp` + `FastQC`
- **Filtering of non-coding RNAs**: `sRNA_counts.py` (alignment with `Bowtie`)
- **Quantification**: manual count + `sqlite3`
- **Differential expression**: `DESeq2`
- **miRNA identification**: `miRNAs_Annotation.sh`, `Group_miRNA_by_family.sh` (alignment with `Bowtie`)

### 3. Long RNA Transcriptome (genes)

- **Preprocessing**: `fastp` + `FastQC`
- **Quantification**: `Salmon`
- **Differential expression**: `tximport` + `DESeq2`
- **GO enrichment**: `ClusterProfiler` using custom `org.CMelo.eg.db` annotation

### 4. Methylome (WGBS)

- **Preprocessing**: `TrimGalore` + `FastQC`
- **Mapping and methylation calling**: `Bismark`
- **DMR identification**: `DMRcaller`
- **Annotation**: `Bedtools` and custom Bash scripts

### 5. Multi-omics Integration

- **Modeling**: `MORE` package
- **Integration via multiple linear regression** using normalized gene expression, miRNA counts, and methylation values.

---

## ⚙️ Software and Versions

| Tool / Package             | Version     |
|----------------------------|-------------|
| Python                     | 3.6.8       |
| fastp                      | 0.23.2      |
| FastQC                     | 0.12.1      |
| Bowtie                     | 1.3.1       |
| sqlite3                    | 3.26.0      |
| R                          | 4.3.3       |
| DESeq2                     | 1.42.1      |
| tximport                   | 1.30.0      |
| ClusterProfiler            | 4.6.2       |
| AnnotationForge            | 1.40.2      |
| Salmon                     | 1.6.0       |
| TrimGalore (includes Cutadapt)| 0.6.10 |
| Cutadapt                   | 4.6         |
| Bismark                    | 0.24.2      |
| DMRcaller                  | 1.34.0      |
| Bedtools                   | 2.30.0      |
| ComplexHeatmap             | 2.14.0      |
| Cytoscape                  | 3.10.2      |
| ggplot2                    | 3.5.1       |
| methylKit                  | 1.28.0      |
| Genomation                 | 1.34.0      |
| MORE                       | 1.0         |

---
