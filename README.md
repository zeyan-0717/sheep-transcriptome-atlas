# A time-resolved multi-omics atlas of transcriptional regulation in response to high-altitude hypoxia across whole-body tissues

## 1. Introduction
This repository provides the scripts used in the paper **A time-resolved multi-omics atlas of transcriptional regulation in response to high-altitude hypoxia across whole-body tissues**. 

High-altitude hypoxia acclimatization requires whole-body physiological regulation in highland immigrants, but the underlying genetic mechanism has not been clarified. Here we used sheep as an animal model for plain-to-plateau translocationplantation. We generated multi-omics data including **49 whole genome sequencing, 1,277 time-resolved bulk RNA-Seq, 66 ATAC-Seq and 6 single-cell RNA-Seq datasets** from multiple tissues as well as **phenotypic data from 20 bio-indicators**. We characterized transcriptional changes of all genes in each tissue, and examined multi-tissue temporal dynamics and transcriptional interactions among genes. In particular, we identified critical functional genes regulating the short response to hypoxia in each tissue. We further identified TAD-constrained cis-regulatory elements, which suppressed the transcriptional activity of most genes under hypoxia. Phenotypic and transcriptional evidence indicated that antenatal hypoxia could improve hypoxia tolerance in offspring. Furthermore, we provided time-series expression data of candidate genes associated with human mountain sickness and high-altitude adaptation.**Our study provides valuable resources and insights for future hypoxia-related studies in mammals**.

The preprint is now on biorxiv: https://www.biorxiv.org/content/10.1101/2023.10.25.563964v1.


## 2. The main contents
### 2.1 Codes:
There are five main parts of analysis in our study, including **WGS, RNA-Seq, scRNA-Seq, ATAC-Seq and Hi-C**. Scripts for scRNA-Seq data analysis have been uploaded to https://github.com/tian6067/WeiwtCAU-sc-RNA_seq.

### 2.2 Genome files:
Genome reference and annotation files used in this study are downloaded from Ensembl genome browser: https://asia.ensembl.org/Ovis_aries_rambouillet/Info/Index (i.e., **Oar_rambouillet_v1.0**).

### 2.3 Data table
To enhance the reusability and reproducibility of our data, and fully facilitate them for secondary analysis and novel discoveries, we uploaded **"Data table.xlsx"** which include **sample collection, information for sequencing sample, and metatables for sequencing datasets**. You can clearly find the details about newly generated sequencing data and public data in this study.

### 2.4 Result files
We provide **integrated expression matrix, consensus peak set with metatables and measurement phenotypic data** in our study, which have been uploaded to Google Drive: https://drive.google.com/file/d/1iTbn5IVPGVGEp7nkeNCgjDvtqjXjPQQ8/view?usp=sharing.


## 3.Data availability 
High-throughput sequencing data in this study are publicly available for download without restriction, which was deposited in the Sequence Read Archive (SRA) database in NCBI under accession numbers: **PRJNA1053506 (WGS), PRJNA1000743 (bulk RNA-Seq), PRJNA1001016 (bulk RNA-Seq) and PRJNA1001505 (scRNA-Seq, ATAC-Seq).** 

Any questions, please feel free to contact yanze@cau.edu.cn.
