# A time-resolved multi-omics atlas of transcriptional regulation in response to high-altitude hypoxia across whole-body tissues

## 1. Introduction
This repository provides the scripts used in the paper **"A time-resolved multi-omics atlas of transcriptional regulation in response to high-altitude hypoxia across whole-body tissues"**. 

High-altitude hypoxia acclimatization requires whole-body physiological regulation in highland immigrants, but the underlying genetic mechanism has not been clarified. Here we used sheep as an animal model for plain-to-plateau translocationplantation. We generated multi-omics data including **49 whole genome sequencing, 1,277 time-resolved bulk RNA-Seq, 66 ATAC-Seq and 6 single-cell RNA-Seq datasets** from multiple tissues as well as **phenotypic data from 20 bio-indicators**. Our study provides valuable resources and insights for future hypoxia-related studies in mammals.

Our work is out now on **Nature Communications**: https://doi.org/10.1038/s41467-024-48261-w.


## 2. The main contents
### 2.1 Codes:
There are five main parts of analysis in our study, including **WGS, RNA-Seq, scRNA-Seq, ATAC-Seq and Hi-C**. Scripts for **scRNA-Seq** data analysis have been uploaded to https://github.com/tian6067/WeiwtCAU-sc-RNA_seq.

### 2.2 Genome files:
Genome reference and annotation files used in this study are downloaded from Ensembl genome browser: https://asia.ensembl.org/Ovis_aries_rambouillet/Info/Index (i.e., **Oar_rambouillet_v1.0**).

### 2.3 Data table
To enhance the reusability and reproducibility of our data, and fully facilitate them for secondary analysis and novel discoveries, we uploaded **"Data table.xlsx"** which includes **sample collection, information for sequencing sample, and metatables for sequencing datasets**. You can find the details about newly generated sequencing data and public data in this study directly.

### 2.4 Result files
We provide **integrated expression matrix, consensus peak set with metatables and measurement phenotypic data** in our study, which have been uploaded to Google Drive: https://drive.google.com/file/d/1iTbn5IVPGVGEp7nkeNCgjDvtqjXjPQQ8/view?usp=sharing.


## 3.Data availability 
High-throughput sequencing data in this study are publicly available for download without restriction, which was deposited in the Sequence Read Archive (SRA) database in NCBI under accession numbers: **PRJNA1053506 (WGS), PRJNA1000743 (bulk RNA-Seq), PRJNA1001016 (bulk RNA-Seq) and PRJNA1001505 (scRNA-Seq, ATAC-Seq).** 

## 4. Citation
Yan, Z., Yang, J., Wei, WT. et al. A time-resolved multi-omics atlas of transcriptional regulation in response to high-altitude hypoxia across whole-body tissues. **Nat. Commun.** 15, 3970 (2024).

Any questions, please feel free to contact yanze@cau.edu.cn.
