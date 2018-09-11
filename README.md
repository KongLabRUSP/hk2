##  Project: Study of Diabetic Nephropathy in HK2 Kidney Cells
### Study ID: 
### Scientist/Data Analysis: David Cheng, Renyi Wu, Davit Sargsyan 
### Created: 05/23/2018

---    

## Table of Contents
[File Legend](#leg)
[Daily Logs](#logs)  

## File Legend<a name="files"></a>
***/source/hk2_rnaseq_degseq_v1.R***: current (06/10/2018) analysis script      
***/source/hk2_dn_rnaseq_allignment_hisat2.txt***: current(06/10/2018) alignment script    
***/docs/GENEWIZ NGS Sample Submission Form 2-1-18 RNA.xlsx***: sample legend     
***/docs/hk2_tiia_rnaseq_v1.pptx***: current (06/10/2018) results        
***/data/featurecounts.results.human.csv***: current (06/10/2018) expressions table     
***/share/Renyi/RNA_gw/human/featurecounts.results.human.csv***: location of original RNA-seq raw counts file. Only need 3 samples: HG.dedup.bam, LG.dedup.bam and MIC1.dedup.bam    
***/datastorage/FastQ_2018/RNA/April***: pair-ended FastQ files (i.e. 2 data files per sample): LG_R1_001.fastq.gz and LG_R2_001.fastq.gz, HG_R1_001.fastq.gz and HG_R2_001.fastq.gz, and MIC1_R1_001.fastq.gz and MIC1_R2_001.fastq.gz       
***/datastorage/Processed_BAM_Files/Renyi/RNA_gw/human***: processed .BAM files LG.sorted.bam, HG.sorted.bam and MIC1.sorted.bam      

## Daily Logs<a name="logs"></a>
### 09/10/2018
* Relabeled samples correctly - LG, HG and MITC, in that order. Rerun the code.

### 06/10/2018
* Added donut hitmap with clustering

### 05/31/2018
* Adopted MES13 MITC script for HK2