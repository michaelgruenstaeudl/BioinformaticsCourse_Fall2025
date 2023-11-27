### Analyzing gene expression in humans using RNA seq
Aim: to infer the expression of a gene of interest in humans under different conditions

#### Overview
1. Download two samples of human RNA-Seq data
2. Download latest assembly of the human genome plus its sequence annotations
3. Perform quality control with 'FastQC' and 'MultiQC'
4. Perform trimming and adapter removal with 'cutadapt'
5. Map reads against specific chromosome of human reference genome with 'STAR'
6. Remove duplicate reads
7. Count number of mapping reads

#### Step 1: Download two samples of human RNA-Seq data
Download the paired-end RNA-Seq from NCBI SRA and split each sample into R1 and R2
Two samples have been selected:
- [SRR10677729 (normoxia state)](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10677729)
- [SRR10677731 (hypoxia state)](https://www.ncbi.nlm.nih.gov/sra/?term=SRR10677731)_
```
module load SRA-Toolkit

# Sample A (i.e., normoxia state)
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/\
sra-pub-zq-16/SRR010/10677/SRR10677729/SRR10677729.lite.1
fastq-dump --split-files --gzip --skip-technical --read-filter pass \
--origfmt --readids --clip SRR10677729.lite.1
mv SRR10677729.lite.1_pass_1.fastq.gz sampleA_R1.fastq.gz
mv SRR10677729.lite.1_pass_2.fastq.gz sampleA_R2.fastq.gz

# Sample B (i.e., hypoxia state)
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/\
sra-pub-zq-16/SRR010/10677/SRR10677731/SRR10677731.lite.1
fastq-dump --split-files --gzip --skip-technical --read-filter pass \
--origfmt --readids --clip SRR10677731.lite.1
mv SRR10677731.lite.1_pass_1.fastq.gz sampleB_R1.fastq.gz
mv SRR10677731.lite.1_pass_2.fastq.gz sampleB_R2.fastq.gz
```
