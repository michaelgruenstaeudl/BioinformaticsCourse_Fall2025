### Analyzing gene expression in humans using RNA seq

#### Step 1. Download samples and split into R1 and R2
##### Sample A
```
module load SRA-Toolkit

wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/\
sra-pub-zq-16/SRR010/10677/SRR10677729/SRR10677729.lite.1

fastq-dump --split-files --gzip --skip-technical --read-filter pass \
--origfmt --readids --clip SRR10677729.lite.1

mv SRR10677729.lite.1_pass_1.fastq.gz sampleA_R1.fastq.gz
mv SRR10677729.lite.1_pass_2.fastq.gz sampleA_R2.fastq.gz
```
##### Sample B
```
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/\
sra-pub-zq-16/SRR010/10677/SRR10677731/SRR10677731.lite.1

fastq-dump --split-files --gzip --skip-technical --read-filter pass \
--origfmt --readids --clip SRR10677731.lite.1

mv SRR10677731.lite.1_pass_1.fastq.gz sampleB_R1.fastq.gz
mv SRR10677731.lite.1_pass_2.fastq.gz sampleB_R2.fastq.gz
```

#### Step 2a. Download the latest human genome assembly
```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/\
Gencode_human/release_40/GRCh38.p13.genome.fa.gz &
# Takes ~30 min; consider prepending 'nohup'

gunzip -c GRCh38.p13.genome.fa.gz > GRCh38.p13.genome.fa
# '-c' keeps the input file
```

#### Step 2b. Download the annotations for the latest human genome assembly
```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/\
Gencode_human/release_40/gencode.v40.annotation.gtf.gz

gunzip -c gencode.v40.annotation.gtf.gz > gencode.v40.annotation.gtf
```
