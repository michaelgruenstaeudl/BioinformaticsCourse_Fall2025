### Assembly of human mitochondrial genome

#### Step 1: Retrieving sequence reads of genome skimming or hybrid capture experiment on human genome
```
# Example: Biosample SAMN08193509
wget -c https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR006/664/SRR6664769.sralite.1

# Splitting raw reads file into R1 and, if present, R2 reads using FASTQ-DUMP
module load SRA-Toolkit
fastq-dump --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip SRR6664769.sralite.1
mv SRR6664769.sralite.1_pass_1.fastq.gz SRR6664769.sralite.1_R1.fastq.gz
```

#### Step 2: Retrieving human mitochondrial reference genome of similar sequence origin
```
# Example: NCBI accession number MG936619.1
efetch -db nucleotide -id MG936619 -format fasta > MG936619.fasta
efetch -db nucleotide -id MG936619 -format genbank > MG936619.gb
```

#### Step 3: Quality filtering of sequence reads using CUTADAPT
```
cutadapt -m 35 -q 30 -o SRR6664769.sralite.1_R1_QC.fastq.gz SRR6664769.sralite.1_R1.fastq.gz
```
