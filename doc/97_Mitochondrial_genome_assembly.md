### Assembly of human mitochondrial genome

#### Step 1: Retrieving sequence reads of genome skimming or hybrid capture experiment on human genome
```
# Example: Biosample SAMN08193509 - SRA number SRR6664769
# Note: You can check if the sequence reads match your aims by evaluating the title of the SRR-file
esearch -db sra -query "SRR6664769" | esummary | xtract -pattern DocumentSummary -element Title

# Downloading raw sequence reads
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

#### Step 4: Mapping the reads to the human mitohcondrial reference genome using BOWTIE2 and SAMTOOLS
```
INF=SRR6664769.sralite.1_R1_QC.fastq.gz
REFGENOME=MG936619.fasta
mkdir -p db 
bowtie2-build $REFGENOME db/myRef > refdb.log 
bowtie2 -x db/myRef -U $INF 2> mapping.err | samtools view -bS > mapping.bam
samtools flagstat mapping.bam > mapping_stats.txt
```

#### Step 5: Extracting the mapped reads using SAMTOOLS and BEDTOOLS
The mapped reads will represent the mitochondrial reads of theinput sequence read set
```
samtools view -b -F 0 -F 16 mapping.bam | samtools sort -n - > extracted_mappedReads.bam  
bedtools bamtofastq -i extracted_mappedReads.bam -fq mapped_reads.fastq
```

#### Step 6: Assemblying the mitochondrial genome suing SPADES
```
# Note: Option '-s' stands for single-end reads
mkdir -p assembly_results
spades.py --isolate -s mapped_reads.fastq -o assembly_results -t 4
```
