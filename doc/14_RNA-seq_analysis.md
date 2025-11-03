### Analyzing gene expression in humans using RNA seq

#### Step 1. Download samples and split into R1 and R2
##### Sample A
```
module load SRA-Toolkit

SRR_NUMBER=SRR10677729
prefetch ${SRR_NUMBER}

# Convert SRA data to FASTQ files with gzip compression
fasterq-dump ${SRR_NUMBER} --split-files --threads 8 -O ./fastq_output
gzip ./fastq_output/${SRR_NUMBER}*.fastq

# Rename output files for clarity
mv ./fastq_output/${SRR_NUMBER}_1.fastq.gz sampleA_R1.fastq.gz
mv ./fastq_output/${SRR_NUMBER}_2.fastq.gz sampleA_R2.fastq.gz
```
##### Sample B
```
module load SRA-Toolkit

SRR_NUMBER=SRR10677731
prefetch ${SRR_NUMBER}

# Convert SRA data to FASTQ files with gzip compression
fasterq-dump ${SRR_NUMBER} --split-files --threads 8 -O ./fastq_output
gzip ./fastq_output/${SRR_NUMBER}*.fastq

# Rename output files for clarity
mv ./fastq_output/${SRR_NUMBER}_1.fastq.gz sampleB_R1.fastq.gz
mv ./fastq_output/${SRR_NUMBER}_2.fastq.gz sampleB_R2.fastq.gz
```

#### Step 2a. Download the latest human genome assembly
```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/\
Gencode_human/release_49/GRCh38.p14.genome.fa.gz &
# Takes ~30 min; consider prepending 'nohup'

gunzip -c GRCh38.p14.genome.fa.gz > GRCh38.p14.genome.fa
# '-c' keeps the input file
```

#### Step 2b. Download the annotations for the latest human genome assembly
```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/\
Gencode_human/release_49/gencode.v49.annotation.gtf.gz

gunzip -c gencode.v49.annotation.gtf.gz > gencode.v49.annotation.gtf
```

#### Step 3a. Perform quality control
```
module load FastQC

mkdir -p fastqc_output
fastqc -t 2 -f fastq -o fastqc_output *.fastq.gz \
1>fastqc_output/fastqc_runtime.log 2>&1 &
# '2>&' merges standard error into standard output
```

#### Step 3b. Summarize QC results

```
module load MultiQC

mkdir -p multiqc_output
multiqc -o multiqc_output ./fastqc_output \
1>multiqc_output/multiqc_runtime.log 2>&1 &
# '2>&' merges standard error into standard output
```

#### Step 3c. Check proper read pairing

```
module load cutadapt

cutadapt -o /dev/null -p /dev/null -j 2 \
sampleA_R1.fastq.gz sampleA_R2.fastq.gz \
1>checkPairing_sampleA_runtime.log 2>&1 &

cutadapt -o /dev/null -p /dev/null -j 2 \
sampleB_R1.fastq.gz sampleB_R2.fastq.gz \
1>checkPairing_sampleB_runtime.log 2>&1 &
```

##### Are reads properly paired? Check by counting reads:
```
zcat sampleA_R1.fastq.gz | grep "^@" | wc -l
zcat sampleA_R2.fastq.gz | grep "^@" | wc -l
```
##### Note that read 447 in sampleA_R1 is missing or misaligned in sampleA_R2:
```
zcat sampleA_R1.fastq.gz | head -n 1800 | tail -n 50
zcat sampleA_R2.fastq.gz | head -n 1800 | tail -n 50
```