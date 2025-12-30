### Analyzing gene expression in humans using RNA seq

#### Step 1. Download samples and split into R1 and R2
##### Sample A
```bash
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
```bash
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
```bash
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/\
Gencode_human/release_49/GRCh38.p14.genome.fa.gz &
# Takes ~30 min; consider prepending 'nohup'

gunzip -c GRCh38.p14.genome.fa.gz > GRCh38.p14.genome.fa
# '-c' keeps the input file
```

#### Step 2b. Download the annotations for the latest human genome assembly
```bash
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/\
Gencode_human/release_49/gencode.v49.annotation.gtf.gz

gunzip -c gencode.v49.annotation.gtf.gz > gencode.v49.annotation.gtf
```

#### Step 3a. Perform quality control
```bash
module load FastQC

mkdir -p fastqc_output
fastqc -t 2 -f fastq -o fastqc_output *.fastq.gz \
1>fastqc_output/fastqc_runtime.log 2>&1 &
# '2>&' merges standard error into standard output
```

#### Step 3b. Summarize QC results

```bash
module load MultiQC

mkdir -p multiqc_output
multiqc -o multiqc_output ./fastqc_output \
1>multiqc_output/multiqc_runtime.log 2>&1 &
# '2>&' merges standard error into standard output
```

#### Step 3c. Check proper read pairing

```bash
module load cutadapt

cutadapt -o /dev/null -p /dev/null -j 2 \
sampleA_R1.fastq.gz sampleA_R2.fastq.gz \
1>checkPairing_sampleA_runtime.log 2>&1 &

cutadapt -o /dev/null -p /dev/null -j 2 \
sampleB_R1.fastq.gz sampleB_R2.fastq.gz \
1>checkPairing_sampleB_runtime.log 2>&1 &
```

##### Are reads properly paired? Check by counting reads:
```bash
zcat sampleA_R1.fastq.gz | grep "^@" | wc -l
zcat sampleA_R2.fastq.gz | grep "^@" | wc -l
```
##### Note that read 447 in sampleA_R1 is missing or misaligned in sampleA_R2:
```bash
zcat sampleA_R1.fastq.gz | head -n 1800 | tail -n 50
zcat sampleA_R2.fastq.gz | head -n 1800 | tail -n 50
```

#### Step 4a. Ensure proper pairing with Trimmomatic
```bash
module load Trimmomatic

INF1=sampleA_R1.fastq.gz
INF2=sampleA_R2.fastq.gz
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
-threads 4 -phred33 $INF1 $INF2 \
${INF1%.fastq.gz}_paired.fastq.gz ${INF1%.fastq.gz}_unpaired.fastq.gz \
${INF2%.fastq.gz}_paired.fastq.gz ${INF2%.fastq.gz}_unpaired.fastq.gz \
MINLEN:25
```

#### Step 4b. Trim reads and remove adapters with cutadapt

##### Adapter sequences:
```
R1:  5'-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA-3'
R2:  5'-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3'
```

##### Cutadapt commands:
```bash
module load cutadapt

cutadapt -m 22 -O 4 -j 2 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-o sampleA_R1_trimmed.fastq.gz sampleA_R1_paired.fastq.gz \
1>cutadapt_sampleA_R1.log 2>&1 &

cutadapt -m 22 -O 4 -j 2 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o sampleA_R2_trimmed.fastq.gz sampleA_R2_paired.fastq.gz \
1>cutadapt_sampleA_R2.log 2>&1 &
```

#### Step 4c. File hygiene
```bash
rm SRR106777??.lite.1_pass_?.fastq.gz

rm sample?_R?_paired.fastq.gz

rm sample?_R?_unpaired.fastq.gz
```

#### Step 4d. Reduce number of input reads to 5 million per sample
```bash
## CURRENT READS
zcat sampleA_R1_trimmed.fastq.gz | grep "^@" | wc -l
zcat sampleA_R2_trimmed.fastq.gz | grep "^@" | wc -l

## REDUCE TO 5 MILLION READS
zcat sampleA_R1_trimmed.fastq.gz | head -n 20000000 | \
gzip > sampleA_R1_trimmed_subset5M.fastq.gz &
zcat sampleA_R1_trimmed_subset5M.fastq.gz | grep "^@" | wc -l

zcat sampleA_R2_trimmed.fastq.gz | head -n 20000000 | \
gzip > sampleA_R2_trimmed_subset5M.fastq.gz &
zcat sampleA_R2_trimmed_subset5M.fastq.gz | grep "^@" | wc -l
```


#### Step 5b. Split human genome by chromosome
```bash
mkdir chromosomes

csplit --digits=2 --prefix=chromosomes/GRCh38.p14.genome.chromosome \
GRCh38.p14.genome.fa "/>chr/" "{*}"

# Remove empty first split
rm chromosomes/GRCh38.p14.genome.chromosome00

# Rename chromosomes 23-25
mv chromosomes/GRCh38.p14.genome.chromosome23 chromosomes/GRCh38.p14.genome.chromosomeX
mv chromosomes/GRCh38.p14.genome.chromosome24 chromosomes/GRCh38.p14.genome.chromosomeY
mv chromosomes/GRCh38.p14.genome.chromosome25 chromosomes/GRCh38.p14.genome.chromosomeMT
```


#### Step 5c. Extract genome annotations for specific chromosomes
```bash
grep "^chr8" gencode.v49.annotation.gtf \
> gencode.v49.annotation_chr8.gtf
```


#### Step 5d. Building index for specific chromosome with STAR
```bash
module load STAR

# Build index
mkdir -p GRCh38_chr8_index
STAR --runThreadN 2 \
--runMode genomeGenerate \
--genomeDir GRCh38_chr8_index \
--genomeFastaFiles chromosomes/GRCh38.p14.genome.chromosome08 \
--sjdbGTFfile gencode.v49.annotation_chr8.gtf \
--genomeSAindexNbases 12  \
--sjdbOverhang 149 --outFileNamePrefix chr8 \
1>STAR_buildingIndex_runtime.log 2>&1 &
```

#### Step 5e. Map reads against specific chromosome with STAR
```bash
STAR --runThreadN 4 \
--genomeDir GRCh38_chr8_index \
--readFilesIn sampleA_R1_trimmed_subset5M.fastq.gz sampleA_R2_trimmed_subset5M.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix GRCh38_chr8_mapping/sampleA_ \
1>STAR_mapping_sampleA_runtime.log 2>&1 &
```

#### Step 5f. File hygiene
```bash
cp GRCh38_chr8_mapping/sampleA_Aligned.sortedByCoord.out.bam \
GRCh38_chr8_mapping/sampleA_mapping.bam
```

#### Step 6. Remove technical duplicates from mapping file
```bash
module load picard

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=GRCh38_chr8_mapping/sampleA_mapping.bam \
O=GRCh38_chr8_mapping/sampleA_mapping_techDuplRemoved.bam \
M=GRCh38_chr8_mapping/sampleA_marked_dupl_metrics.txt \
REMOVE_DUPLICATES=true
```

#### Step 7a. Count unique reads per gene
```bash
module load Subread

featureCounts -T 2 -s 2 -p \
-a gencode.v49.annotation_chr8.gtf \
-o GRCh38_chr8_sampleA_featCounts.txt \
GRCh38_chr8_mapping/sampleA_mapping_techDuplRemoved.bam \
1>featCounts_sampleA.log 2>&1 &
```

#### Step 7b. Report count for gene of interest
```bash
EnsemblID="ENSG00000104237"
# All information
grep $EnsemblID GRCh38_chr8_sampleA_featCounts.txt

# Count of mapped reads
grep $EnsemblID GRCh38_chr8_sampleA_featCounts.txt | \
   awk -F'\t' '{print $NF}'
```
