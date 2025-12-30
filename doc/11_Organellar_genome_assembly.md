### Introduction to organellar genome assembly

#### Extracting organellar genome reads -- Mapping and extracting mapped reads
```bash
# Load Bowtie2 and SAMtools
module load Bowtie2
module load SAMtools

# Input FASTQ files and reference genome
INF1=${SRR_NUMBER}_small_1.fastq.gz
INF2=${SRR_NUMBER}_small_2.fastq.gz
REF=related_genomes_plastid.fasta
INDEX=db/myRef

# Build Bowtie2 index
mkdir -p db
bowtie2-build $REF $INDEX > refdb.log 2>&1

# Map reads and generate sorted BAM file
bowtie2 -x $INDEX -1 $INF1 -2 $INF2 2> mapping.err | \
samtools view -bS - | samtools sort -o mapping.sorted.bam
samtools index mapping.sorted.bam

# Generate mapping statistics
samtools flagstat mapping.sorted.bam > mapping_stats.txt

# Select properly paired mapped reads and convert to FASTQ
# -F 12: exclude unmapped reads and mates
samtools view -b -F 12 mapping.sorted.bam | \
samtools sort -n -o mapped_F12_namesorted.bam

# Convert directly to paired FASTQ files
samtools fastq \
  -1 mapped_R1.fastq \
  -2 mapped_R2.fastq \
  mapped_F12_namesorted.bam

# The files mapped_R1.fastq and mapped_R2.fastq are ready for downstream use.
```

#### Quality check
```bash
# Count mapped read pairs:
grep "with itself and mate mapped" mapping_stats.txt | \
awk '{print $1}'

# Count reads in paired FASTQ files:
# Each read starts with '@' on the header line
grep "^@" mapped_R1.fastq | wc -l
grep "^@" mapped_R2.fastq | wc -l
```

#### Plastid genome assembly with SPAdes
```bash
# Start the assembly:
module load SPAdes
mkdir -p SPAdes_results
spades.py --isolate -1 mapped_R1.fastq -2 mapped_R2.fastq -o SPAdes_results -t 8

# Evaluate the output of the assembly:
head -n100 SPAdes_results/contigs.fasta
grep "^>" SPAdes_results/contigs.fasta
```

#### Installing and testing GetOrganelle
```bash
# Installing GetOrganelle - Option 1
module unload Python  # GetOrganelle needs its own version of Python
conda create -n getorganelle_env -c conda-forge python=3.9 requests -y
conda activate getorganelle_env
conda install -c bioconda getorganelle -y

# Installing GetOrganelle - Option 2
module unload Python  # GetOrganelle needs its own version of Python
conda create -n getorganelle_env -c conda-forge -c bioconda python=3.10 getorganelle -y
conda activate getorganelle_env

# Testing if GetOrganelle works
get_organelle_from_reads.py --help

# Downloading GetOrganelle databases
module load BLAST+
get_organelle_config.py -a all
```

#### Plastid genome assembly with GetOrganelle
```bash
# Run GetOrganelle on filtered reads
get_organelle_from_reads.py \
  -1 mapped_R1.fastq -2 mapped_R2.fastq \
  -o getorganelle_results \
  -F embplant_pt -R 15 -k 21,45,65,85,105 \
  -t 8

# Key options:
# -F embplant_pt    # Target: plastid genome of seed plants
# -R 15             # Maximum extension rounds
# -k 21,45,65,85,105 # K-mer sizes for iterative assembly
# -t 8              # Number of CPU threads
```

#### Evaluate GetOrganelle output
```bash
# Evaluate the output of the assembly:
head -n100 getorganelle_results/extended_spades/contigs.fasta
grep "^>" getorganelle_results/extended_spades/contigs.fasta
```
