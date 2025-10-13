### Downloading short-read sequence data

#### Logging on to Beocat and installing edirect
```
# Log on to Beocat
ssh username@headnode.beocat.ksu.edu

# Installing edirect on Beocat
conda install bioconda::entrez-direct

# Installing some Perl dependencies on Beocat
module load Perl
module load Perl-bundle-CPAN
cpanm --local-lib=~/perl5 Time::HiRes
export PERL5LIB=$HOME/perl5/lib/perl5:$PERL5LIB
export PATH=$HOME/perl5/bin:$PATH
source ~/.bashrc
```

#### Download final genome and reference genomes
```
# Download genome in FASTA and GenBank format
efetch -db nuccore -id ${ACC_NUMBER} -format fasta \
  > complete_genome_${ACC_NUMBER}.fasta

efetch -db nuccore -id ${ACC_NUMBER} -format gb \
  > complete_genome_${ACC_NUMBER}.gb

esearch -db nucleotide -query "$RELATED_GENOMES" | 
	efetch -format fasta -stop 5 > related_genomes.fasta
```

#### Download short-read sequence data
```
# Check if SRA Toolkit is available:
module avail sra-toolkit
module load SRA-Toolkit		# Notice the capitalization !

# Confirming that it loaded correctly
fasterq-dump --version

# Download the SRA data into cache
prefetch ${SRR_NUMBER}

# Convert the SRA data to FASTQ with subsequent gzip compression
fasterq-dump ${SRR_NUMBER} --split-files --threads 8 -O ./fastq_output
gzip ./fastq_output/${SRR_NUMBER}*.fastq

# Note: This may take an hour or two!
```

### Quality control of short-read sequence data

#### Checking number of reads
```
ls -al fastq_output/${SRR_NUMBER}_?.fastq.gz

INF1=fastq_output/${SRR_NUMBER}_1.fastq.gz
INF2=fastq_output/${SRR_NUMBER}_2.fastq.gz

zcat $INF1 | grep "^@" | wc -l   # Will take 1-2 min!
zcat $INF2 | grep "^@" | wc -l   # Will take 1-2 min!
```

#### Look inside a FASTQ file
```
zcat $INF1 | head -n8
echo ""
zcat $INF2 | head -n8
```

#### Checking read lengths
```
zcat $INF1 | awk '{if(NR%4==2) print length($1)}' | \
sort | uniq -c     # Will take 2-3 min!
zcat $INF2 | awk '{if(NR%4==2) print length($1)}' | \
sort | uniq -c     # Will take 2-3 min!
```

#### Reducing read numbers
```
INFS1=${SRR_NUMBER}_small_1.fastq.gz
INFS2=${SRR_NUMBER}_small_2.fastq.gz

zcat $INF1 | head -n 2000000 | gzip > $INFS1
zcat $INF2 | head -n 2000000 | gzip > $INFS2

zcat $INFS1 | grep "^@" | wc -l
zcat $INFS2 | grep "^@" | wc -l
```

#### Checking read quality of FASTQ file
```
module load FastQC
mkdir -p fastqc_results
fastqc -o fastqc_results $INFS1 $INFS2

module load MultiQC
multiqc fastqc_results -o multiqc_report
```

#### Quality-filtering reads of FASTQ file
```
module load cutadapt

cutadapt -m 140 -q 30 \
-o ${SRR_NUMBER}_small_R1_QCchecked.fastq.gz \
-p ${SRR_NUMBER}_small_R2_QCchecked.fastq.gz \
$INFS1 $INFS2
```
