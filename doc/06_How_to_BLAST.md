## Introduction to BLAST

### Installing BLAST and seqkit
```
# On WSL / Linux:
sudo apt install ncbi-blast+
sudo apt install seqkit

# On Mac:
brew install blast
brew install seqkit
```


### Extracting the spike glycoprotein from the SARS-CoV-2 genome
```
# Step 1: Download the FASTA file of SARS-CoV-2 record 'NC_045512.2'

# Step 2: Extract the spike glycoprotein from the FASTA 
# sequence of the SARS-CoV-2 genome

# Option 1: Using seqkit
seqkit subseq -r 21563:25384 NC_045512.2.fasta > spike_glycoprotein.fasta

# Option 2: Using EMBOSS
extractseq -sequence NC_045512.2.fasta -start 21563 \
  -end 25384 -outseq spike_gene.fasta

# Looking at the DNA sequence
cat spike_glycoprotein.fasta
```


### Using BLAST to analyze a SARS-CoV-2 genome fragment
```
# Step 3: Run a basic BLAST search locally or remotely
blastn -query spike_glycoprotein.fasta -db nt -remote \
  -outfmt '6 std stitle' \
  -out BLAST_results.tsv
# This will take a minute or two

# blastn: nucleotide BLAST
# -query: input sequence
# -db nt: search NCBI nucleotide database
# -remote: run on NCBI’s servers
# -out BLAST_results.tsv: save results
# -outfmt '6 std stitle': tabular output, including hit title
```


### Evaluating the BLAST results
```
# Step 4: Inspect the results
head BLAST_results.tsv
# Displays first few hits: accession, description, alignment stats

# Step 5: Sort hits by sequence similarity
sort -k12,12nr BLAST_results.tsv | head
# Column 12 is the bit score; higher = better

# Step 6: Extract top 3 accession numbers
awk '{print $2}' BLAST_results.tsv | head -3
# Column 2 contains accession IDs
```


### BLAST spike glycoprotein excluding coronaviruses
```
# For coronaviruses ('Coronaviridae'), the NCBI taxid is 11118.

blastn \
  -query spike_glycoprotein.fasta \
  -db nt \
  -remote \
  -entrez_query "all[filter] NOT txid11118[Organism:exp]" \
  -outfmt '6 std stitle' \
  -out BLAST_non_coronaviruses.tsv
```
