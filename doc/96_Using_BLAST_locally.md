### Using BLAST locally

#### Comparing nucleotide sequences directly
Assume that we wish to identify the largest nucleotide difference (in bp) between the human and the chimpanzee COI gene
```
# We download 10 human mitochondrial genomes and set up a local BLAST database of these genome records.
esearch -db nucleotide -query "complete mitochondrial genome[TITL] AND homo sapiens[ORGN]" | efetch -stop 10 -format fasta > human_mitos.fasta
makeblastdb -in human_mitos.fasta -title 10_human_mitos -dbtype nucl -out human_mitos_db

# We download a random copy of a chimpanzee mitochondrial COI gene.
esearch -db nucleotide -query "COI[TITL] AND pan troglodytes[ORGN]" | efetch -stop 1 -format fasta > chimpanzee_COI.fasta

# We BLAST the nucleotide sequence of the chimpanzee mitochondrial COI gene against the local database.
blastn -query chimpanzee_COI.fasta -db human_mitos_db > blast_out.txt
cat blast_out.txt | grep "Identities"  # largest nucleotide difference encountered: 59 bp
```
