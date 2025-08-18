### Assembly of plant plastid genome

#### Step 1. Extracting all plastome reads from a genome skimming read set
Using blastn against a local database of existing plastid genomes to extract those reads from the genome skimming read set that represent the plastid genome

#### Step 2. Assembly with different plastome assemblers
Different assemblers generate slightly different assembly sequences
```
# Plastid genome assembly with NOVOPlasty
perl NOVOPlasty.pl -c config.txt

# Plastid genome assembly with GetOrganelle
get_organelle_from_reads.py -1 mapped_F12_R1.fastq -2 mapped_F12_R2.fastq -t 1 -o mapped_F12G.plastome -F embplant_pt -R 10

# Plastid genome assembly with SPAdes
spades.py --isolate -1 mappedF12_reads_R1.fastq -2 mappedF12_reads_R2.fastq -o spades_sbatch -t 2
```

#### Step 3. Comparison of assembly results via blastn
The plastome assembly via NOVOPlasty is designated the reference assembly and used to generate a local BLAST database
```
# Make local BLAST database
makeblastdb -in mapped_F12.fasta -title mapped_NOVO -dbtype nucl -out mapped_NOVO_db

# Compare other assembly results to the NOVOPlasty assembly via blastn
blastn -query mapped_F12G.fasta -db mapped_NOVO_db > blast_out.txt
cat blast_out.txt | grep "Identities" | sed 's/Identities =//g' | sort -n |tail -n 4
```
