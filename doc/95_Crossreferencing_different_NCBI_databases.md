### Cross-linking different NBCI databases

#### Finding links between databases with 'elinks'
```
# Generating (abbreviated) UID list
MYQUERY="(mitochondrion[TITLE] OR mitochondrial[TITLE]) \
    AND complete genome[TITLE] \
    AND (human[TITLE] or homo sapiens[TITLE])"
    
esearch -db nuccore -query "$MYQUERY" | efetch -format uid | \
    head -n1000 > human_mitos_uids.list

# Confirming that UIDs actually our target genomes
cat human_mitos_uids.list | epost -db nuccore | \
    esummary | xtract -pattern DocumentSummary -element Title

# Finding actual links (this step takes 15-20 min.)
cat human_mitos_uids.list | epost -db nuccore | \
    elink -target biosample | \
    xtract -pattern ENTREZ_DIRECT -element Count
```

#### Finding human mitochondrial genomes crosslinked between NCBI Nucleotide and NCBI SRA
Are there complete human mitochondrial genomes that are BOTH present as a regular GenBank record on NCBI Nucleotide, and present as that GenBank record's short sequence reads on NCBI SRA.
```
# SOLUTION 1 (no hits found): NCBI SRA --> NCBU Nucleotide
# Step 1:
esearch -db sra -query "genome[TITLE] AND (human[ORGN] OR Homo sapiens[ORGN]) AND (mitochondrion[TITLE] OR mitochondrial[TITLE] OR mito[TITLE])" | esummary | xtract -pattern DocumentSummary -element Biosample > bios_biop_links.txt &
# Step 2:
esearch -db sra -query "genome[TITLE] AND (human[ORGN] OR Homo sapiens[ORGN]) AND (mitochondrion[TITLE] OR mitochondrial[TITLE] OR mito[TITLE])" | esummary | xtract -pattern DocumentSummary -element Bioproject >> bios_biop_links.txt &
# Step 2b 
cat bios_biop_links.txt | sort | uniq > bios_biop_links_uniq.txt
# Step 3:
for i in $(cat bios_biop_links_uniq.txt); do echo $i >> almost_final_output.txt; esearch -db nucleotide -query "$i" >> almost_final_output.txt; done &
# Step 4
cat almost_final_output.txt | grep  "Count\|PRJ"
cat almost_final_output.txt | grep -B5 "<Count>1" | grep "PRJNA" (needs minor modified)
# Step 5
for i in $(cat almost_final_output.txt); do echo $i; esearch -db nucleotide -query "$i" | esummary >> almost_final_output.txt; done 
```

```
# SOLUTION 2 (2-5 hits found): NCBU Nucleotide --> NCBI SRA
# Initial esearch that downloads accession numbers
esearch -db nucleotide -query "(human[ORGN] OR homo sapiens[ORGN]) AND (mitochondria[TITLE] OR mitochondrial[TITLE] OR mitochondrion[TITLE]) AND complete genome[TITLE]" | efetch -format acc > accession_numbers.txt
# Takes accession numbers and individually searches through and greps the SAMN and/or PRJNA numbers
for i in $(cat accession_numbers.txt); do 
	echo $i; 
	echo $i >> PRJNA_SAMN_nucleotide.txt; 
	efetch -db nucleotide -id $i -format gb | grep "SAMN\|PRJNA" >> PRJNA_SAMN_nucleotide.txt; 
done
# Grepping so that only the PRJNA and SAMN numbers go into the new file
# I had to manually go in and remove what was left of the information since I used .*
grep -o "PRJNA.*\|SAMN.*" PRJNA_SAMN_nucleotide.txt > PRJNA_SAMN.txt
# To show link between nucleotide/nuccore and sra
for i in $(cat PRJNA_SAMN.txt); do 
	esearch -db sra -query "$i"; 
done
# If you want to efetch the files... this might be a lot of files considering some of the project numbers have multiple hits associated with them
for i in $(cat PRJNA_SAMN.txt); do 
	esearch -db sra -query "$i" | efetch;
done
```

```
# SOLUTION 3 (2-5 hits found): NCBU Nucleotide --> NCBI SRA
MYQUERY="complete genome[TITLE] AND (human[ORGN] OR Homo sapiens[ORGN]) AND (mitochondrial[TITLE] OR mitochondrion[TITLE] OR mitochondria[TITLE])"
esearch -db nucleotide -query "$MYQUERY" | efetch -format acc > accession_num.txt
nohup sh -c ' for i in $(cat accession_num.txt); do echo $i >> dirtyDB.txt; efetch -db nucleotide -id $i -format gb | grep "SAMN\|PRJNA" >> dirtyDB.txt; done' &
grep -o "SAMN........" dirtyDB.txt > cleanDB.txt
grep -o "PRJNA.*" dirtyDB.txt >> cleanDB.txt
esearch -db sra -query "PRJNA422662" | esummary | xtract -pattern DocumentSummary -element Title
esearch -db sra -query "PRJNA422662" | esummary | xtract -pattern DocumentSummary -element Run@acc
```

```
# OPTIMAL SOLUTION (5-10 hits found): NCBU Nucleotide --> NCBI SRA
esearch -db nucleotide -query "(human[ORGN] OR homo sapiens[ORGN]) AND (mitochondrial[TITLE] OR mitochondrion[TITLE]) AND complete genome[TITLE]" | efetch -format acc > accession_numbers.txt &

nohup sh -c 'for i in $(cat accession_numbers.txt); do
  echo $i >> simple_search_output.txt;
  efetch -db nucleotide -id $i -format gb | grep "SAMN\|PRJNA" >> simple_search_output.txt;
  echo "" >> simple_search_output.txt;
done ' &

grep -B1 -A1 "DBLINK" PRJNA_SAMN_nucleotide.txt

esearch -db sra -query "PRJNA422662" | esummary | xtract -pattern DocumentSummary -element Run@acc

% Example Output:
% MG936619.1
% DBLINK      BioProject: PRJNA422662
%             BioSample: SAMN08193509
% SRR6664768 .. SRR6664785
```

