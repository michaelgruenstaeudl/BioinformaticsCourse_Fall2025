### Pattern recognition among local GenBank records using grep

### Downloading a single genome record from GenBank using the Bash shell
```
# Download the full reference genome of SARS-CoV-2 from GenBank

i=NC_045512.2
curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${i}&rettype=gb&retmode=txt">$i.gbk
```
### Downloading various additional genome records from GenBank using Bash shell
```
# Downloading the genome records of several additional SARS-CoV-2 amples from GenBank:

for i in MT079851.1 MZ472096.1 OK439973.1 MZ353007.1 MW194121.1 \
MT412312.1 OU171384.2 OR477016.1 MZ544366.1 MW876953.1 OY715744.1 \
LR992043.1 OX648098.1 OU801385.1 OX446419.1 MZ579386.1 OK546254.1;
  do curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${i}&rettype=gb&retmode=txt">$i.gbk;
done
```

### Using grep to extract relevant information from local GenBank records
```
# Grep can be executed on multiple files simultaneously

grep "Sequencing Technology" *.gbk   # print all lines with the key phrase
                                     # "Sequencing Technology" across all 
                                     # input files

grep "ACCESSION" *.gbk | wc -l  # count the number of correctly formatted 
                                # GenBank files you are working with


# Illustrating that not all of these SARS-CoV-2 genomes have the same 
# nucleotide sequence in the first 60 nucleotides

grep -h -A1 --no-group-separator "^ORIGIN" *.gbk | grep -v "ORIGIN"

```
