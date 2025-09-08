## Grep with approximate matching

### Introduction to agrep
```
# Approximate text file filtering
agrep    # like grep, just with the ability to conduct
         # approximate matching

# General usage
agrep -# "PATTERN" <file_to_be_filtered>  # with '#' being the number 
                                          # of maximal differences allowed
```

### Using agrep
```
# Extracting all matches that contain the keyword "SARS-CoV-2" exactly
agrep -0 "SARS-CoV-2" *.gbk   # is equal to: grep "SARS-CoV-2" *.gbk

# Extracting all matches that contain the keyword "SARS-CoV-2" in an 
# altered form, namely with up to three letters being different
agrep -3 "SARS-CoV-2" *.gbk | grep -v "SARS-CoV-2"

# Extracting all matches that contain the keyword "SARS-CoV-2" in an 
# altered form, namely with up to five letters being different
agrep -5 "SARS-CoV-2" *.gbk | grep -v "SARS-CoV-2"
```

### Producing our first realistic genomic dataset
```
# I am searching for *all* SARS-CoV-2 genomes that were submitted
# to GenBank in the first half of 2020
esearch -db nucleotide -query \
"Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]"
```

### If you have difficulties with edirect in 2025:
```
# Create a new conda environment
conda create -n edirect_openssl1.1 perl openssl=1.1 curl
# Activate it
conda activate edirect_openssl1.1
# Install edirect manually inside this environment
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
# Add to PATH
export PATH=$HOME/edirect:$PATH

esearch -db nucleotide -query \
"Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]"
```

### Downloading a LARGE genomic dataset
```
# Downloading all 5809 genome files and combining them into one file
esearch -db nucleotide -query \
"Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]
"  | efetch -format gb >> SARS-CoV-2_collected_genomes_2020H2.gbk &

# Compressing the file to reduce its file size
gzip SARS-CoV-2_collected_genomes_2020H2.gbk
```

### Confirming that full file contains correct number of records
```bash
# Let us verify that all 5,809 SARS-CoV-2 genomes are saved in the 
# compressed file 'SARS-CoV-2_collected_genomes_2020H2.gbk.gz'

# Checking for beginning of each file
zgrep "^LOCUS" SARS-CoV-2_collected_genomes_2020H2.gbk.gz | wc -l

# Checking for end of each file
zgrep "^//" SARS-CoV-2_collected_genomes_2020H2.gbk.gz | wc -l
```

### Conducting data mining on complete dataset
```
# Using grep to summarize what kind of sequencing technologies 
# were used to produce these genome records
zgrep -h "Sequencing Technology" SARS-CoV-2_collected_genomes_2020H2.gbk.gz | \
sort | uniq -c | sort -n

# How many different spelling versions of 'Illumina' are there among
# these 5,809 files?
zcat SARS-CoV-2_collected_genomes_2020H2.gbk.gz | agrep -h -2 | grep -o "Il......" | \
sort | uniq -c

# What is the distribution of countries that these genomes come from?
zgrep -h "geo_loc_name=" SARS-CoV-2_collected_genomes_2020H2.gbk.gz | grep -o '".*"' | \
tr -d '"' | sort | uniq -c | sort -n
```
