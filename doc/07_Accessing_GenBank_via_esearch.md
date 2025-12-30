### Accessing sequence records via Entrez

#### Entrez search interface - Example query 1
```bash
Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]
```

#### Entrez search interface - Example query 2
```bash
complete genome[TITLE]
AND (chloroplast[TITLE] OR plastid[TITLE])
AND 2000/01/01:2020/12/31[PDAT]
AND 50000:250000[SLEN]
NOT unverified[TITLE]
NOT partial[TITLE]
AND Magnoliopsida[ORGN]
```

#### Installing Entrez Direct
```bash
# Linux or emulated Linux
apt install ncbi-entrez-direct

# MacOS
brew install ncbi-entrez-direct

# Otherwise, please follow the install instr. under:
# https://www.ncbi.nlm.nih.gov/books/NBK179288/#_chapter6_Getting_Started_
```

#### If you receive the error 'curl: (35) OpenSSL/3.0.17 ...' when operating edirect
```bash
# Create a new conda environment
conda create -n edirect_openssl1.1 perl openssl=1.1 curl
# Activate it
conda activate edirect_openssl1.1
# Install edirect manually inside this environment
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
# Add to PATH
export PATH=$HOME/edirect:$PATH
```

#### Defining the search query
```bash
MYQUERY="Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]
"
```

#### Conducting esearch and extracting XML output
```bash
esearch -db nucleotide -query "$MYQUERY"

esearch -db nucleotide -query "$MYQUERY" | \
xtract -pattern ENTREZ_DIRECT -element Count
```

#### Tips & Tricks
```bash
einfo -db nucleotide -fields

efetch -format gb -stop 1
```

#### Obtaining list of available feature keys
```bash
# Let us download the document and then grep it!
curl -s https://www.insdc.org/submitting-standards/feature-table/#7.2 > \
feature-table.info

# Get list of available feature keys:
grep "^Feature Key" feature-table.info

# Get more info on feature keys:
grep -v "<.*>" feature-table.info | grep -A6 "^Feature Key"
```

#### Esearch using indexed field 'FKEY'
```bash
# Searching for any complete bacterial genome published this year that
# has a replication origin specified:

MYQUERY="gbdiv_bct[PROP] \
AND complete genome[TITLE] \
AND 2023/01/01:2023/12/31[PDAT]\
AND rep_origin[FKEY]"

esearch -db nucleotide -query "$MYQUERY"

# Let us have a look at that single record:
esearch -db nucleotide -query "$MYQUERY" | esummary
```
