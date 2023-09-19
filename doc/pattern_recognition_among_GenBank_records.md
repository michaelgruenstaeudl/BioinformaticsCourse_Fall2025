### Pattern recognition among local GenBank records using grep

### Downloading a single genome record from GenBank using the Bash shell
```
# Download the full reference genome of SARS-CoV-2 from GenBank

i=NC_045512.2
curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${i}&rettype=gb&retmode=txt">$i.gbk
```


