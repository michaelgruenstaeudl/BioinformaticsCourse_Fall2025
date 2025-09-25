#### Tools that make your life easier
```
# Making life easier with 'esummary'
# 'esummary' is the same as 'efetch -format docsum'

esearch -db nucleotide -query "$MYQUERY" | efetch -format docsum
esearch -db nucleotide -query "$MYQUERY" | esummary
```

```
# Making life easier with 'efilter'
# 'efilter' provides additional filtering options (beyond what the query allows)

esearch -db nucleotide -query "$MYQUERY" |
efilter -division bct -days 260
```

```
# Making life easier with 'xtract'
# Use xtract to parse out specific information from the XML output

esearch -db nucleotide -query "$MYQUERY" |
esummary |
xtract -pattern DocumentSummary -element Title
```



#### Extracting and reusing id numbers
```
esearch -db nucleotide -query "$MYQUERY" |
efetch -format uid > uid_numbers.list

# Note: uid number is different from accession number 
cat uid_numbers.list

# Using the accession numbers again
cat uid_numbers.list |
epost -db nucleotide | esummary |
xtract -pattern DocumentSummary -element Title
```


#### Extracting taxonomic information
```
# You can query NCBI Taxonomy just like any other NCBI database

esearch -db taxonomy -query "Limnothrix" |
efetch -format xml |
xtract -pattern Taxon -element Lineage
```



#### How to known which indexed field are available?
```
# Use command:
einfo -db database_name -fields

# Example: Indexed fields available for database 'nuccore'
einfo -db nuccore -fields

# Example: Indexed fields available for database 'taxonomy'
einfo -db taxonomy -fields
```


#### Fundamental questions
```
# What is the the longest and the shortest complete bacterial genome sequence ever published on GenBank? (Note: Test initially with 100 records.)

MYQUERY="gbdiv_bct[PROP] AND complete genome[TITLE]"
esearch -db nucleotide -query "$MYQUERY" |
esummary -stop 100 |
xtract -pattern DocumentSummary -element Slen -element Title |
sort -k1,1n
```


```
# What is the most common host organism of all complete viral genome sequences ever published on GenBank? (Note: Test initially with 100 records.)

MYQUERY="gbdiv_vrl[PROP] AND complete genome[TITLE] AND source[FKEY]"
esearch -db nucleotide -query "$MYQUERY" |
efetch -format gb -stop 100 |
 grep "/host=\"" | sort |  uniq -c | sort -nr
```



### Finding similar entries in other NCBI databases

#### Which databases are there and how are they linked?
```
# There are various databases accessible via edirect
einfo -dbs

# Each database has numerous links to other databases (at least in theory!)
einfo -db nuccore -links
einfo -db bioproject -links
einfo -db biosample -links
einfo -db sra -links
```


#### Finding SRA entries with connections to nuccore
```
# Aim: Finding short-read datasets stored on NCBI SRA that match genome sequences on NCBI Nucleotide.}
# Critical aim for re-assembly of genomes!

esearch -db sra -query "Limnothrix[ORGN]" | 
elink -target nuccore

esearch -db sra -query "Pseudanabaenales[ORGN]" | 
elink -target nuccore

# Does that mean that no short-read dataset stored on NCBI SRA that matches genome sequences on NCBI GenBank for the genus \textit{Limnothrix}?
```


#### Single elink connections do not produce hits
```
# Evidence of links between 'biosample' and 'nuccore' for \textit{Limnothrix} cyanobacteria
esearch -db biosample -query "Limnothrix[ORGN]" |
elink -target nuccore | esummary | 
xtract -pattern DocumentSummary -element Title

# Evidence of links between 'biosample' and 'SRA' for \textit{Limnothrix} cyanobacteria
esearch -db biosample -query "Limnothrix[ORGN]" |
elink -target sra | esummary | 
xtract -pattern DocumentSummary -element Title
```


#### How to find SRA datasets that match GenBank genomes?
```
# Possible solution: from nuccore via bioproject to SRA

esearch -db nuccore -query "Limnothrix[ORGN] AND 100000:10000000[SLEN]" | 
elink -target bioproject | 
elink -target sra | 
esummary |
xtract -pattern DocumentSummary -element Title -ACC @acc -block DocumentSummary -element "&ACC"
```


```
# Possible solution--inverse: from SRA via bioproject to nuccore

esearch -db sra -query "Limnothrix[ORGN] AND (Nanopore [TITLE] OR Illumina [TITLE] OR WGS[TITLE] OR complete genome[TITLE])" | 
elink -target bioproject | 
elink -target nuccore | 
esummary |
xtract -pattern DocumentSummary -element Title -ACC @acc -block DocumentSummary -element "&ACC"

# IMPORTANT: Any NCBI nuccore record now needs to be matched with the correct short-read dataset on NCBI SRA, for example by comparing the metadata (e.g., strain name, lab name, country name, publication title, etc.).
```