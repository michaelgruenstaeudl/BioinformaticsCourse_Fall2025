### Using edirect to conduct data mining on GenBank

#### Conducting a simple search via esearch
```
# Defining a simple search query and saving it as a variable
MYQUERY="Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]
"

# Conducting the actual search via esearch
esearch -db nucleotide -query "$MYQUERY"

# Extracting the number of records found (i.e., count)
esearch -db nucleotide -query "$MYQUERY" | \
    xtract -pattern ENTREZ_DIRECT -element Count
```

#### Conducting a more complex search via esearch
```
# Defining a search query that also uses PROP (i.e., properties) 
# and FKEY (i.e., feature keys):

# Note: The query searches for any complete bacterial genome 
# published in 2023 that has a replication origin specified
MYQUERY="gbdiv_bct[PROP] \
         AND complete genome[TITLE] \
         AND 2023/01/01:2023/12/31[PDAT]\
         AND rep_origin[FKEY]"

# Summarizing the results via esummary:
esearch -db nucleotide -query "$MYQUERY" | esummary
```

#### Simplifying a complex search via efilter
```
# efilter allows the user to apply simpler, more general
# search queries and filter specific results after the initial esearch

# Note: The query searches for any complete bacterial genome 
# published in 2023 that has a replication origin specified

# WITHOUT efilter:
MYQUERY="gbdiv_bct[PROP] \
         AND complete genome[TITLE] \
         AND 2023/01/01:2023/12/31[PDAT]\
         AND rep_origin[FKEY]"
esearch -db nucleotide -query "$MYQUERY"

# WITH efilter (assuming a search on 01-Jan-2024):
MYQUERY="
         AND complete genome[TITLE] \
         
         AND rep_origin[FKEY]"
esearch -db nucleotide -query "$MYQUERY" | \
    efilter -division bct -days 365
```
