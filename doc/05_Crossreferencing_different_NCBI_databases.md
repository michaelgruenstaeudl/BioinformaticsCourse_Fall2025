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
