### Introducing grep
#### Text file filtering
```bash
grep    # Global search for Regular Expression and Print the results
        # i.e., printing lines that match patterns
```
## Simple operations with 'grep'

### Introduction to grep
```bash
# General usage
$ grep "PATTERN" <file_to_be_searched>

# Example usage with options
$ grep -i "PATTERN" <file_to_be_searched>

# To make out life easier, let us define the input file as a variable
$ INF=NC_045512.2.gbk
```

### Using grep on a GenBank formatted flatfile
```bash
$ cat $INF       # Display the entire file

$ grep "LOCUS" $INF              # print all lines containing 'LOCUS'

$ grep --color "DEFINITION" $INF # same, with keyword highlighted

$ grep -v "ACCESSION" $INF       # print lines without 'ACCESSION'

$ grep 'FEATURES\|ORIGIN' $INF   # lines with 'FEATURES' or 'ORIGIN'

$ grep "^ORIGIN" $INF            # lines starting with 'ORIGIN'
```

### Pipeling grep results
```bash
# lines with 'DEFINITION' but not 'ACCESSION'
$ grep -v "ACCESSION" $INF | grep "DEFINITION"

$ grep "PUBMED" $INF | wc -l	# count all lines with 'PUBMED'
```

### Grepping with context
```bash
$ grep -A2 "^DEFINITION" $INF	# the DEFINITION line + 2 following lines

$ grep -B3 "^FEATURES" $INF		# the FEATURES header + 3 preceding lines

$ grep -B1 -A2 "ORGANISM" $INF	# ORGANISM line + 1 before, 2 after
```

### Grepping with regular expressions
```bash
$ grep "gene.*ORF" $INF      	# full lines containing 'gene' and later 'ORF'

$ grep -o "/product=.*" $INF 	# print only everything after '/product='

$ grep "/product=" $INF | sort -u 	# print only unique '/product=' lines

$ grep "[3,5].UTR" $INF        	# any line containing 'gene' followed by numbers
```

### Exercise: Downloading additional GenBank records
```bash
# Download full genomes of 17 different SARS-CoV2 samples

$ for i in MT079851.1 MZ472096.1 OK439973.1 MZ353007.1 MW194121.1 \
MT412312.1 OU171384.2 OR477016.1 MZ544366.1 MW876953.1 OY715744.1 \
LR992043.1 OX648098.1 OU801385.1 OX446419.1 MZ579386.1 OK546254.1;
  do echo "Downloading $i"; 
  curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ \
  efetch.fcgi?db=nucleotide&id=${i}&rettype=gb&retmode=txt">$i.gbk;
done

# Note: No line breaks allowed in the curl command.
```

### Grepping across multiple files
```bash
# print all sequencing technologies used across the samples
$ grep "Sequencing Technology" *.gbk

# print all isolate information across the samples
$ grep "\isolate=" *.gbk

# print all cases where a journal with the title "Virology" has referenced that sequence
$ grep "JOURNAL .*Virology" *.gbk

# print all protein products that contain the keyword 'spike' across all samples
$ grep -o "/product=.*spike.*" *.gbk
```
