### Using edirect to conduct data mining on GenBank

#### Conducting a simple search via esearch
```
# Defining the search query and saving it as a variable
MYQUERY="Severe acute respiratory syndrome coronavirus 2[TITLE] \
AND complete genome[TITLE] \
AND 29500:30500[SLEN] \
AND 2020/01/01:2020/06/31[PDAT]
"

# Conducting the actual search via esearch
esearch -db nucleotide -query "$MYQUERY"
```
