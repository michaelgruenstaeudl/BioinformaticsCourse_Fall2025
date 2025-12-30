### Introduction to Python for Biologists

#### Starting the Python IDE IDLE

```bash
# Start IDLE from the command line
idle3 &

# Alternative
python3 -m idlelib.idle
```

#### Basics: Variables and Strings

```python
# Assignment of string to variables
species = "Homo sapiens"
seq = "ATGCGTACGTTAG"

# Accessing string length
print(len(seq))

# String indexing/slicing
print(seq[0:3])   # First codon
print(seq[-3:])   # Last codon
```

#### Loops and Conditionals

```python
seq = "ATGCGTACGTTAG"

# Loop through seq
for base in seq:
    if base == "G":
        print("Found a G!")
```

#### Functions in Python

```python
def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq) * 100

dna = "ATGCGTACGTTAG"
print(gc_content(dna))
```


### First use of Biopython

#### Using Biopython: The Seq Object

```python
from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)

# Test the various inherited functions
my_seq.complement())
my_seq.reverse_complement()
my_seq.transcribe()
my_seq.translate()
```

#### Using Biopython: The Seq Object (functions of Seq)

```python
# ONLY FUNCTION NAMES
dir(my_seq)

# VERBOSE INFORMATION
help(my_seq)

# HELPFUL OVERVIEW
import inspect; from pprint import pprint
info = {
    name: inspect.getdoc(obj).splitlines()[0] if inspect.getdoc(obj) else ""
    for name, obj in inspect.getmembers(my_seq)
    if callable(obj) and not name.startswith("_")
}
pprint(info)
```


### An example with SARS-CoV-2

#### Reading FASTA Files

```python
# DETERMINING WORKING DIRECTORY
import os
print(os.getcwd())
os.chdir("/home/mi")

# IMPORTING FILE INTO PYTHON
from Bio import SeqIO
covid19_genome = SeqIO.read("NC_045512.fasta", "fasta")
covid19_genome.id
covid19_genome.seq
```


#### Calculating GC Content with Biopython

```python
from Bio.SeqUtils import gc_fraction

gc = gc_fraction(covid19_genome.seq) * 100

print(covid19_genome.id, "GC content:", round(gc, 2))
```


#### Back-translation of the Furin cleavage site

```python
# STEP 1. BACK-TRANSLATION
from Bio.Seq import Seq

site_motif_AA = Seq("PRRA")
codon_choice = {
    "P": "CCT",
    "R": "CGG",
    "A": "GCA",
}
site_motif_DNA = "".join(codon_choice[aa] for aa in site_motif_AA)
print(site_motif_DNA)
```


#### Test for presence of Furin cleavage site

```python
# STEP 2. TEST FOR PRESENCE
seq_str = str(covid19_genome.seq)

if site_motif_DNA in seq_str:
    print(covid19_genome.id, "contains motif", site_motif_DNA)
```


#### Test for location of Furin cleavage site

```python
# STEP 3. TEST FOR LOCATION
pos = covid19_genome.seq.find(site_motif_DNA)

if pos > 0:
    print("Furin cleavage site found at position:", pos)
else:
    print("Furin cleavage site not found.")
```
