### Introduction to SQLite for Biologists

#### Prerequisites
1. Download the annotations of the [SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) from NCBI in GFF3 format
2. Install the package `gffutils`to convert the GFF3 file to a database file and do the file conversion
```
pip install gffutils
gffutils-cli create NC_045512.2.gff3
```

#### Opening input files and start working with it
```
sqlite3 NC_045512.2.gff3.db

.tables
SELECT COUNT(*) FROM features;
```
