### Introduction to SQLite for Biologists

#### Prerequisites
1. Download the annotations of the [SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) from NCBI in GFF3 format
2. Install the package `gffutils`to convert the GFF3 file to a database file and do the file conversion
```bash
pip install gffutils
gffutils-cli create NC_045512.2.gff3
```

#### Opening input files and start working with it
```bash
sqlite3 NC_045512.2.gff3.db
```

#### Exploring the dataset
```sqlite
--Printing the number of sequence features
SELECT COUNT(*) FROM features;

--Printing all features that are genes
SELECT * FROM features
WHERE featuretype = 'gene'

--Printing all features that are genes
SELECT * FROM features
WHERE featuretype = 'gene'
   OR featuretype = 'CDS'

--Printing all features that are NOT genes
SELECT * FROM features
WHERE featuretype != 'gene'

--Printing all features that are NEITHER genes NOR CDS
SELECT * FROM features
WHERE featuretype != 'gene'
  AND featuretype != 'CDS'
```
