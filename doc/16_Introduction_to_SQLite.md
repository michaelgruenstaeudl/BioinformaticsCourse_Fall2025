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
```sql
--Print the attributes that the table 'features' contains
PRAGMA table_info(features);

--Print the number of features in table 'features'
SELECT COUNT(*) FROM features;

--Print the first 10 features
SELECT * FROM features
LIMIT 10;

--Print all features that are genes
SELECT * FROM features
WHERE featuretype = 'gene';

--Print all features that are genes
SELECT * FROM features
WHERE featuretype = 'gene'
   OR featuretype = 'CDS';

--Print all features that are NOT genes
SELECT * FROM features
WHERE featuretype != 'gene';

--Print all features that are NEITHER genes NOR CDS
SELECT * FROM features
WHERE featuretype != 'gene'
  AND featuretype != 'CDS';
```

#### Performing specific selections
```sql
--Select all genes of a specific genome region (i.e., the region between position 10,000 and 30,000)
SELECT featuretype, start, end, attributes FROM features
WHERE featuretype = 'gene'
  AND start >= 10000
  AND end   <= 30000
ORDER BY start;

--Extract all gene names
/* Note:
Data in the column 'attributes' are stored in JSON-format: {"ID":["gene-GU280_gp01"],"Dbxref":["GeneID:43740578"],"gbkey":["Gene"],"gene":["ORF1ab"],...}
*/
SELECT
  json_extract(attributes, '$.gene[0]') AS gene_name
FROM features
WHERE featuretype = 'gene';

--Extract substrings following keyword 'a_string='
/* Note:
instr(X, Y): returns position where string Y starts inside X
substr(X, start, length): returns substring of X
*/
SELECT
  featuretype,
  start,
  end,
  substr(
    attributes,
    instr(attributes, 'a_string=') + 8,
    instr(substr(attributes, instr(attributes, 'a_string=') + 8), ';') - 1
  ) AS my_substring
FROM features
WHERE attributes LIKE '%a_string=%'
  AND featuretype = 'gene';
```

#### Altering the table
```sql
--Extract all gene names from attributes and add them as a new column

--- Step 1. Define new (but empty) column
ALTER TABLE features
ADD COLUMN gene_name TEXT;

--- Step 2. Populate new column with gene names
/* Note:
Command 'SET' sets NULL where no 'gene' key exists
*/
UPDATE features
SET gene_name = json_extract(attributes, '$.gene[0]');

--- Step 3. Verify successful extraction
SELECT featuretype, gene_name
FROM features
WHERE gene_name IS NOT NULL
```
