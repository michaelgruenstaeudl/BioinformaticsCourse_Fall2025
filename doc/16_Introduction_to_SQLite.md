### Introduction to SQLite for Biologists

#### Prerequisites
Preparing the necessary software tools and annotated sequence data for working with SQL. Specifically, we will work with the curated genome annotation for the SARS-CoV-2 genome.
1. Download the annotations of the [SARS-CoV-2 genome](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) from NCBI in GFF3 format
2. Install the package `gffutils`to convert the GFF3 file to a database file and do the file conversion
   
```bash
pip install gffutils
gffutils-cli create NC_045512.2.gff3
```

#### Opening input files and start working with it
Opening the newly created SQLite database that stores the genome annotations.

```bash
sqlite3 NC_045512.2.gff3.db
```

#### Exploring the dataset - Part 1
Introducing basic queries to inspect the structure and content of the annotation database

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

#### Exploring the dataset - Part 2
Introducing basic queries to count the number of all or specific features as well as group or order them

```sql
--Count the number of features and group by featuretype
SELECT featuretype, COUNT(*) AS n
FROM features
GROUP BY featuretype
ORDER BY n DESC;

--Count the number of specific features and group by featuretype
SELECT featuretype, COUNT(*) AS n
FROM features
WHERE featuretype = 'gene'
   OR featuretype = 'CDS'
GROUP BY featuretype
ORDER BY n DESC;

--Identify genes with more than one feature associated with it
SELECT
  json_extract(attributes, '$.gene[0]') AS gene,
  COUNT(*) AS total_features_in_table
FROM features
WHERE json_extract(attributes, '$.gene[0]') IS NOT NULL
GROUP BY gene
ORDER BY total_features_in_table DESC;

--For each gene, count the different feature classes and keep only those genes with at least one CDS
SELECT
  json_extract(attributes, '$.gene[0]') AS gene,
  COUNT(*) AS total_features,
  SUM(CASE WHEN featuretype = 'CDS'  THEN 1 ELSE 0 END) AS n_cds,
  SUM(CASE WHEN featuretype = 'mRNA' THEN 1 ELSE 0 END) AS n_mrna,
  SUM(CASE WHEN featuretype LIKE '%UTR%' THEN 1 ELSE 0 END) AS n_utr
FROM features
WHERE json_extract(attributes, '$.gene[0]') IS NOT NULL
GROUP BY gene
HAVING n_cds > 0
ORDER BY n_cds DESC, total_features DESC;
```

#### Performing simple selections
Filtering/selecting genomic features based on their genomic coordinates and biological meaning

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
```

#### Altering a table / adding new column
Extracting gene names from the annotation metadata and storing them in a dedicated column for easier access

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

#### Performing advanced selections
Filtering/selecting biologically meaningful information embedded within free-text annotations

```sql

--Selects all features whose attribute item 'note' contains the phrase “produced by” and extracts the subsequent text into a new column
/* Note:
substr(X, start pos, length): returns substring of X, continues until end if length not specified
instr(X, Y): returns position where string Y starts inside X
*/
SELECT
  featuretype,
  gene_name,
  substr(
    json_extract(attributes, '$.Note[0]'),
    instr(json_extract(attributes, '$.Note[0]'), 'produced by') + 12
  ) AS produced_by
FROM features
WHERE json_extract(attributes, '$.Note[0]') LIKE '%produced by%';
```

#### Joining features
Adding new columns from a subset of the original data

--For all features of type 'CDS', make a new column and add it ("join it") to the selected data under the column name 'synonym' as long as none of the attribute items are 'intron'
```sql
SELECT
  f.start,
  f.end,
  f.strand,
  json_extract(f.attributes, '$.gene[0]') AS gene,
  f.featuretype AS original,
  x.new_class   AS synonym
FROM features AS f
JOIN (
  SELECT 'CDS' AS original_class, 'exon' AS new_class
) AS x
ON f.featuretype = x.original_class
WHERE lower(f.attributes) NOT LIKE '%intron%';
```
