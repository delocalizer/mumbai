# Goal

Faster position-based slice or lookup on multiple bams.

# Background

You have dozens or more indexed bams sitting on fast disk and you wish
to query them all by position. Possible solutions include:

1. Query all the bams in parallel and aggregate.
2. Create an aggregate data structure out of all the bams (e.g. in HDF5) and
   query that instead.

They both have some appeal — the first needs no pre-processing time or extra
storage; the second requires less compute at runtime and should be faster. Here
we try to improve upon the most naïve implementation of the first approach by
creating an aggregate _index_ from the .bai files, so that instead of reading
all indices to locate all the offsets we make one query to the aggregate index.

# Implementation

1. A database of bam index ('bai') files that can be queried by genomic coordinates
   and returns virtual file offsets to the indexed bams. This allows lookup on hundreds
   or thousands of bams without first reading all their indexes. 

2. A client that queries the database by genomic coordinates and returns records
   that overlap.

# Install

```bash
git clone https://github.com/delocalizer/mumbai
cd mumbai
pip install -r requirements.txt --user .
```

# Test
(optional, requires tox)
```bash
tox
```

# Usage

0. get help

```bash
mumbai_db --help  # for help on the db create and load tool
mumbai --help     # for help on the db client
```

1. create a bam index database for bams aligned to GRCh38

```bash
mumbai_db create newdb GRCh38
```

2. load the database with bam indexes

```bash
mumbai_db load newdb /path/to/first.bam /path/to/second.bam ...
```

3. query the bams by position and return overlapping records in SAM format

```bash
mumbai sam newdb chr1 1000000 1000100
```

4. query the bams by position and count overlapping records

```bash
mumbai count newdb chr1 1000000 1000100
```

5. query the bams by position and visualize the overlapping region

```bash
mumbai tview newdb chr1 1000000 1000100
```

5. query the bams by position and visualize the overlapping region as a pileup

```bash
mumbai pileup newdb chr1 1000000 1000100
```
