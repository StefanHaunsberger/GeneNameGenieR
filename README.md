

# _GeneNameGenieR_

The _GeneNameGenieR_ R package has been developed to provide functions for molecular identifier conversion, such as gene symbols, aliases, transcripts and miRNAs. The base thereby is a Neo4j graph database containing identifiers from Ensembl v91 and 22 miRBase release versions ([miRBase](http://www.mirbase.org)).

* [Introduction](#introduction)
* [Use Cases](#use-cases)
   - [How to translate molecular identifiers](#how-to-translate-molecular-identifiers)
       - [Get official gene symbol](#get-official-gene-symbol)
       - [Convert identifiers from and to different databases](#convert-identifiers-from-and-to-different-databases)
   - [How to translate miRNAs](#how-to-translate-mirnas)
       - [Translate input Ids to current version and retrieve metadata](#translate-input-ids-to-current-version-and-retrieve-metadata)
       - [Translate mature miRNA names to different miRBase release versions](#translate-mature-mirna-names-to-different-mirbase-release-versions)
   - [List valid databases and supported parameter values](#list-valid-databases-and-supported-parameter-values)
   - [Performing Reactome Pathway Enrichment Analysis](#performing-reactome-pathway-enrichment-analysis)
   - [Additional information](#additional-information)

*Availability:* 

The package can be installed using the following R command from the [devtools](https://github.com/r-lib/devtools) package:

```r
devtools::install_github("stefanhaunsberger/GeneNameGenieR")
```

# Introduction

The _GeneNameGenieR_ package provides a programmatic interface to the GeneNameGenie
Neo4j graph database (GDB). The GeneNameGenie Neo4j GDB can either be used as standalone
or also as part of miRGIK, in which GeneNameGenieR is fully integrated.
GeneNameGenieR provides functions for translating any general molecular identifier which 
is supported by Ensembl v91. Functions include translation to the official gene symbol, 
retrieving metadata for Ensembl genes and transcripts and translating between different
database formats. Thereby, providing the source database of the input identifier is optional.

Moreover GeneNameGenieR supports 22 different miRBase release versions, with the most recent one 
from March 2018, miRBase version 22. miRNAs can be translated to the most recent as well as
any other supported miRBase version version and metadata, such as sequence, previous 
identifiers or type, can also be included.

Reactome pathway enrichment analysis can be performed on Reactome v.63.

The main methods are: 

- `getOfficialGeneSymbol`: Convert input id to its corresponding official gene symbol,
- `convertFromTo`: Convert molecular identifier to target identifiers,
- `getValidDatabases`: List all supported databases
- `convertToCurrentMirbaseVersion`: Get current miRNA name for mature input miRNAs
- `convertMatureMirnasToVersions`: Translate mature miRNA names to different versions.
- `getValidMirnaMetadataValues`: List all valid miRNA metadata parameters.
- `rpea`: Reactome pathway enrichment analysis.

 
To load the package and gain access to the functions just run the 
following command:

```r
library(GeneNameGenieR)
```

## Database information

GeneNameGenieR depends on a GeneNameGenie Neo4j graph database instance. This can either 
be locally or online. The default URL is set to `http://localhost:7474/db/data/`.
Hence, the default values , `host = http://localhost`, `port = 7474` and `path = db/data/`.

If a different host address needs to be set, the `setNeo4jConnection` function can be used:
```r
# setNeo4jConnection(host = "<url>", port = <port>, path = "<path>")
```

# Use Cases

## How to translate molecular identifiers

### Get official gene symbol

To retrieve the official gene symbol for any given molecular identifier, the 
`getOfficialGeneSymbol` function can be used.
Let us assume we want to find the official gene symbol for the following molecular 
identifiers: 'AMPK', 'Bcl-2', '596' and 'NM_000657'. We can use the `getOfficialGeneSymbol` 
function to retrieve the official gene symbol for each input identifier respectively.
```r
ids = c('P10415', 'AMPK', 'ENSP00000381185', 'ENST00000589955', 'A0A024R2B3', '596');
getOfficialGeneSymbol(ids);
```
The output is a data.frame object. Thereby, if warnings appear on the console, affected rows are printed at the top of the results. 
The actual result is the dataframe containing the three columns `InputId`, `InputSourceDb` and `OfficialGeneSymbol`.

```r
# A tibble: 2 x 2
     InputId     n
       <chr> <int>
1        596     2
2 A0A024R2B3     2
3     P10415     2
# --------
          InputId             InputSourceDb OfficialGeneSymbol
1 ENSP00000381185 Ensembl Human Translation               BCL2
2 ENST00000589955  Ensembl Human Transcript               BCL2
3             596                  WikiGene               BCL2
4             596                 NCBI gene               BCL2
5          P10415       UniProtKB Gene Name               BCL2
6          P10415      UniProtKB/Swiss-Prot               BCL2
7      A0A024R2B3       UniProtKB Gene Name               BCL2
8      A0A024R2B3          UniProtKB/TrEMBL               BCL2
9            AMPK         Gene Symbol Alias             PRKAA2
Warning message:
In postCheckInput(x) :
  Some input identifiers match to more than one input source database
```
Because we have not provided the source database parameter we get a warning that some input 
identifiers were actually mapped to multiple databases, such as P10415, which is a UniProtKB 
gene name as well as a UniProtKB/Swiss-Prot protein name.
We can run the same command again with the `Uniprot_gn`, for the UniProtKB Gene Name database, 
as input source database parameter. 

```r
ids = c('P10415', 'AMPK', 'ENSP00000381185', 'ENST00000589955', 'A0A024R2B3', '596');
getOfficialGeneSymbol(ids, sourceDb = 'Uniprot_gn');
```
Outputs:

```r
     InputId       InputSourceDb OfficialGeneSymbol
1     P10415 UniProtKB Gene Name               BCL2
2 A0A024R2B3 UniProtKB Gene Name               BCL2
```
This time we only received the official gene symbol for input identifiers which are from the 
UniProtKB gene name database.

The function `getValidDatabases` can be used to look up valid `sourceDb` values.

### Convert identifiers from and to different databases

With the `convertFromTo` function we can convert any given identifier to any supported 
target database. In the following example we are want to get all identifiers from 
`"Uniprot/SWISSPROT"`, `"Uniprot/SPTREMBL"`, `"EntrezGene"` for the `"BCL2"` and `"AMPK"`.

```r
convertFromTo(c("BCL2", "AMPK"), c("Uniprot/SWISSPROT", "Uniprot/SPTREMBL", "EntrezGene"));
```
Outputs:

```r
  InputId        InputSourceDb             TargetDb   TargetId
1    BCL2 Official Gene Symbol     UniProtKB/TrEMBL A0A024R2C4
2    BCL2 Official Gene Symbol UniProtKB/Swiss-Prot     P10415
3    BCL2 Official Gene Symbol     UniProtKB/TrEMBL A0A024R2B3
4    BCL2 Official Gene Symbol            NCBI gene        596
5    AMPK    Gene Symbol Alias     UniProtKB/TrEMBL A0A087WXX9
6    AMPK    Gene Symbol Alias UniProtKB/Swiss-Prot     P54646
7    AMPK    Gene Symbol Alias            NCBI gene       5563
```

## How to translate miRNAs

### Translate input Ids to current version and retrieve metadata

Most of the time all one needs is to retrieve the miRNA name of the current 
miRBase release version. The `convertToCurrentMirbaseVersion` function 
takes a single or multiple input values and returns the current miRNA name for 
each input value respectively (where possible). The following example illustrates 
a simple use case where a single miRNA is converted to the current miRBase version.

```r
> convertToCurrentMirbaseVersion('hsa-miR-29a');

      InputId    Accession   CurrentMirna CurrentVersion
1 hsa-miR-29a MIMAT0000086 hsa-miR-29a-3p             22
```

Moreover, metadata for mature as well as precursor miRNAs, such as comments and reads, 
can be retrieved by passing the parameter value `metadata` to the function. 
The function accepts mature miRNA identifiers (miRNA name and MIMAT-accession) as well as 
precursor identifiers (pre-miRNA name and MI-accession). Hence, some metadata values 
only apply to precursor miRNAs, such as genomic location and strand.

The following example returns the current miRNA name and metadata for a the mature 
miRNA name `hsa-miR-29a` and the precursor miRNA name `hsa-mir-29a`.

```r
> convertToCurrentMirbaseVersion(c('hsa-miR-29a', 'hsa-mir-29a'), metadata = c('nExperiments', 'strand', 'reads'));
      InputId    Accession   CurrentMirna CurrentVersion nExperiments Strand   Reads
1 hsa-miR-29a MIMAT0000086 hsa-miR-29a-3p             22           77   <NA> 2045292
2 hsa-mir-29a    MI0000087    hsa-mir-29a             22           77      - 2055333
```
As we can see, the strand information is only available for the precursor miRNA. 
It is to mention that not all metadata information is available for each miRNA. So, for 
example, the mature miRNA-X could return a value for the `communityAnnotation` metadata value 
whereas miRNA-Y does not just because one has a value in the database and the other one does not.

*Some more examples:*

```r
> convertToCurrentMirbaseVersion('hsa-mir-29a', metadata = c('sequence', 'type', 'previousIds'));
      InputId Accession CurrentMirna CurrentVersion                      Sequence      Type PreviousIds
1 hsa-mir-29a MI0000087  hsa-mir-29a             22 AUGACUGAUUUCUU...GAAAUCGGUUAU antisense  hsa-mir-29
```

```r
convertToCurrentMirbaseVersion(c('hsa-miR-29a', 'MI0000087'), metadata = c('nExperiments', 'evidenceType', 'reads'));
      InputId    Accession   CurrentMirna CurrentVersion nExperiments EvidenceType   Reads
1   MI0000087    MI0000087    hsa-mir-29a             22           77         <NA> 2055333
2 hsa-miR-29a MIMAT0000086 hsa-miR-29a-3p             22           77 experimental 2045292
```

### Translate mature miRNA names to different miRBase release versions

In cases where we want to use tools which require the miRNA from a certain miRBase release 
version one can use the `convertMatureMirnasToVersions` function.

If no target versions are provided, the name from the most recent supported miRBase 
release version is returned.

```r
> convertMatureMirnasToVersions('hsa-miR-29a');
         InputId MatureAccession miRBaseVersion    TargetMirna
   1 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p
```

The following code returns the names for miRBase version 17, 21 and 22 for the 
mature miRNA `hsa-miR-29a`.

```r
> convertMatureMirnasToVersions('hsa-miR-29a', c(17, 21, 22))
        InputId MatureAccession miRBaseVersion    TargetMirna
  1 hsa-miR-29a    MIMAT0000086             17    hsa-miR-29a
  2 hsa-miR-29a    MIMAT0000086             21 hsa-miR-29a-3p
  3 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p
```

In the case where also the sequence for a specific version is required we can 
pass the parameter `sequence = TRUE` to the function.

```r
convertMatureMirnasToVersions('hsa-miR-29a', c(17, 21, 22), sequence = TRUE)
       InputId MatureAccession miRBaseVersion    TargetMirna         TargetSequence
 1 hsa-miR-29a    MIMAT0000086             17    hsa-miR-29a UAGCACCAUCUGAAAUCGGUUA
 2 hsa-miR-29a    MIMAT0000086             21 hsa-miR-29a-3p UAGCACCAUCUGAAAUCGGUUA
 3 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p UAGCACCAUCUGAAAUCGGUUA
```

## List valid databases and supported parameter values

### `getValidDatabases`: for supported database values

The `getValidDatabases` function returns a list of supported database values. The values 
contained in the `DatabaseId` column are valid as `targetDb` and `sourceDb` parameter values.

```r
getValidDatabases();
```

Using the `head` function we can show only the first couple of entries:
```r
head(getValidDatabases());
                   DatabaseDisplayName                     DatabaseId
1                                 CCDS                           CCDS
2                               ChEMBL                         ChEMBL
3           Clone-based (Ensembl) gene       Clone_based_ensembl_gene
4     Clone-based (Ensembl) transcript Clone_based_ensembl_transcript
5 DataBase of Aberrant 3' Splice Sites                         DBASS3
6 DataBase of Aberrant 5' Splice Sites                         DBASS5
```

### `getValidMirnaMetadataValues`: List valid miRNA metadata parameter values

The `getValidMirnaMetadataValues` function returns a list of supported miRNA metadata parameter values. 
Some of the parameters apply only to mature miRNAs whereas others only return values for precursor miRNAs.

```r
getValidMirnaMetadataValues()
          ParameterKey                                     Description
1           confidence                           Annotation confidence
2                 type                                   Type of miRNA
3             sequence                             Nucleobase sequence
4             comments     Comments, such as general miRNA information
5          previousIds                            Previous miRBase IDs
6                  url                            URL to miRBase entry
7           chromosome                                   Chromosome ID
8          regionStart                    Location start on chromosome
9            regionEnd                      Location end on chromosome
10              strand                                    Strand (+/-)
11 communityAnnotation     Annotations and comments from the community
12        nExperiments                Number of validation experiments
13               reads Number of reads from the validation experiments
14        evidenceType          Type of evidence, such as experimental
```

### `getCurrentMirbaseVersion`: Get information on the latest miRBase release version supported by the package

```r
getCurrentMirbaseVersion()
[1] "miRBase Release 22, March 12, 2018"
```

## Performing Reactome Pathway Enrichment Analysis

```r
GeneNameGenieR::rpea(c('PRKAA1', 'PRKAB2', 'PRKAG2', 'STK11'), pValueCutoff = 0.00001)
     reactomeID                                         description   bgRatio geneRatio                    geneIds       pValue pValueAdjusted
1  R-HSA-380972    Energy dependent regulation of mTOR by LKB1-AMPK  29/10592       4/4 PRKAA1,PRKAB2,PRKAG2,STK11 4.531342e-11   6.797013e-10
2  R-HSA-165159                                     mTOR signalling  40/10592       4/4 PRKAA1,PRKAB2,PRKAG2,STK11 1.743590e-10   2.615385e-09
3 R-HSA-6804756 Regulation of TP53 Activity through Phosphorylation  92/10592       4/4 PRKAA1,PRKAB2,PRKAG2,STK11 5.330843e-09   7.996264e-08
4 R-HSA-5633007                         Regulation of TP53 Activity 160/10592       4/4 PRKAA1,PRKAB2,PRKAG2,STK11 5.016583e-08   7.524874e-07
  numberOfPopulationSuccesses populationSize sampleSize numberOfSampleSuccesses
1                          29          10592          4                       4
2                          40          10592          4                       4
3                          92          10592          4                       4
4                         160          10592          4                       4
```

## Additional information

