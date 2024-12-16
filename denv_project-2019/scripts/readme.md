## (I) Diversity index<br/><br/>

Script to obtain the shannon diversity index and the number of polymorphisms in the genome sequence of amino acids. 

**Outputs:**
- *myaa.csv* (transformed tcv data to csv).
- *mydiversity.csv* (shannon diversity index, normalized shannon diversity index and the number of polymorphisms in the genome  sequences by position).

**Requirements:**
- *denv4_aa.tab* stored in *data* folder (genomic aminoacid sequences in tcv format).
- R dependencies: *vegan*, *ggplot2*, *stringr*.

`we did the same for denv4_table_genotypeI.csv and denv4_table_genotypeII.csv`<br/><br/>

**Outputs:**
- *mydiversity_genI.csv* (shannon diversity index, normalized shannon diversity index and the number of polymorphisms in the genome  sequences by position).
- *mydiversity_genII.csv* (shannon diversity index, normalized shannon diversity index and the number of polymorphisms in the genome  sequences by position).

**Requirements:**
-  *denv4_table_genotypeI.csv* stored in *data* folder (genomic aminoacid sequences of genotype I in tcv format).
-  *denv4_table_genotypeI.csv* stored in *data* folder (genomic aminoacid sequences of genotype II in tcv format).
- R dependencies: *vegan*, *ggplot2*, *stringr*.

## (II) Epitopes<br/><br/>

Script to obatin how many epitopes are in each position of the genome sequence and its frequency. 

**Outputs:**
- *epitopesfreq.csv* (epitopes frequency at each position of denv genome sequences).

**Requirements:**
- *epitopes.csv* stored in *data* folder (epitopes only with starting and ending position).
- R dependencies: *dplyr*, *ggplot2*.

## (III) Denv4 table plots<br/><br/>

Script to general informative plots like collection date of the genomic sequences per country by genotype and number of genomic sequences per year by genotype.

**Requirements:**
- *denv4_table.csv* stored in *data* folder.
- R dependencies: *ggplot2*, *reshape2*.
