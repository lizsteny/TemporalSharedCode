## (I) Data origin <br/><br/>

All genome sequences were obtained from genbank (https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/Database/nph-select.cgi).

**Outputs:**
- *denv_nu_raw.fasta*(nucleotide genome sequences).
- *nucleotideaccesionlist.fasta* (genbank accessions).
- *denv4_tableraw.csv* (353 unique nucleotide sequences with id, location, description and collection date).
  


## (II) Sequence removal with unknown information<br/><br/>

We removed the genomic sequences whose year and sampling site was not described.
Also in *Dengue Virus Typing Tool v3.80* we identify the genotypes corresponding to each genomic sequence and eliminate those that could not be identified.

**Outputs:**
- *denv4_genotypes.csv* (genotypes assigned by Dengue Virus Typing Tool).
- *denv4_table.csv* (351 unique nucleotide sequences with id, location, description, collection date and genotype).
- *denv4_seqs.fas*(351 unique genomic nucleotide sequences).

**Requirements:**
- *Dengue Virus Typing Tool v3.80* (can be used in https://www.genomedetective.com/app/typingtool/dengue/introduction).

## (III) Alignment <br/><br/>

With *MUSCLE v3.8.31* alignments were performed.
In *Mesquite v3.61* sequences with missing data were removed and manual adjustments were made. Also a protein matrix was created from DNA.

**Outputs:**
- *denv4_nu.fasta*(341 aligned genomic nucleotide sequences).
- *denv4_aa.fasta*(341 aligned genomic aminoacid sequences).

**Requirements:**
- *MUSCLE v3.8.31* (https://www.drive5.com/muscle/downloads.htm).
- *denv4_seqs.fas*(351 genomic nucleotide sequences).

## (IV) Epitopes<br/><br/>

All epitopes were obtained from The Immune Epitope Database (IEDB) (https://www.iedb.org/home_v3.php).

**Outputs:**
- *epitopes_raw.csv*(epitopes table description).

