# Advance Bioinformatics Course

About the project: The final group project for this course aimed to assess whether uploading the TAIR10 reference genome to Araport11 unveils additional highly differentially expressed genes, providing insights into the potential benefits of updating the genome annotation.

## Project structure

- YAML *[config](https://git.wur.nl/groen-group/arabidopsis-project/-/tree/main/config?ref_type=heads)* contains the *[config.yaml](https://git.wur.nl/groen-group/arabidopsis-project/-/blob/main/config/config.yaml?ref_type=heads)* with the indices, genomes, metadata and sample calls for the *[snakefile](https://git.wur.nl/groen-group/arabidopsis-project/-/blob/main/workflow/Snakefile?ref_type=heads)* workflow 

- In *[workflow](https://git.wur.nl/groen-group/arabidopsis-project/-/tree/main/workflow?ref_type=heads)* file you can find the *[snakefile](https://git.wur.nl/groen-group/arabidopsis-project/-/blob/main/workflow/Snakefile?ref_type=heads)*, *[R_files](https://git.wur.nl/groen-group/arabidopsis-project/-/tree/main/workflow/R_files?ref_type=heads)*, *[python_files](https://git.wur.nl/groen-group/arabidopsis-project/-/tree/main/workflow/python_files?ref_type=heads)* and *[transcriptome](https://git.wur.nl/groen-group/arabidopsis-project/-/blob/main/workflow/transcriptome?ref_type=heads)* data used in this project.

## Library installation for shared usage

- In *[install_shared_packages.R](https://git.wur.nl/groen-group/arabidopsis-project/-/blob/main/workflow/R_files/install_shared_packages.R)* file you can find the scripts to install the shared packages needed for the different scripts in Rstudio

## Data files

- *[metadata.txt](https://git.wur.nl/groen-group/arabidopsis-project/-/blob/main/metadata.txt?ref_type=heads)* contains the metadata conditions for treatments used by the authors
