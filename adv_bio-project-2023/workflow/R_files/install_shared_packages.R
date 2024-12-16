## LIBRARY INSTALLATION FOR SHARED USAGE

# !! please don't run this code again, they are installed now !!

# INSTALL TXIMPORT AND EDGER
#BiocManager::install("tximportData")
BiocManager::install("tximport")
BiocManager::install("edgeR")
BiocManager::install('readr')
BiocManager::install('rhdf5')

Path_shared_dir <- "/prog/BIF30806/project/groen_team/R_shared_libraries"
dir.exists(Path_shared_dir) # should return true
.libPaths(c(Path_shared_dir, .libPaths())) # set path as first in row

## BIOCONDUCTOR ITSELF

# check if Bioconductor itself is installed
isitinstalled <- require("BiocManager", lib.loc = Path_shared_dir, 
                         quietly = TRUE)
isitinstalled

# if not(!) isitinstalled:
if(!isitinstalled){
  install.packages("BiocManager", lib=Path_shared_dir)
}

## BIOCONDUCTOR SUBSIDIARIES
# re-install(force=TRUE) packages to the shared directory(lib) 

# doable programmes
BiocManager::install(version = "3.18", force = TRUE, lib=Path_shared_dir)
BiocManager::install("tximport", force=TRUE, lib=Path_shared_dir)
BiocManager::install("edgeR", force=TRUE, lib=Path_shared_dir)

# large programmes
BiocManager::install('readr', force=TRUE, lib=Path_shared_dir)
BiocManager::install('rhdf5', force=TRUE, lib=Path_shared_dir)
BiocManager::install("tximportData", force=TRUE, lib=Path_shared_dir)
