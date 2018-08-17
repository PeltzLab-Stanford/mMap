
https://rstudio.github.io/reticulate/articles/introduction.html

  File tree
  ./mmap/
     mmapR.py
     user's input data-file.txt
     /Functional datafiles/Functional-data.txt
     /Functional datafiles/Phenotype.txt
     /Regulatory datafiles/Vista_enhancers_flankinGenes.txt
     /Regulatory datafiles/EPDnew_promoter.txt
     /Regulatory datafiles/TSS_promoter.txt
     /Regulatory datafiles/cpgIsland.txt
     /Regulatory datafiles/insulator-CTCF-binding-sites.txt

library(reticulate)

os <- import("os")

#Path to the folder
os$chdir("/path/to/folder/mMap")

#import package
mod <- import_from_path("mmapR", path = ".", convert = TRUE)

mod$gene_mutation('genes.txt') #see example-data folder for gene.txt file format

mod$nocode_mutation("regulatory.txt") #see example-data folder for regulatory.txt file format

PS: See file tree to run mMap successfully
