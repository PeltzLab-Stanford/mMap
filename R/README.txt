
https://rstudio.github.io/reticulate/articles/introduction.html

library(reticulate)

os <- import("os")

#Path to the folder
os$chdir("/path/to/folder/mMap")

#import package
mod <- import_from_path("mmapR", path = ".", convert = TRUE)

mod$gene_mutation('genes.txt') #see example-data folder for gene.txt file format

mod$nocode_mutation("regulatory.txt") #see example-data folder for regulatory.txt file format

PS: See file tree to run mMap successfully
