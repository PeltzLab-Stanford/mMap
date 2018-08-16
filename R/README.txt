
https://rstudio.github.io/reticulate/articles/introduction.html

library(reticulate)

os <- import("os")

#Path to the folder
os$chdir("/path/to/folder/mMap")

#import package
mod <- import_from_path("mmapR.py", path = ".", convert = TRUE)

mod$gene_mutation('genes.txt')

mod$nocode_mutation("regulatory.txt")

PS: See file tree to run mMap successfully
