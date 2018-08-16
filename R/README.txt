
https://rstudio.github.io/reticulate/articles/introduction.html

library(reticulate)

os <- import("os")

#Path to the folder
os$chdir("/Users/ahmedarslan/Downloads/updated_MacWin")

#import package
mod <- import_from_path("mousepycage", path = ".", convert = TRUE)

mod$gene_mutation('genes.txt')

mod$nocode_mutation("regulatory.txt") TODO files-move