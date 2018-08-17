# Introduction

   Python implementation of Mouse Genotype to Phenotype mapping method for the prioritizing the genetic candidates based on the protein functional and genomic regulatory features

# use

    The package can be used to analyse the impact of mutations on the protein functional regions like domains and PTM sites. Also, for regulatory regions like  Promoter, enhancer etc.
    
      File tree
      ./mmap/
         mmap.py
         user's input data-file.txt
         /Functional datafiles/Functional-data.txt
         /Functional datafiles/Phenotype.txt
         /Regulatory datafiles/Vista_enhancers_flankinGenes.txt
         /Regulatory datafiles/EPDnew_promoter.txt
         /Regulatory datafiles/TSS_promoter.txt
         /Regulatory datafiles/cpgIsland.txt
         /Regulatory datafiles/insulator-CTCF-binding-sites.txt
            
   # instructions
   python3 -g 'filename.txt' #in the Exmaple datafiles "gene.txt" format
   
   python3 -nc 'filename.txt'  #in the Exmaple datafiles "regulatory.txt" format
   
   PS: 
   
   The Functional data folder contains a help file for more information
      
   The Regulatory data folder contains a help file for more information
         
# compatibility

    The package is compatibility with both Windows and Mac setup. 
    
# requirements

    python 3 
    
    R (see R folder for more information)

# contact

   email: aarslan@stanford.edu 
   
   visit: http://med.stanford.edu/peltzlab/
