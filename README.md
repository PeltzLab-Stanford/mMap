# Introduction

   Python implementation of Mouse Genotype to Phenotype mapping method for the prioritizing the genetic candidates based on the protein functional and genomic regulatory features

# use

    The package can be used to analyse the impact of mutations on the protein functional regions like domains and PTM sites. Also, for regulatory regions like  Promoter, enhancer etc.
    
      File tree
      ./mmap/
         mmap.py
         Functional.txt
         Phenotype.txt
         ./Regulatory/
            Vista_enhancers_flankinGenes.txt
            EPDnew_promoter.txt
            TSS_promoter.txt
            cpgIsland.txt
            insulator-CTCF-binding-sites.txt
            
   # instructions
   python3 -g 'filename.txt' #in the Exmaple datafiles "gene.txt" format
   
   python3 -nc 'filename.txt'  #in the Exmaple datafiles "regulatory.txt" format
   
   PS: The Functional data folder contains a help file for more information
      
      The Regulatory data folder contains a help file for more information
         
# compatibility

    The package is compatibility with both Windows and Mac setup. 
    
# requirements

    python 3 
    
    R (see R folder for more information)

# contact

   email: aarslan@stanford.edu 
   
   visit: http://med.stanford.edu/peltzlab/
