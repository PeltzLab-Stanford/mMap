![](Example%20results/mmap.png)


# Introduction

   Python implementation of Mouse Genome Genotype to Phenotype mapping method for the prioritizing the genetic candidates based on the protein functional and          genomic regulatory features

# installation
   
   pip install -i https://test.pypi.org/simple/ mgMap==0.0.2 --target=/path/where/to/install/ && cd mgmap
            
            
# usage
            
   python3 mgmap.py -g filename.txt (see *Exmaple* datafiles "gene.txt" for exact format)
   
            
   python3 mmgap.py -nc filename.txt (see *Exmaple* datafiles "regulatory.txt" for exact format)
   
   PS: 
   
   The Functional data folder contains a help file for more information
      
   The Regulatory data folder contains a help file for more information

# usage
   For details of output, see [Wiki](https://github.com/PeltzLab-Stanford/mMap/wiki/working-with-mgMap)  

# compatibility

   The package is compatibility with both Windows and Mac setup. 
    
# requirements

   python 3 
    
   R (see R folder for more information)

# tree

   The package can be used to analyse the impact of mutations on the protein functional regions like domains and PTM sites. Also, for regulatory regions                like Promoter, enhancer etc.
    
      File tree
      ./mgmap/
         mgmap.py
         user's input data-file.txt
         /Functional datafiles/Functional-data.txt
         /Functional datafiles/conservation/
         /Functional datafiles/numbers/
         /Functional datafiles/Phenotype.txt
         /Regulatory datafiles/Vista_enhancers_flankinGenes.txt
         /Regulatory datafiles/EPDnew_promoter.txt
         /Regulatory datafiles/TSS_promoter.txt
         /Regulatory datafiles/cpgIsland.txt
         /Regulatory datafiles/insulator-CTCF-binding-sites.txt

# contact

   email: aarslan@stanford.edu 
   
   visit: https://mmap.stanford.edu/
