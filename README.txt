This is the IterativeRank python library, written by Matthew Z. Tien,

Reference:
Matthew Z. Tien, Aretha Fiebig, Sean Crosson (2017).
Gene network analysis identifies a central post-transcriptional regulator of cellular stress survival
bioRxiv 212902; doi: https://doi.org/10.1101/212902 

The most up-to-date version of this software is
available at https://github.com/mtien/IterativeRank.

The files 'RhoNetwork.py' and 'IterativeRank.py' constitute the IterativeRank
library. They need to be placed into the search path of your python
installation or locally.

The file 'RhoNetwork.py' sets up the correlation coefficient matrix file from
either tab-delimited files or a dictionary. In the case of the reference paper the keys are the
annotated genes inputed into the Rockhopper (Tjaden 2015) analysis file and the values are
the "Expression" values of RNA-seq data from 5 time points post-synchrony of triplicate
Caulobacter crescentus cultures (Fang et. al. 2013).

The file 'make_Caulobacter_cell_cycle_network.py' is a brief example script demonstrating
basic use of the RhoNetwork.py library. It takes in the output of a standard Rockhopper (Tjaden 2015)
"verbose output" file, creates a dictionary data structure where genes are the keys and the 
"Expression" values are the nodes. The example data in the "test_data/" folder contains such a 
Rockhopper output file. The dictionary is then used to create a correlation coefficient
matrix (symmetrical) where the rows/columns of the matrix correspond to a key list object.
Creation of the matrix object requires that a matrix file and key file be created and
subsequently read. Once a matrix object and key file are created, they do not need
to be created again.

The file 'predict_phyR_parameter_exploration.py' is another brief example of how
the parameters of iterative rank can affect the outcomes of which genes are
strongly associated with the gene of interest. The output file is tab-delimited 
and can easily be open in a spreadsheet program.

The file 'GSR_PageRank.py' is another brief example of how the iterative rank
algorithm was used to predict the role of gsrN in the general stress response
of Caulobacter crescentus.


