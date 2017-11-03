#!usr/bin/python

import numpy
import RhoNetwork
import IterativeRank

file_seed= "test_output/Caulobacter_Fang"
matrix_file_name= file_seed+ "_matrix_file.txt"
key_file_name= file_seed+ "_key_file.txt"
##if correlation_matrix_file does not exist
##correlation_matrix, keys= RhoNetwork.makeCorrelationMatrixFromDelimFile(test_data, True, [], matrix_file_name, key_file_name)

##if correlation_matrix_file does exist
correlation_matrix, keys=RhoNetwork.readCorrelationMatrix(matrix_file_name)

##create weight_matrix
GSR_genes=["sigT","phyR", "phyK", "lovK", "lovR", "nepR", "sigU"]
weight_matrix= IterativeRank.makeWeightMatrixFromKeys(keys)
weight_matrix= IterativeRank.markWeightMatrix(weight_matrix, keys, GSR_genes)

##prepare cell cycle network for PageRank
threshold= 0.9
correlation_matrix= RhoNetwork.thresholdMatrix(correlation_matrix, threshold)
correlation_matrix= RhoNetwork.normalizeCorrelationMatrix(correlation_matrix)

##perform iterative Rank
alpha= 0.50
IR=IterativeRank.PageRank(weight_matrix, correlation_matrix, alpha)
IR_out_filename= file_seed+ "PageRank_GSR_alpha50_threshold90.txt"
IterativeRank.writeOutRanks(IR_out_filename, keys, IR)
