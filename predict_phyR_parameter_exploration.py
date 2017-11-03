#!usr/bin/python

import RhoNetwork
import IterativeRank
	
def chooseParameters(inquiry, mark_genes, corr_keys, corr_matrix, thresholds, alphas, out_file):
	weight_matrix= IterativeRank.makeWeightMatrixFromKeys(corr_keys)
	weight_matrix= IterativeRank.markWeightMatrix(weight_matrix, corr_keys, mark_genes)
	out_file=open(out_file,"w")
	out_file.write("threshold\talpha_value\trank\tscore\n")
	for threshold in thresholds:
		threshold_matrix= corr_matrix.copy()
		threshold_matrix= RhoNetwork.thresholdMatrix(threshold_matrix, threshold)
		threshold_matrix= RhoNetwork.normalizeCorrelationMatrix(threshold_matrix)
		for alpha in alphas:
			alpha_IR=IterativeRank.PageRank(weight_matrix, threshold_matrix, alpha)
			alpha_dic=IterativeRank.getRanks(inquiry, corr_keys, alpha_IR)
			dic_keys= alpha_dic.keys()
			for d_k in dic_keys:
				rank= alpha_dic[d_k][0]
				score= alpha_dic[d_k][1]
				out_file.write(str(threshold) +"\t"+ str(alpha) +"\t"+ str(rank) +"\t"+ str(score) +"\n")
	out_file.close()
	
if __name__ == "__main__":
	file_seed= "test_output/Caulobacter_Fang"
	matrix_file_name= file_seed+ "_matrix_file.txt"
	key_file_name= file_seed+ "_key_file.txt"
	correlation_matrix, keys= RhoNetwork.readCorrelationMatrix(matrix_file_name, key_file_name)
	
	alpha_values=[.99, .95, .9, .8, .75, .7, .6, .5, .4, .3, .25, .2, .1, .05, .01, .001]
	threshold_values= [0.0, .05, .25, .5, .75, .9, .95]
	marked_genes=["sigT", "phyK", "lovK", "lovR", "nepR", "sigU"]
	inquiry_gene=["phyR"]
	
	file_out_name= file_seed + "_predict_phyR_parameter_exploration.txt"
	chooseParameters(inquiry_gene, marked_genes, keys, correlation_matrix, threshold_values, alpha_values, file_out_name)
