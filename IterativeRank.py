#!usr/bin/python
'''This module is part of the IterativeRank library,
written by Matthew Z. Tien

The IterativeRank module contains code to compute
the outcomes of implementing the iterative rank algorithm
from objects created in RhoNetwork. It requires the RhoNetwork module
(also part of the IterativeRank library), which contains
the necessary tools to create a correlation coefficient matrix
and a corresponding key list.
This module also requires the numpy.linalg module from
numpy, for matrix computations.'''

import numpy
import RhoNetwork
from numpy.linalg import inv

def makeWeightMatrixFromMatrix(matrix):
	length=math.sqrt(matrix.size)
	return numpy.zeros((length,1))
	
def makeWeightMatrixFromKeys(keys):
	return numpy.zeros((len(keys),1))

def markWeightMatrix(weight_matrix, keys, marked_keys):
	for i in range(0,len(keys)):
		key= keys[i]
		if(key in marked_keys):
			weight_matrix.itemset(i, 1)
	return weight_matrix
	
def	writeOutRanks(file_out, keys, key_weights):
	temp_dic= {}
	for i in range(0,len(keys)):
		key= keys[i]
		key_weight= key_weights.item(i)
		temp_dic[key]= key_weight
		
	sorted_keys= sorted( temp_dic, key=lambda x: temp_dic[x], reverse=True) 
	
	file_out= open(file_out, "w")
	file_out.write("Rank\tKeyName\tValue\n")
	rank=1
	for s in sorted_keys:
		file_out.write(str(rank) + "\t" + s + "\t" + str(temp_dic[s]) + "\n")
		rank+=1
	file_out.close()
	
def	getRanks(inquiry_keys, keys, key_weights):
	temp_dic= {}
	for i in range(0,len(keys)):
		key= keys[i]
		key_weight= key_weights.item(i)
		temp_dic[key]= key_weight
		
	sorted_keys= sorted( temp_dic, key=lambda x: temp_dic[x], reverse=True) 
	
	rank=1
	inquiry_ranks={}
	for s in sorted_keys:
		if(s in inquiry_keys):
			inquiry_ranks[s]= [rank, temp_dic[s]]
		rank+=1
	return inquiry_ranks

def PageRank(initial_weight, correlation_matrix, alpha):
	identity_matrix= numpy.identity(initial_weight.size)
	scaled_correlation_matrix= (1.0-alpha)*correlation_matrix
	difference= inv(identity_matrix-scaled_correlation_matrix)
	return alpha*difference*initial_weight
