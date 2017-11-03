#!usr/bin/python
'''This module is part of the IterativeRank library,
written by Matthew Z. Tien

The RhoNetwork module sets up the correlation coefficient matrix file from
either tab-delimited files or a dictionary. It also creates a corresponding 
key list that corresponds the to rows/columns of the correlation coefficient
matrix.'''

from __future__ import print_function
import numpy
import scipy.stats

def computePearsonCorrelation(gene_array_1, gene_array_2):
	"""
	Returns float:	Pearson correlation coefficient between two input gene arrays scipy.stats.pearsonr( val_1, val_2 )[0]
							returns a float from -1 to 1 if pearson correlation coefficient can be calculated
							returns a float of 0.0 if pearson correlation coefficient cannot be calculated
	
			gene_array_X:	arrays of float, string, or char types
				
			description:	The function will first sift out compatible values in the index of the arrays. 
							The correlation coefficient will only consider numeric values
							If values for either gene arrays are non-numeric, the index from both arrays will not be used to calculate the correlation coefficient
	"""
	arr_1_length= len(gene_array_1)
	arr_2_length= len(gene_array_2)
	if( arr_1_length == arr_2_length):
		val_1=[]
		val_2=[]
		for i in range(0, arr_1_length):
			if( isinstance(gene_array_1[i], float) and isinstance(gene_array_2[i], float)):
				val_1.append(gene_array_1[i])
				val_2.append(gene_array_2[i])
		if(len(val_1) > 1 ):
			pearson= scipy.stats.pearsonr( val_1, val_2 )[0]
			if(pearson < 0):
				return 0.0
			elif( numpy.isnan(pearson)):
				return 0.0
			else:
				return pearson
		else:
			return 0.0
	else:
		return 0.0
	
def makeKeyFile(keys=[],key_file_name= "temporary_key_file.txt"):
	temp_key_file= open(key_file_name, "w")
	for k in keys:
		temp_key_file.write(k+ ",")
	temp_key_file.close()
	return len(keys), key_file_name

def readKeyFile(key_file_name):
	key_file= open(key_file_name)
	keys= key_file.readline().strip().split(",")[:-1]
	key_file.close()
	return keys, len(keys)

def getKeyFileFromMatrixFile(matrix_file_name):
	file_in=open(matrix_file_name)
	header=file_in.readline()
	header_info= header.strip().split()
	number_of_keys= int(header_info[0].split(",")[0])
	key_file_name= header_info[0].split(",")[1]
	file_in.close()
	return key_file_name, number_of_keys

def makeDictionaryFromDelimFile(file_name, column_header=True):
	file_input= open(file_name)
	if(column_header):
		header=file_input.readline()
	output_dictionary={}
	for line in file_input:
		information= line.strip().split("\t")
		key= information[0]
		values=[]
		for i in range(1,len(information)):
			entry= information[i]
			try:
				values.append(float(entry))
			except ValueError:
				values.append( entry )
		output_dictionary[key]= values
	return output_dictionary

def readCorrelationMatrix(matrix_file_name, key_file_name=""):
	if(key_file_name==""):
		key_file_name, number_of_keys= getKeyFileFromMatrixFile(matrix_file_name)
	keys, number_of_keys= readKeyFile(key_file_name)
	file_in=open(matrix_file_name)
	header=file_in.readline()
	
	matrix_string= ""
	counter=1
	for line in file_in:
		n=0
		info= line.strip().split("\t")
		for entry in info:
			n+=1
			matrix_string+=entry+ " "
		zeros= numpy.zeros(number_of_keys-counter)
		for zero in zeros:
			n+=1
			matrix_string+=str(zero) + " "
		matrix_string= matrix_string[:-1] + ";"
		counter+=1
	file_in.close()
	print("finished reading file")
	
	lower_triangle=numpy.matrix(matrix_string[:-1])
	upper_triangle= lower_triangle.T
	complete_matrix= lower_triangle+upper_triangle - numpy.identity(number_of_keys)
	print("created correlation matrix")
	return complete_matrix, keys

def makeArrayIntoMatrixRow(input_array, delim= "\t"):	
	n= len(input_array)
	output_string=""
	for i in range(0,n-1):
		output_string+= str(input_array[i]) + delim
	output_string+= str(input_array[n-1]) + "\n"
	return output_string
	
def makeCorrelationMatrixFromDictionary(value_dictionary, keys=[], matrix_file_name="temporary_matrix_file.txt", key_file_name="temporary_key_file.txt"):
	"""
	Returns numpy.matrix:	matrix object of correlation coefficients, symmetric matrix with a diagonal of ones
	Returns list:	list of keys that correspond to the rows/columns of keys
	
	value_dictionary:	{key: [values]}
	keys:	list option where the user may specify a subset of keys in which to make the correlation coefficient matrix from the delim_file_name
	matrix_file_name:	string containing the path to construct the correlation matrix 
						matrix stored as just the lower triangle of a symmetical matrix
						each line is tab-delimited and corresponds to the correlation coefficients of one key to another as listed in the key_file_name
						key_file_name should be in the first line of the matrix_file, where the first entry is the number of keys and the second is the key_file_path 
						line1: number_of_keys,key_file_name
	key_file_name:	string containing the path to the key to each line/column in the correlation matrix
	
	description:	function returns a matrix of correlations between all keys to all keys
	"""
	if(len(keys)==0):
		keys= list(value_dictionary.keys())
	temp_matrix_file= open(matrix_file_name, "w")
	number_of_keys, temp_key_file_name= makeKeyFile(keys, key_file_name)
	temp_matrix_file.write(str(number_of_keys) + "," +temp_key_file_name +"\n")
	for k in range(0, len(keys)):
		main_key=keys[k]
		main_value= value_dictionary[main_key]
		correlation_array=[]
		for k2 in range(0,k+1):
			compared_key= keys[k2]
			compared_value= value_dictionary[compared_key]	
			correlation_array.append(computePearsonCorrelation(main_value, compared_value))
		if(sum(correlation_array)==0.0):
			correlation_array[-1]=1.0
		temp_matrix_file.write(makeArrayIntoMatrixRow(correlation_array))
	temp_matrix_file.close()
	
	return readCorrelationMatrix(matrix_file_name, key_file_name)

def makeCorrelationMatrixFromDelimFile(delim_file_name, header=False, keys=[], temp_matrix_file_name="temporary_matrix_file.txt", temp_key_file_name="temporary_key_file.txt"):
	"""
	Returns numpy.matrix:	matrix object of correlation coefficients, symmetric matrix with a diagonal of ones
	Returns list:	list of keys that correspond to the rows/columns of keys
	
	delim_file_name:	string containing the path to construct the correlation matrix. 
						file should be tab-delimited.
						first column of the delimited file is assumed to be the KEYS of the correlation matrix.
						subsequent columns are the VALUES to compute the correlation coefficient
						If values for certain genes are missing within a experimental set (data within a column), a non-numeric character is entered in it's place for that column
	header:	boolean option indicates if the delim_file_name has a header
	keys: list option where the user may specify a subset of keys in which to make the correlation coefficient matrix from the delim_file_name
	temp_matrix_file_name:	string of path in which the matrix_file_name will be written
	temp_natrix_key_name:	string of path in which the key_file_name will be written  
	"""
	print("make dictionary")
	val_dictionary= makeDictionaryFromDelimFile(delim_file_name, header)
	print("make Matrix")
	return makeCorrelationMatrixFromDictionary(val_dictionary, keys, temp_matrix_file_name, temp_key_file_name)
	
def writeOutCorrelationMatrixFile(numpy_matrix, output_matrix_file_name, key_file_name):
	"""
	Returns string:	path of file where matrix was written
					matrix stored as just the lower triangle of a symmetical matrix
					each line is tab-delimited and corresponds to the correlation coefficients of one key to another as listed in the key_file_name
					key_file_name should be in the first line of the matrix_file, where the first entry is the number of keys and the second is the key_file_path 
					line1: number_of_keys,key_file_name
	
	numpy_matrix:	matrix object to write out
	output_matrix_file_name:	file path in which to write the matrix file
	key_file_name:	file path in which keys will be associated with matrix file
	"""
	keys, number_of_keys= readKeyFile(key_file_name)
	file_out= open(output_matrix_file_name, "w")
	file_out.write(str(number_of_keys) + "," + key_file_name + "\n")
	for i in range(number_of_keys):
		for j in range(0,i):
			file_out.write( str(numpy_matrix.item((i,j))) + "\t" )
		file_out.write(str(numpy_matrix.item((i,i))) + "\n")
	file_out.close()
	return output_matrix_file_name

def normalizeCorrelationMatrix( correlation_matrix ):
	col_sums= correlation_matrix.sum(axis=0)
	col_norm= correlation_matrix/col_sums
	row_sums= col_norm.sum(axis=1)
	norm_matrix= col_norm/row_sums
	return norm_matrix
	
def thresholdMatrix(correlation_matrix, threshold):
	correlation_matrix[correlation_matrix < threshold]= 0.0
	return correlation_matrix
