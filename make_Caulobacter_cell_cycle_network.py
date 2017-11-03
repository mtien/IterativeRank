#!usr/bin/python

import RhoNetwork

def writeOutDictionary( file_output_path, headers, dictionary ):
	file_o= open( file_output_path, "w")
    
	for h in range(0, len(headers)-1):
		file_o.write( headers[h] +"\t" )
	file_o.write( headers[-1] + "\n" )

	key_list= dictionary.keys()
	for key in key_list:
 		file_o.write(key)
		for key_comps in dictionary[key]:
			file_o.write( "\t" + key_comps )
		file_o.write("\n")

	file_o.close()
	return file_output_path
	
def normalizeValues(dictionary):
	keys= list(dictionary.keys())
	new_dictionary={}
	for k in keys:
		values= dictionary[k]
		float_values=[]
		for v in values:
			float_values.append(float(v))
		norm_value= sum(float_values)
		if(norm_value > 0.0 ):
			new_values=[]
			for val in float_values:
				new_values.append(str(val/norm_value))
			new_dictionary[k]=new_values
	return new_dictionary

def parseByRhopperHeader(Rhopper_file_name, header_criterium):
	file_in= open(Rhopper_file_name)
	header=file_in.readline().split("\t")
	index_counter=0
	indicies=[]
	out_header=["CCNA"]
	for h in header:
		for c in header_criterium:
			if(h.find(c)>-1):
				indicies.append(index_counter)
				out_header.append(h.strip())
				break
		index_counter+=1
	
	output_dictionary={}
	for line in file_in:
		info=line.split("\t")
		if(len(info) > 1):
			name=info[5]
			CCNA=info[6]
			if( CCNA.find("CCNA") >-1):
				if(name!="-"):
					CCNA= name
				values=[]
				for index in indicies:
					values.append(info[index].strip())
				output_dictionary[CCNA]=values
	file_in.close()
	return out_header, output_dictionary

if __name__== "__main__":
	Rockhopper_output_file="test_data/CP001340_transcripts_Fang.txt"
	criterium= ["Expression"]
	
	temp_header, temp_dic= parseByRhopperHeader(Rockhopper_output_file, criterium)
	temp_dic=normalizeValues(temp_dic)
	
	norm_expression_file= "test_data/Fang_Rhopper_NormExpression.txt"
	writeOutDictionary(norm_expression_file, temp_header, temp_dic)
	
	file_seed= "test_output/Caulobacter_Fang"
	matrix_file_name= file_seed+ "_matrix_file.txt"
	key_file_name= file_seed+ "_key_file.txt"
	
	RhoNetwork.makeCorrelationMatrixFromDelimFile(norm_expression_file, True, [], matrix_file_name, key_file_name)
	
	##alternatively
	##RhoNetwork.makeCorrelationMatrixFromDictionary(temp_dic, [], matrix_file_name, key_file_name)
