#!/root/miniconda3/bin/python
import pickle
import optparse
from optparse import OptionParser
import os
import re


parser= OptionParser("Usage:%prog --help")
parser.add_option("-i","--input",dest="input",help="input file: GO_annotation")
parser.add_option("-o","--out",dest="out",help="output file: file name,eg:GO.annotation.stast")
(options,args) = parser.parse_args()
data = open(options.input,"r")
out = open(options.out,"w")
go_des = open("/media/ruizhi/database/GO/go_des_dict","rb")
go_term = open("/media/ruizhi/database/GO/go_term_dict","rb")
go_ont = open("/media/ruizhi/database/GO/go_ont_dict","rb")
go_secondary = open("/media/ruizhi/database/GO/go_secondary_set","rb")
go_secondary_out = open("GO.secondary.annotation","w")
go_secondary_out.write("GO_ID,GENE_NUM,TOTAL_GENE_NUM,DESCRIPTION,ONTOLOGY,TERM,GENE_ID\n")

go_des_dict = pickle.load(go_des)
go_term_dict = pickle.load(go_term)
go_ont_dict = pickle.load(go_ont)
go_secondary.set = pickle.load(go_secondary)

go_id_dict = {}
n = 0
for i in data:
	n += 1
	i = i.strip().split("\t")
	gene_name = i[0]
	go_id = i[1].split(",")
	for j in go_id:
		if j not in go_id_dict.keys():
			go_id_dict[j] = gene_name
		else:
			go_id_dict[j] += "|" + gene_name
data.close()
out.write("GO_ID,GENE_NUM,TOTAL_GENE_NUM,DESCRIPTION,ONTOLOGY,TERM,GENE_ID\n")
for i in go_id_dict.keys():
	gene_list = go_id_dict[i].split("|")
	gene_num = len(gene_list)
	if i in go_des_dict.keys():
		description = go_des_dict[i]
	else:
		description = "NA"
	if i in go_ont_dict.keys():
		ontology = go_ont_dict[i]
	else:
		ontology = "NA"
	if i in go_term_dict.keys():
		term = go_term_dict[i]
	else:
		term = "NA"
	out.write(i + "," + str(gene_num) + "," + str(n) + "," + description + "," + ontology + "," + term + "," + go_id_dict[i] + "\n")
	if i in go_secondary.set:
		go_secondary_out.write(i + "," + str(gene_num) + "," + str(n) + "," + description + "," + ontology + "," + term + "," + go_id_dict[i] + "\n")
out.close()		

os.system('cut -d "," -f 1-6 GO.secondary.annotation > GO.secondary.annotation.csv')
os.system('cut -d "," -f 1-6 GO.annotation.stast > GO.annotation.stast.csv')







