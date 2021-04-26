#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import optparse
import re
from optparse import OptionParser
import pickle



parser = OptionParser("Usage: %prog -i output_file.emapper.annotations -k kegg -l kegg.list -o kegg.ko.list")
parser.add_option("-i", "--input",dest="input",help="Input file: output_file.emapper.annotations")
parser.add_option("-k", "--kegg",dest="kegg",help="tmp file: KEGG; GENE_id:ko_id; 1:>=1")
parser.add_option("-l", "--list",dest="list",help="tmp file: format KEGG; GENE_id:ko_id; 1:1")
parser.add_option("-o","--out",dest="out",help="Output file: GO_forstat")
(options,args) = parser.parse_args()

file = open(options.input, "r")
kegg = open(options.kegg, "w")
lines = file.readlines()
for line in lines:
    line=line.strip()  
    if line.startswith("#"):
        continue
    else:
        tmp=line.split("\t")
        if tmp[9] == "":
            continue
        else :
             kegg.write(tmp[0]+"\t"+tmp[9] + "\n")

file.close()
kegg.close()
def before_GOstat(f1,f2):
	for i in f1.readlines():
		j = i.split("\t")
		for k in j[1].split(","):
			if re.match("map",k):
				pass
			else:
				m = j[0] +"\t" +k
				if (m[-1] != "\n"):
					m += "\n"
					f2.write(m)
f1 = open(options.kegg,"r")
f2 = open(options.list,"w")
before_GOstat(f1,f2)
f1.close()
f2.close()

GENE_id_path_id_list = open(options.list,"r")
out = open(options.out,"w")
data = open("/media/ruizhi/database/KEGG/kegg/genes/fasta/script/path_descrip_dict","rb")
data_dict = pickle.load(data)
out.write("gene_id\tpathway_id\tpathway_description\n")
for i in GENE_id_path_id_list:
	line = i.strip()
	i = i.strip().split("\t")
	ko_id = i[1]
	if ko_id in data_dict.keys():
		out.write(line + "\t" + data_dict[ko_id] + "\n")

GENE_id_path_id_list.close()
out.close()
data.close()


