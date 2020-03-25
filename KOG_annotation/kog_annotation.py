#!/root/miniconda3/bin/python

""" stast COG """
""""    ignoring partial data      """

__date__ = "2020-3-24"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"

import re
import os
import optparse
import pickle
from optparse import OptionParser
from collections import Counter


parser = OptionParser("Usage:%prog --help")

parser.add_option("-d","--database",dest="kog",help="kog dict")
parser.add_option("-f","--fun",dest="fun",help="fun.dict")
parser.add_option("-b","--blast",dest="blast",help="blast output")
parser.add_option("-o","--output",dest="out",help="output file")
(options,args) = parser.parse_args()
kog = open(options.kog,"rb")
fun = open(options.fun,"rb")

kog_dict = pickle.load(kog)
fun_dict = pickle.load(fun)
#print(fun_dict)
out = open(options.out,"w")
out.write('Query_id\tSubject_id\tIdentity(%)\tAlign_length\tMismatch\tGap\tQuery_start\tQuery_end\tSubject_start\tSubject_end\tE-value\tScore\tclassification\tdescription\n')
blast = open(options.blast,"r")
blast_id = set()
for i in blast:
	line = i.strip()
	i = i.strip().split("\t")
	classification = ''
	description = ''
	if i[0] not in blast_id:
		blast_id.add(i[0])
		if i[1] in kog_dict.keys():
			for i in kog_dict[i[1]]:
				classification += i 
				description += fun_dict[i] + ';'
		else:
			classification = "None"
			description = "None"
		out.write(line + '\t' + classification + '\t' + description + '\n' )
	else:
		continue
blast.close()		
out.close()
data = open(options.out,"r")
result = open("stast",'w')
result.write("mark" + "\t" + "description" + "\t" + "number" + "\n")
num_stat = []
data.readline()
for i in data:
	i = i.strip().split("\t")
	j = i[13].split(";")
	for i in j:
		if i != '':

			num_stat.append(i)
data.close
des = set(num_stat)
mark = ord("A")
for i in des:
	num = num_stat.count(i)
	result.write(chr(mark) + "\t" + chr(mark) + ":" + i + "\t" + str(num) + "\n")
	mark += 1
result.close()


