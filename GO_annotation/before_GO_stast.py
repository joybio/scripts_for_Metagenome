#!/root/miniconda/bin/python
import optparse
from optparse import OptionParser

parser = OptionParser("Usage: %prog -i GO_annotation -o GO_forstat")
parser.add_option("-i", "--input",dest="input",help="Input file: GO_annotation")
parser.add_option("-o","--out",dest="out",help="Output file: GO_forstat")
(options,args) = parser.parse_args()

def before_GOstat(f1,f2):
	for i in f1.readlines():
		j = i.split("\t")
		for k in j[1].split(","):
			m = j[0] +"\t" +k
			if (m[-1] != "\n"):
				m += "\n"
				f2.write(m)
f1 = open(options.input,"r")
f2 = open(options.out,"w")
before_GOstat(f1,f2)
f1.close()
f2.close()

