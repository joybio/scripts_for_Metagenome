#!/root/miniconda3/bin/python
import pickle
import re

data = open("kog","r")
kog_dict = {}
for i in data:
	if i.startswith("["):
		i = i.strip().split("[")
		j = i[1].split("]")
		value = j[0]
	elif re.search(":",i):
		i = i.strip().split(":")
		j = i[1].lstrip().split("_")
		k = j[0]
		kog_dict[k] = value
	else:
		continue
data.close()
pickle_file =open("kog.dict","wb")
pickle.dump(kog_dict,pickle_file)
pickle_file.close()

