#!/root/miniconda3/bin/python
import pickle
import re

data = open("whog","r")
whog_dict = {}
for i in data:
        if i.startswith("["):
                i = i.strip().split("[")
                j = i[1].split("]")
                value = j[0]
        elif re.search(r":",i):
                i = i.strip().split(":  ")
                j = i[1].split(" ")
                for k in j:
                        whog_dict[k] = value
        elif re.search(r"        ",i):
                i = i.lstrip("        ").split("\t")
                for k in j:
                        whog_dict[k] = value
        else:
                continue
data.close()
pickle_file =open("whog.dict","wb")
pickle.dump(whog_dict,pickle_file)
pickle_file.close()

