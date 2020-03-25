#!/root/miniconda3/bin/python
import pickle
import re
fun_dict = {}
fun = open("fun.txt","r")
for i in fun:
        if i.startswith(" "):
                i = i.strip().split('] ')
                key = str(i[0][1])
                value = i[1]
                fun_dict[key] = value
#print(fun_dict)
fun.close()
pickle_file =open("fun.dict","wb")
pickle.dump(fun_dict,pickle_file)
pickle_file.close()


