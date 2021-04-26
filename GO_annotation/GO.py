#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os
x = sys.argv[1]
file = open(x, "r")
lines = file.readlines()
for line in lines:
    line=line.strip()  
    if line.startswith("#"):
        continue
    else:
        tmp=line.split("\t")
        if tmp[6] == "":
            continue
        else :
             mystr=tmp[0]+"\t"+tmp[6]
             print (mystr)

