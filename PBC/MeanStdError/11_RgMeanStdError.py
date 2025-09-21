#Computes mean standard deviation and standard error
#Requires rend data for 36 chains

#Load math, pandas for reading csv, string, os, numpy
import math
import pandas as pd
from string import *
import os
import numpy as np

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='rg'

f=open('mderg.csv',"w")

for i in range (0,21):
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,header=0)
	rg=pos.iloc[:,1]
	size=np.size(rg)
	dim=rg.ndim
	print(size,dim)
	mean=np.mean(rg,axis=0)
	std=np.std(rg,axis=0)
	err=std/math.sqrt(size)
	print(mean,std,err)
	f.write(f"{mean},{std},{err}\n")
print("Exit")
f.close();
