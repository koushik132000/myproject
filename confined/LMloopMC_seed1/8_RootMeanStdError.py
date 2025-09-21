#Computes mean standard deviation and standard error
#Requires rend data for 36 chains
#output = R_rms, std_Rrms, err_Rrms for each file 
#Modified Date:24Aug2025
	#21-08-2025: Corrected to Root mean square for end to end distance.
	#24-08-2025: Corrected standard deviation and standard error (R_rms)

#Load math, pandas for reading csv, string, os, numpy
import math
import pandas as pd
from string import *
import os
import numpy as np

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='rend'

f=open('mderend.csv',"w")

for i in range (0,41):
    # Reading the rend data
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,header=None)
	rend=pos[0] # extract rend column as a series
	size=np.size(rend)
	dim=rend.ndim
	print(size,dim)

	# calculating Root-mean-square(R_rms)
	mean_sq = np.mean(rend**2, axis=0)
	R_rms = math.sqrt(mean_sq)

	# calculating standard deviation and standard error of R_rms
	std_R2 = np.std(rend**2, axis=0)# std deviation of R^2
	std_Rrms = std_R2 /(2 * R_rms)
	err_sq = std_R2 / math.sqrt(size)# standard error of R^2
	err_Rrms = err_sq /(2 * R_rms)
	print(R_rms,std_Rrms,err_Rrms)
	f.write(f"{R_rms:.6f},{std_Rrms:.6f},{err_Rrms:.6f}\n")
print("Exit")
f.close();