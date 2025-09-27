#Computes mean standard deviation and standard error
#Requires rend data for 36 chains
#output = R_rms, std_Rrms, err_Rrms for each file 
#Modified Date:24Aug2025
	#21-08-2025: Corrected to Root mean square for end to end distance.
	#24-08-2025: Corrected standard deviation and standard error (R_rms)
	#26-09-2025: correcting the pandas to read header of rend files and using iloc for slicing
	#27-09-2025: included function to compute R_rms, error and its components

#Load math, pandas for reading csv, string, os, numpy
import math
import pandas as pd
from string import *
import os
import numpy as np

# subroutine
# rootmeansq - to compute root mean square and its error
# rmscomp - to compute component rms and its error

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='rend'

f=open('mderend.csv',"w")
f.write(f"R_rms,Rx_rms,Ry_rms,Rz_rms,err_Rrms,err_Rx,err_Ry,err_Rz\n")

# to compute R_rms and its error
def rootmeansq(rend):
    size = len(rend)
    mean_sq = np.mean(rend**2)
    Rrms = math.sqrt(mean_sq)

    std_R2 = np.std(rend**2,ddof=1)
    err_sq = std_R2 / math.sqrt(size)
    err_Rrms = err_sq / (2 * Rrms)

    return Rrms, err_Rrms

# to compute Rx_rms, Ry_rms, Rz_rms and their errors
def rmscomp(r2):
    size = len(r2)
    mean_r2 = np.mean(r2)
    Rcomp = math.sqrt(mean_r2)

    std_r2 = np.std(r2, ddof=1)
    err_r2 = std_r2 / math.sqrt(size)
    err_Rcomp = err_r2 / (2 * Rcomp)

    return Rcomp, err_Rcomp


for i in range (0,41):
    # Reading the rend data
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname)
# extract rend column as a series of float
	rend=pos.iloc[0:36,0] 
	rend_x2=pos.iloc[0:36,1] 
	rend_y2=pos.iloc[0:36,2] 
	rend_Z2=pos.iloc[0:36,3] 

	size=np.size(rend)
	dim=rend.ndim
	print(size,dim)

	R_rms, err_Rrms = rootmeansq(rend)
	Rx_rms, err_Rx = rmscomp(rend_x2)
	Ry_rms, err_Ry = rmscomp(rend_y2)
	Rz_rms, err_Rz = rmscomp(rend_Z2)
	print(R_rms,Rx_rms,Ry_rms,Rz_rms,err_Rrms,err_Rx,err_Ry,err_Rz)
	f.write(f"{R_rms:.6f},{Rx_rms:.6f},{Ry_rms:.6f},{Rz_rms:.6f},{err_Rrms:.6f},{err_Rx:.6f},{err_Ry:.6f},{err_Rz:.6f}\n")
print("Exit")
f.close();