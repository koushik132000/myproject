#Computes mean standard deviation and standard error
#Requires rend data for 36 chains
#output = Rg_rms, std_Rgrms, err_Rgrms for each file 
#Modified Date:24Aug2025
	#21-08-2025: Corrected to Root mean square for Radius of gyration.
	#24-08-2025: Corrected standard deviation and standard error (Rg_rms)
	#26-09-2025: correcting the pandas to read header of rg files and using iloc for slicing
	#27-09-2025: included function to compute Rg_rms, error and its components

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

ipname='rg'

f=open('mderg.csv',"w")
f.write(f"Rg_rms,Rgx_rms,Rgy_rms,Rgz_rms,err_Rgrms,err_Rgx,err_Rgy,err_Rgz\n")

# to compute Rg_rms and its error
def rootmeansq(rg):
    size = len(rg)
    mean_sq = np.mean(rg**2)
    Rgrms = math.sqrt(mean_sq)

    std_R2 = np.std(rg**2,ddof=1)
    err_sq = std_R2 / math.sqrt(size)
    err_Rgrms = err_sq / (2 * Rgrms)

    return Rgrms, err_Rgrms

# to compute Rgx_rms, Rgy_rms, Rgz_rms and their errors
def rmscomp(r2):
    size = len(r2)
    mean_r2 = np.mean(r2)
    Rgcomp = math.sqrt(mean_r2)

    std_r2 = np.std(r2, ddof=1)
    err_r2 = std_r2 / math.sqrt(size)
    err_Rgcomp = err_r2 / (2 * Rgcomp)

    return Rgcomp, err_Rgcomp


for i in range (0,41):
    # Reading the rend data
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname) #default reads header
# extract rend column as a series of float
	rg=pos.iloc[0:36,0] 
	rg_x2=pos.iloc[0:36,1] 
	rg_y2=pos.iloc[0:36,2] 
	rg_Z2=pos.iloc[0:36,3] 

	size=np.size(rg)
	dim=rg.ndim
	print(size,dim)

	Rg_rms, err_Rgrms = rootmeansq(rg)
	Rgx_rms, err_Rgx = rmscomp(rg_x2)
	Rgy_rms, err_Rgy = rmscomp(rg_y2)
	Rgz_rms, err_Rgz = rmscomp(rg_Z2)
	print(Rg_rms,Rgx_rms,Rgy_rms,Rgz_rms,err_Rgrms,err_Rgx,err_Rgy,err_Rgz)
	f.write(f"{Rg_rms:.6f},{Rgx_rms:.6f},{Rgy_rms:.6f},{Rgz_rms:.6f},{err_Rgrms:.6f},{err_Rgx:.6f},{err_Rgy:.6f},{err_Rgz:.6f}\n")
print("Exit")
f.close();