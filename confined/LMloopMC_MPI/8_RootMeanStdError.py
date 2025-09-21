#Computes mean standard deviation and standard error
#Requires rend data for 36 chains
#Modified Date: 24Aug2025
		#10July2025: Added MPI parallelization
		#05Aug2025: Corrected input and output file names to include rank
		#21-08-2025: Corrected to Root mean square for end to end distance.
		#24-08-2025: Corrected standard deviation and standard error (R_rms)
#Load math, pandas for reading csv, string, os, numpy
import math
import pandas as pd
from string import *
import os
import numpy as np

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='rend_rank'

f=open('mderend'+str(rank)+'.'+'csv',"w")

for i in range (0,41):
    # Reading the rend data
	infname=ipname+str(rank)+'_'+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,header=None)
	rend=pos[0]
	npoints=np.size(rend)
	dim=rend.ndim
	print(npoints,dim)
	mean_sq = np.mean(rend**2, axis=0)# mean of squares
	R_rms = math.sqrt(mean_sq)# root-mean-square
	std_R2 = np.std(rend**2, axis=0)# std deviation of R^2
	std_Rrms = std_R2 /(2 * R_rms)# std deviation of R_rms
	err_sq = std_R2 / math.sqrt(npoints)# standard error of R^2
	err_Rrms = err_sq /(2 * R_rms)# standard error of R_rms 
	print(R_rms,std_Rrms,err_Rrms)
	f.write(f"{R_rms:.6f},{std_Rrms:.6f},{err_Rrms:.6f}\n")
print("Exit")
f.close();
