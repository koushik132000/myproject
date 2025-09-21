#Computes end-to-end distance after translating all coordinates by 0.5
#The analysis is for confined chains
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains
#Modified Date: 05Aug2025
		#10July2025: Added MPI parallelization
		#05Aug2025: Corrected input and output file names to include rank
#Load math, pandas for reading csv, string, os
import math
import pandas as pd
from string import *
import os

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='mchains_rank'
opname='rend_rank'

# Reading data from the input files
for i in range (0,41):
	infname=ipname+str(rank)+'_'+str(i)+'.'+'csv'
	pos=pd.read_csv(infname)
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5

    # Creating a file 
	opfname=opname+str(rank)+'_'+str(i)+'.'+'csv'
	f=open(opfname,"w")
  
    # Calculating the end-to-end distance
	for c in range(36):
		x0=pos.at[20*c,'x']
		y0=pos.at[20*c,'y']
		z0=pos.at[20*c,'z']
		print(x0,y0,z0)
		xl=pos.at[20*c+19,'x']
		yl=pos.at[20*c+19,'y']
		zl=pos.at[20*c+19,'z']
		print(xl,yl,zl)

		xdiff=xl-x0
		ydiff=yl-y0
		zdiff=zl-z0

		sq=math.sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)
  
        # writing the end-to-end distance to the file 
		f.write(f"{sq}\n")
		print(sq)
	print("Exit")
	f.close();
