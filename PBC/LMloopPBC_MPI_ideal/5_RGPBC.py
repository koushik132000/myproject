#Computes radius of gyration after translating the coordinated by 0.5
#The analysis is for free chains
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains

#Load math, pandas for reading csv, string, os
import math
import pandas as pd
from string import *
import os
import array

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='mchains_rank'
icom='com_rank'
opname='rg_rank'

xmax=7.0
ymax=7.0
zmax=20.0

xms=xmax+0.5
yms=ymax+0.5
zms=zmax+0.5

for i in range (0,41):
    # Intialize arrays of doubles
	sxs=array.array('d',(0.0 for c in range(0,36)))
	sys=array.array('d',(0.0 for c in range(0,36)))
	szs=array.array('d',(0.0 for c in range(0,36)))
	crg=array.array('d',(0.0 for c in range(0,36)))
	
    # Read the center of mass and position files
	infcom=icom+str(rank)+'_'+str(i)+'.'+'csv'
	com=pd.read_csv(infcom,dtype={"c":int,"x":float,"y":float,"z":float})
	comx=com[['x']]
	comy=com[['y']]
	comz=com[['z']]
	infname=ipname+str(rank)+'_'+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,dtype={"x":float,"y":float,"z":float})
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5
	
    # Open the output file for writing
	opfname=opname+str(rank)+'_'+str(i)+'.'+'csv'
	f=open(opfname,"w")
	f.write(f"c,Rg\n")
 
    # Calculating the radius of gyration for each chain
	for c in range(36):
		for b in range(20):
			xdiff=comx.iloc[c,0]-xs.iloc[c*20+b,0]
			ydiff=comy.iloc[c,0]-ys.iloc[c*20+b,0]
			zdiff=comz.iloc[c,0]-zs.iloc[c*20+b,0]

            # Applying periodic boundary conditions
			if (xdiff>xms*0.5):
				xdiff-=xms
			if (xdiff<=-xms*0.5):
				xdiff+=xms
			if (ydiff>yms*0.5):
				ydiff-=yms
			if (ydiff<=-yms*0.5):
				ydiff+=yms
			if (zdiff>zms*0.5):
				zdiff-=zms
			if (zdiff<=-zms*0.5):
				zdiff+=zms

			sxs[c]+=xdiff*xdiff
			sys[c]+=ydiff*ydiff
			szs[c]+=zdiff*zdiff
		print(sxs[c],sys[c],szs[c])			
	for c in range(36):
		crg[c]=math.sqrt(sxs[c]/20+sys[c]/20+szs[c]/20)
        # Writing the results to the output file
		f.write(f"{c},{crg[c]}\n")
	print("Exit")
	f.close();
