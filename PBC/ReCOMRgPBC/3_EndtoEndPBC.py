#Computes end-to-end distance after translating all coordinates by 0.5
#The analysis is for a code with periodic boundary conditions
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains

#Load math, pandas for reading csv, string, os
import math
import pandas as pd
from string import *
import os

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='mchains'
opname='rend'

xmax=7.0
ymax=7.0
zmax=20.0

# Reading data from the input files
for i in range (0,41):
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname)
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5

    # creating a file
	opfname=opname+str(i)+'.'+'csv'
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

        # Applying periodic boundary conditions for x,y and z.
		if (xdiff>-xmax*0.5):
			xdiff=xdiff-xmax
		if (xdiff<=xmax*0.5):
			xdiff=xdiff+xmax

		if (ydiff>-ymax*0.5):
			ydiff=ydiff-ymax
		if (ydiff<=ymax*0.5):
			ydiff=ydiff+ymax

		if (zdiff>-zmax*0.5):
			zdiff=zdiff-zmax
		if (zdiff<=zmax*0.5):
			zdiff=zdiff+zmax
		
		sq=math.sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)
        # writing the result to the file 
		f.write(f"{sq}\n")
		print(sq)
	print("Exit")
	f.close();
