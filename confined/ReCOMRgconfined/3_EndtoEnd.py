#Computes end-to-end distance after translating all coordinates by 0.5
#The analysis is for confined chains
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains
#Modified Date:29-09-2025
	#29-09-2025: included header in rend files for rms calculation and its components

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

# Reading data from the input files
for i in range (0,41):
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname)
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5
    
    # Creating a file 
	opfname=opname+str(i)+'.'+'csv'
	f=open(opfname,"w")
	f.write(f"Re,Rx2,Ry2,Rz2\n")

    # Calculating the end-to-end distance
	for c in range(36):
		xs=pos.at[20*c,'x']
		ys=pos.at[20*c,'y']
		zs=pos.at[20*c,'z']
		print(xs,ys,zs)
		xl=pos.at[20*c+19,'x']
		yl=pos.at[20*c+19,'y']
		zl=pos.at[20*c+19,'z']
		print(xl,yl,zl)

		xdiff=xl-xs
		ydiff=yl-ys
		zdiff=zl-zs

		sq=math.sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)
		Rx2=xdiff*xdiff
		Ry2=ydiff*ydiff
		Rz2=zdiff*zdiff

		# writing the end-to-end distance to the file 
		f.write(f"{sq},{Rx2},{Ry2},{Rz2}\n")
		print(sq)
	print("Exit")
	f.close();
