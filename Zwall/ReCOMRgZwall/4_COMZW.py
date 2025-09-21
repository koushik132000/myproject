#Computes center of mass after translating the coordinated by 0.5
#The analysis is for chains with PBC in x and y and wall in the z direction
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains

#Load math, pandas for reading csv, string, os
import math
import pandas as pd
from string import *
import os
import array

#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='mchains'
opname='com'

xmax=7
ymax=7

# adjusting the box size for periodic boundary conditions 
xms=xmax+0.5
yms=ymax+0.5

PI=math.pi

for i in range (0,41):
	#Intialize arrays of doubles
	sxsc=array.array('d',(0.0 for c in range(0,36)))
	sxss=array.array('d',(0.0 for c in range(0,36)))
	sysc=array.array('d',(0.0 for c in range(0,36)))
	syss=array.array('d',(0.0 for c in range(0,36)))

	txsc=array.array('d',(0.0 for c in range(0,36)))
	txss=array.array('d',(0.0 for c in range(0,36)))
	tysc=array.array('d',(0.0 for c in range(0,36)))
	tyss=array.array('d',(0.0 for c in range(0,36)))
		
	szs=array.array('d',(0.0 for c in range(0,36)))
	
	cxs=array.array('d',(0.0 for c in range(0,36)))	
	cys=array.array('d',(0.0 for c in range(0,36)))	
	czs=array.array('d',(0.0 for c in range(0,36)))
	
	# Reading data from the input file and translating coordinates
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,dtype={"x":float,"y":float,"z":float})
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5
	
	# creating and write to a file
	opfname=opname+"_"+str(i)+'.'+'csv'
	f=open(opfname,"w")
	f.write(f"c,x,y,z\n")
	
	for c in range(36):
		for b in range(20):
			mxs=xs.iloc[c*20+b,0]
			mys=ys.iloc[c*20+b,0]
			mzs=zs.iloc[c*20+b,0]

			# convert coordinates to angles
			thetax=2*PI*(mxs/xms)	
			thetay=2*PI*(mys/yms)	
	
			# calculating the cosine and sine of the angles
			sxsc[c]+=math.cos(thetax)
			sysc[c]+=math.cos(thetay)
	
			sxss[c]+=math.sin(thetax)
			syss[c]+=math.sin(thetay)
	
			szs[c]+=mzs

		print(c,sxsc[c],sxss[c],sysc[c],syss[c],szs[c],"\n")		

	for c in range(36):
		# calculating the average values
		txsc[c]=sxsc[c]/20	
		tysc[c]=sysc[c]/20
		
		txss[c]=sxss[c]/20
		tyss[c]=syss[c]/20
  
		# calculating the average angles
		tavgx=math.atan2(-txss[c],-txsc[c])+PI	
		tavgy=math.atan2(-tyss[c],-tysc[c])+PI	
			
		# converting angles back to center of mass coordinates
		cxs[c]=round(xms*tavgx/(2*PI),5)
		cys[c]=round(yms*tavgy/(2*PI),5)
		czs[c]=szs[c]/20

		# writing the center of mass values to the file
		f.write(f"{c},{cxs[c]},{cys[c]},{czs[c]}\n")
	print("Exit")
	f.close();
