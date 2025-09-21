#Computes center of mass after translating the coordinated by 0.5
#The analysis is for confined chains
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

for i in range (0,21):
	#Intialize arrays of doubles
	sxs=array.array('d',(0.0 for c in range(0,36)))
	sys=array.array('d',(0.0 for c in range(0,36)))
	szs=array.array('d',(0.0 for c in range(0,36)))
	cxs=array.array('d',(0.0 for c in range(0,36)))
	cys=array.array('d',(0.0 for c in range(0,36)))
	czs=array.array('d',(0.0 for c in range(0,36)))

	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,dtype={"x":float,"y":float,"z":float})
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5
	
	opfname=opname+str(i)+'.'+'csv'
	f=open(opfname,"w")
	f.write(f"c,x,y,z\n")
	
	for c in range(36):
		for b in range(20):
			sxs[c]+=xs.iloc[c*20+b,0]
			sys[c]+=ys.iloc[c*20+b,0]
			szs[c]+=zs.iloc[c*20+b,0]
		print(sxs[c]/20,sys[c]/20,szs[c]/20)			
	for c in range(36):
		cxs[c]=sxs[c]/20
		cys[c]=sys[c]/20
		czs[c]=szs[c]/20
		f.write(f"{c},{cxs[c]},{cys[c]},{czs[c]}\n")
	print("Exit")
	f.close();
