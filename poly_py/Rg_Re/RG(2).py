#Computes radius of gyration after translating the coordinated by 0.5
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
icom='com'
opname='rg'

for i in range (0,21):
	#Intialize arrays of doubles
	sxs=array.array('d',(0.0 for c in range(0,36)))
	sys=array.array('d',(0.0 for c in range(0,36)))
	szs=array.array('d',(0.0 for c in range(0,36)))
	crg=array.array('d',(0.0 for c in range(0,36)))
	
	infcom=icom+str(i)+'.'+'csv'
	com=pd.read_csv(infcom,dtype={"c":int,"x":float,"y":float,"z":float})
	comx=com[['x']]
	comy=com[['y']]
	comz=com[['z']]
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,dtype={"x":float,"y":float,"z":float})
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5
	
	opfname=opname+str(i)+'.'+'csv'
	f=open(opfname,"w")
	f.write(f"c,Rg\n")
	for c in range(36):
		for b in range(20):
			sxs[c]+=(comx.iloc[c,0]-xs.iloc[c*20+b,0])*(comx.iloc[c,0]-xs.iloc[c*20+b,0])
			sys[c]+=(comy.iloc[c,0]-ys.iloc[c*20+b,0])*(comy.iloc[c,0]-ys.iloc[c*20+b,0])
			szs[c]+=(comz.iloc[c,0]-zs.iloc[c*20+b,0])*(comz.iloc[c,0]-zs.iloc[c*20+b,0])
		print(sxs[c],sys[c],szs[c])			
	for c in range(36):
		crg[c]=math.sqrt(sxs[c]/20+sys[c]/20+szs[c]/20)
		f.write(f"{c},{crg[c]}\n")
	print("Exit")
	f.close();
