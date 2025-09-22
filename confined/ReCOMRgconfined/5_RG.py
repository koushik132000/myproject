#Computes radius of gyration after translating the coordinated by 0.5
#The analysis is for confined chains
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains
#Modified Date:15Sep2025
		#15-09-2025: Initialized array of doubles for Rgx,Rgy,Rgz and writing them to the file

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

for i in range (0,41):
	#Intialize arrays of doubles
	sxs=array.array('d',(0.0 for c in range(0,36)))
	sys=array.array('d',(0.0 for c in range(0,36)))
	szs=array.array('d',(0.0 for c in range(0,36)))
	crg=array.array('d',(0.0 for c in range(0,36)))
	Rgx=array.array('d',(0.0 for c in range(0,36)))
	Rgy=array.array('d',(0.0 for c in range(0,36)))
	Rgz=array.array('d',(0.0 for c in range(0,36)))
	
    # Reading the center of mass files 
	infcom=icom+"_"+str(i)+'.'+'csv'
	com=pd.read_csv(infcom,dtype={"c":int,"x":float,"y":float,"z":float})
	comx=com[['x']]
	comy=com[['y']]
	comz=com[['z']]
 
    # Reading the bead coordinates from the input file
	infname=ipname+str(i)+'.'+'csv'
	pos=pd.read_csv(infname,dtype={"x":float,"y":float,"z":float})
	xs=pos[['x']]+0.5
	ys=pos[['y']]+0.5
	zs=pos[['z']]+0.5
	
   # creating and writing to a file 
	opfname=opname+str(i)+'.'+'csv'
	f=open(opfname,"w")
	f.write(f"c,Rg,Rgx,Rgy,Rgz\n")
 
    # calculating the radius of gyration for each chain
	for c in range(36):
		for b in range(20):
			sxs[c]+=(comx.iloc[c,0]-xs.iloc[c*20+b,0])*(comx.iloc[c,0]-xs.iloc[c*20+b,0])
			sys[c]+=(comy.iloc[c,0]-ys.iloc[c*20+b,0])*(comy.iloc[c,0]-ys.iloc[c*20+b,0])
			szs[c]+=(comz.iloc[c,0]-zs.iloc[c*20+b,0])*(comz.iloc[c,0]-zs.iloc[c*20+b,0])
		print(sxs[c],sys[c],szs[c])

	for c in range(36):
		crg[c]=math.sqrt(sxs[c]/20+sys[c]/20+szs[c]/20)
		Rgx[c] = math.sqrt(sxs[c]/20)
		Rgy[c] = math.sqrt(sys[c]/20)
		Rgz[c] = math.sqrt(szs[c]/20)
		# writing the results to the output file
		f.write(f"{c},{crg[c]},{Rgx[c]},{Rgy[c]},{Rgz[c]}\n")
	print("Exit")
	f.close();
