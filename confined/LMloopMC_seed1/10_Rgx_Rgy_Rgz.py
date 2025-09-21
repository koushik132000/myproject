#writing code to calculate the ensemble average of Rgx^2, Rgy^2, Rgz^2 without using numpy 
import math
import array
import os
import pandas as pd

# setting current directory
cwd = os.getcwd()
print(cwd)

#creating a file 
f = open("Rgxyz.csv","w")
f.write("Rgx,Rgy,Rgz,diff\n")

#Running through each file
for i in range(41):
    # creating arrays to store Rgxyz calculations
    sxs = array.array('d',(0.0 for c in range(36)))
    sys = array.array('d',(0.0 for c in range(36)))
    szs = array.array('d',(0.0 for c in range(36)))

    # Read the file
    pos = pd.read_csv(f"rg{i}.csv",dtype = {"c":int,"Rg":float,"Rgx":float,"Rgy":float,"Rgz":float})
    posx = pos[["Rgx"]] #extract Rgx series and keeps the dataframe structure for further calculations
    posy = pos[["Rgy"]]
    posz = pos[["Rgz"]]
    
    #Running through each chain and finding out the squares
    for c in range(36):
        sxs[c] = posx.iloc[c,0]**2
        sys[c] = posy.iloc[c,0]**2
        szs[c] = posz.iloc[c,0]**2
    
    #finding the sum
    sumx = 0
    sumy = 0
    sumz = 0
    for j in range(36):
        sumx += sxs[j];sumy += sys[j];sumz += szs[j]   

    #averaging over 36 chains and writing to the file
    Rgx = sumx/36
    Rgy = sumy/36
    Rgz = sumz/36
    diff = abs(Rgx - Rgy)
    f.write(f"{Rgx},{Rgy},{Rgz},{diff}\n")
f.close()
print("exit")