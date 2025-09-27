#writing code to calculate the ensemble average of Rgx^2, Rgy^2, Rgz^2 without using numpy 
import math
import array
import os
import pandas as pd
import numpy as np

# setting current directory
cwd = os.getcwd()
print(cwd)

#creating a file 
f = open("Rgxyz.csv","w")
f.write("avg_Rgx2,avg_Rgy2,avg_Rgz2,SEx,SEy,SEz,avg_Rg2/3\n")

#Running through each file
for i in range(41):
    # creating arrays to store Rgxyz calculations
    sxs = np.zeros(36)
    sys = np.zeros(36)
    szs = np.zeros(36)

    # Read the file
    pos = pd.read_csv(f"rg{i}.csv",dtype = {"c":int,"Rg":float,"Rgx":float,"Rgy":float,"Rgz":float})
    posx = pos[["Rgx"]] #extract Rgx series and keeps the dataframe structure for further calculations
    posy = pos[["Rgy"]]
    posz = pos[["Rgz"]]

    # finding out the squares
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

    #averaging over 36 chains and finding standard error
    Rgx = sumx/36
    Rgy = sumy/36
    Rgz = sumz/36
    diff = abs(Rgx - Rgy)
    Rg = (Rgx + Rgy + Rgz)/3
    SEx = np.std(sxs,ddof =1)/math.sqrt(36)
    SEy = np.std(sys,ddof=1)/math.sqrt(36)
    SEz = np.std(szs,ddof=1)/math.sqrt(36)
    
    f.write(f"{Rgx},{Rgy},{Rgz},{SEx},{SEy},{SEz},{Rg}\n")
f.close()
print("exit")

# to find the average of standard errors
df = pd.read_csv("Rgxyz.csv")
avgSEx = df["avg_Rgx2"].std(ddof=1)/math.sqrt(41)
avgSEy = df["avg_Rgy2"].std(ddof=1)/math.sqrt(41)
avgSEz = df["avg_Rgz2"].std(ddof=1)/math.sqrt(41)
print("overall_SEX",avgSEx)
print("overall_SEY",avgSEy)
print("overall_SEZ",avgSEz)