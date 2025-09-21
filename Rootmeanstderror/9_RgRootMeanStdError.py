# Computes mean, standard deviation and standard error for Radius of Gyration
# Requires rg data for chains and finding one average over 36 chains
#Modified Date:24Aug2025
	#21-08-2025: Corrected to Root mean square for radius of gyration.
	#24-08-2025: Corrected standard deviation and standard error (Rg_rms)

#Load math, pandas for reading csv, string, os, numpy
import math
import pandas as pd
import os
import numpy as np

cwd = os.getcwd()
print("1", cwd)

ipname = 'rg'
f = open('mderg.csv', "w")

for i in range(0, 41):
    # Read rg data 
    infname = ipname + str(i) + '.csv'
    pos = pd.read_csv(infname, header=0)
    rg = pos.iloc[:, 1] # extract Rg column as a series
    size = np.size(rg)

    # Calculate Root-mean-square of Rg
    mean_sq = np.mean(rg**2, axis=0) #ensemble average of Rg^2
    Rg_rms = math.sqrt(mean_sq)

    # calculate standard deviation and standard error of Rg_rms
    std_Rg2 = np.std(rg**2, axis=0)
    std_Rg = std_Rg2 / (2 * Rg_rms)
    err_Rg2 = std_Rg2 / math.sqrt(size)
    err_Rg = err_Rg2 /(2 * Rg_rms)
    print(Rg_rms,std_Rg,err_Rg)
    f.write(f"{Rg_rms:.6f},{std_Rg:.6f},{err_Rg:.6f}\n")
print("Exit")
f.close()
