# creating a file to calculate the convergence of root mean square end to end distance
import pandas as pd
import numpy as np
import math
from pathlib import Path

cwd = Path.cwd()
print("1", cwd)
# creating a list to store values
result = []
SE = []

filename = "mderend"
for i in range(0,7):
    outfname = filename + str(i) + '.' + 'csv'
    df = pd.read_csv(outfname,header=None)
    #neglects the 1st value that is array index zero and starts from 2nd value till array index of 40
    z01 = df.iloc[1:41,0] #same for all the 7 files
    z0 = np.mean(z01)
    result.append(z0)

fmean = np.mean(result)
stderr = np.std(result, ddof=1) / np.sqrt(len(result))
print(f"Final Re = {fmean:.3f} ± {stderr:.3f}")
#df.to_csv("SE.csv", index=False, header=False)
df = pd.DataFrame([[fmean,stderr]],columns=["Mean_Re", "StdErr"] ) # converting SE into list 
df.to_csv("avgRe.csv", index=False)