# Computes mean, standard deviation and standard error for Radius of Gyration
# Requires rg data for chains
#Modified Date: 24Aug2025
		#10July2025: Added MPI parallelization
		#05Aug2025: Corrected input and output file names to include rank
        #21-08-2025: Corrected to Root mean square for Radius of Gyration.
        #24-08-2025: Corrected standard deviation and standard error 
import math
import pandas as pd
import os
import numpy as np
from string import *

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
#Identifies current working directory
cwd=os.getcwd()
print("1", cwd)

ipname='rg_rank'
f=open('mderg'+str(rank)+'.'+'csv',"w")

for i in range (0,41):
    infname = ipname+str(rank)+'_'+str(i)+'.csv'
    pos = pd.read_csv(infname, header=0)
    rg = pos.iloc[:, 1]   # Rg data
    npoints = np.size(rg)
    # Calculate <Rg^2>
    mean_sq = np.mean(rg**2, axis=0)
    Rg_rms = math.sqrt(mean_sq)# sqrt(<Rg^2>)

    # Standard deviation of Rg^2
    std_Rg2 = np.std(rg**2, axis=0)

    # Std of Rg (error propagation from variance)
    std_Rg = std_Rg2 /(2 * Rg_rms)

    # Standard error of <Rg^2>
    err_Rg2 = std_Rg2 / math.sqrt(npoints)

    # Error in Rg_rms (propagated)
    err_Rg = err_Rg2 /(2 * Rg_rms)
    print(Rg_rms,std_Rg,err_Rg)
    f.write(f"{Rg_rms:.6f},{std_Rg:.6f},{err_Rg:.6f}\n")
print("Exit")
f.close()
