#Computes end-to-end distance after translating all coordinates by 0.5
#The analysis is for a code with periodic boundary conditions
#Box size is 7x7x20
#Requires modified bead coordinate file for 36 chains
#modified date: 28Aug2025
    #17-08-2025: corrected to proper minimum image convention for each bond vector
    #28-08-2025: Corrected to unwrap coordinates for end to end distance
    #27-09-2025: writing comments for clarity of pandas indexing, csv file indexing and adding header to rend files
#Load math, pandas for reading csv, string, os
import math
import pandas as pd
import os

#Identifies current working directory
cwd = os.getcwd()
print("Working in", cwd)

ipname = 'mchains'
opname = 'rend'

xmax, ymax, zmax = 7.0, 7.0, 20.0
nb = 20
nchains = 36

for i in range(0, 41): #i from 0 to 40
    infname=ipname+str(i)+'.'+'csv'
    pos = pd.read_csv(infname) # default removes header and row 2 as index 0 till index 719
    xs=pos[['x']]+0.5
    ys=pos[['y']]+0.5
    zs=pos[['z']]+0.5

    # creating a file
    opfname=opname+str(i)+'.'+'csv'
    f = open(opfname, "w")
    f.write(f"Re,Rx2,Ry2,Rz2\n")

    # Calculating the end-to-end distance
    for c in range(nchains): # c from 0 to 35
        # wrapped first bead for unwrapping reference
        x0 = pos.at[20*c, 'x']
        y0 = pos.at[20*c, 'y']
        z0 = pos.at[20*c, 'z']

        # start unwrapped coords at the first bead
        xu = x0
        yu = y0
        zu = z0

        # moving through the bonds, using MIC for each bond vector
        px = x0
        py = y0
        pz = z0 

        for k in range(1, nb): # k from 1 to 19
            xk = pos.at[20*c + k, 'x']
            yk = pos.at[20*c + k, 'y']
            zk = pos.at[20*c + k, 'z']

            dx = xk - px
            dy = yk - py
            dz = zk - pz

            # minimum image for each bond vector
            dx = dx - xmax * round(dx / xmax)
            dy = dy - ymax * round(dy / ymax)
            dz = dz - zmax * round(dz / zmax)

            # adding up to build unwrapped coordinates
            xu += dx
            yu += dy
            zu += dz
            px = xk; py = yk; pz = zk

        # end-to-end vector from unwrapped endpoints
        Rx = xu - x0
        Ry = yu - y0
        Rz = zu - z0

        sq = math.sqrt(Rx*Rx + Ry*Ry + Rz*Rz)
        Rx2 = Rx*Rx
        Ry2 = Ry*Ry
        Rz2 = Rz*Rz
        # writing the result to the file 
        f.write(f"{sq},{Rx2},{Ry2},{Rz2}\n")
        print(sq)
    print("Exit")
    f.close();

