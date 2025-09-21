import numpy as np
import math
import time
import random
import csv

# Constants
nc = 36      # number of chains
bc = 20      # number of beads per chain
lx = 7       # box size in x-direction
ly = 7       # box size in y-direction
lz = 20      # box size in z-direction
ns = 980     # number of sites

nb = bc - 1  # index of last bead
pb = bc - 2  # index of penultimate bead
cb = bc - 3  # index of bead connected to penultimate bead

Eb = -0.457  # bead attraction energy
Ex = 200.0   # bead overlap energy

# Structures
class Vec:
    def __init__(self, ex=0, ey=0, ez=0):
        self.ex = ex
        self.ey = ey
        self.ez = ez

class Pos:
    def __init__(self):
        self.x = [0]*bc
        self.y = [0]*bc
        self.z = [0]*bc

class Site:
    def __init__(self):
        self.sx = [0]*ns
        self.sy = [0]*ns
        self.sz = [0]*ns

# Chain beads position for nc chains
beads = [Pos() for _ in range(nc)]

# Lattice configuration and state for 7x7x20 sites
lcs = Site()

# Defining unit vectors for the directions in the cubic lattice
dir = [Vec(1,0,0), Vec(-1,0,0), Vec(0,1,0), Vec(0,-1,0), Vec(0,0,1), Vec(0,0,-1)]
fpos = Vec()
npos = Vec()
kpos = Vec()

scalc = 0
pcalc = 0
fcalc = 0
lcalc = 0
dEf = 0.0
dEl = 0.0
dE = 0.0

# Random seed
random.seed(int(time.time()))
np.random.seed(int(time.time()))
seed = -1 - int(time.time())

def ran2():
    # Simple wrapper for random.random() for demonstration
    return random.random()

# Determines lattice site details
def sstate():
    for k in range(ns):
        lcs.sx[k] = 0
        lcs.sy[k] = 0
        lcs.sz[k] = 0
    for i in range(nc):
        for j in range(bc):
            xval = beads[i].x[j]
            yval = beads[i].y[j]
            zval = beads[i].z[j]
            kidx = xval + lx * yval + lx * ly * zval
            lcs.sx[kidx] = i
            lcs.sy[kidx] = j
            lcs.sz[kidx] += 1

def fmoves(cnum, r):
    global fpos
    fvec = Vec(
        beads[cnum].x[1] - beads[cnum].x[0],
        beads[cnum].y[1] - beads[cnum].y[0],
        beads[cnum].z[1] - beads[cnum].z[0]
    )
    move = 0 if (fvec.ex == dir[r].ex and fvec.ey == dir[r].ey and fvec.ez == dir[r].ez) else 1
    if move == 1:
        fpos.ex = beads[cnum].x[1] - dir[r].ex
        fpos.ey = beads[cnum].y[1] - dir[r].ey
        fpos.ez = beads[cnum].z[1] - dir[r].ez
    else:
        fpos.ex = beads[cnum].x[0]
        fpos.ey = beads[cnum].y[0]
        fpos.ez = beads[cnum].z[0]
    # Periodic boundary condition
    fpos.ex = fpos.ex % lx
    fpos.ey = fpos.ey % ly
    fpos.ez = fpos.ez % lz

def nmoves(cnum, r):
    global npos
    nvec = Vec(
        beads[cnum].x[nb] - beads[cnum].x[pb],
        beads[cnum].y[nb] - beads[cnum].y[pb],
        beads[cnum].z[nb] - beads[cnum].z[pb]
    )
    move = 0 if (nvec.ex == dir[r].ex and nvec.ey == dir[r].ey and nvec.ez == dir[r].ez) else 1
    if move == 1:
        npos.ex = beads[cnum].x[pb] + dir[r].ex
        npos.ey = beads[cnum].y[pb] + dir[r].ey
        npos.ez = beads[cnum].z[pb] + dir[r].ez
    else:
        npos.ex = beads[cnum].x[nb]
        npos.ey = beads[cnum].y[nb]
        npos.ez = beads[cnum].z[nb]
    # Periodic boundary condition
    npos.ex = npos.ex % lx
    npos.ey = npos.ey % ly
    npos.ez = npos.ez % lz

def kmoves(cnum, k):
    global kpos
    vec1 = Vec(
        beads[cnum].x[k+1] - beads[cnum].x[k],
        beads[cnum].y[k+1] - beads[cnum].y[k],
        beads[cnum].z[k+1] - beads[cnum].z[k]
    )
    vec2 = Vec(
        beads[cnum].x[k] - beads[cnum].x[k-1],
        beads[cnum].y[k] - beads[cnum].y[k-1],
        beads[cnum].z[k] - beads[cnum].z[k-1]
    )
    dota = vec1.ex * vec2.ex + vec1.ey * vec2.ey + vec1.ez * vec2.ez
    move = 0 if dota != 0 else 1
    if move == 1:
        kpos.ex = beads[cnum].x[k-1] + vec1.ex
        kpos.ey = beads[cnum].y[k-1] + vec1.ey
        kpos.ez = beads[cnum].z[k-1] + vec1.ez
    else:
        kpos.ex = beads[cnum].x[k]
        kpos.ey = beads[cnum].y[k]
        kpos.ez = beads[cnum].z[k]
    # Periodic boundary condition
    kpos.ex = kpos.ex % lx
    kpos.ey = kpos.ey % ly
    kpos.ez = kpos.ez % lz

def deltaE(olds, news):
    olist = [0]*6
    nlist = [0]*6
    Ediff = 0.0
    if lcs.sz[news] >= 1:
        return lcs.sz[news] * Ex
    # Neighbors of the old site
    olist[0] = olds+1 if (olds+1)%7 != 0 else olds-(lx-1)
    olist[1] = olds-1 if olds%7 != 0 else olds+(lx-1)
    olist[2] = olds+lx if (olds+1)%49 != 0 else olds-lx*(ly-1)
    olist[3] = olds-lx if olds%49 != 0 else olds+lx*(ly-1)
    olist[4] = olds+lx*ly if olds < 931 else olds-lx*ly*(lz-1)
    olist[5] = olds-lx*ly if olds > 48 else olds+lx*ly*(lz-1)
    # Neighbors of the new site
    nlist[0] = news+1 if (news+1)%7 != 0 else news-(lx-1)
    nlist[1] = news-1 if news%7 != 0 else news+(lx-1)
    nlist[2] = news+lx if (news+1)%49 != 0 else news-lx*(ly-1)
    nlist[3] = news-lx if news%49 != 0 else news+lx*(ly-1)
    nlist[4] = news+lx*ly if news < 931 else news-lx*ly*(lz-1)
    nlist[5] = news-lx*ly if news > 48 else news+lx*ly*(lz-1)
    oldn = sum(lcs.sz[i] for i in olist if 0 <= i < 980)
    newn = sum(lcs.sz[i] for i in nlist if 0 <= i < 980)
    Ediff = (newn - oldn) * Eb
    return Ediff

def fmeval(cnum, dindex, fcalc):
    global dEf
    fmoves(cnum, dindex)
    scalc = fpos.ex + lx * fpos.ey + lx * ly * fpos.ez
    dEf = deltaE(fcalc, scalc)

def lmeval(cnum, dindex, lcalc):
    global dEl
    nmoves(cnum, dindex)
    scalc = npos.ex + lx * npos.ey + lx * ly * npos.ez
    dEl = deltaE(lcalc, scalc)

def metrop(delE):
    if delE <= 0:
        pij = 1.0
        acc = 1
    else:
        pij = math.exp(-delE)
        rij = ran2()
        acc = 0 if pij < rij else 1
    return acc

def accmov(cnum, bnum, pcalc, scalc, mpos):
    lcs.sx[pcalc] = 0
    lcs.sy[pcalc] = 0
    lcs.sz[pcalc] -= 1
    beads[cnum].x[bnum] = mpos.ex
    beads[cnum].y[bnum] = mpos.ey
    beads[cnum].z[bnum] = mpos.ez
    lcs.sx[scalc] = cnum
    lcs.sy[scalc] = bnum
    lcs.sz[scalc] += 1

# Main simulation
def main():
    global dEf, dEl, dE, fcalc, lcalc, pcalc, scalc
    ccount = 0
    bcount = 0
    read = 0
    sucm = 0
    totm = 0
    maxm = int(100E6)
    fn = 0

    # Reading data from the input file
    with open("chains.csv", "r") as fptr:
        csvreader = csv.reader(fptr)
        header = next(csvreader)
        for row in csvreader:
            if len(row) == 3:
                beads[ccount].x[bcount] = int(row[0])
                beads[ccount].y[bcount] = int(row[1])
                beads[ccount].z[bcount] = int(row[2])
                bcount += 1
                if bcount == bc:
                    ccount += 1
                    bcount = 0

    # Determining lattice state in the simulation box
    sstate()
    for sindex in range(ns):
        if lcs.sz[sindex] == 1:
            print(f"site: {sindex} {lcs.sx[sindex]} {lcs.sy[sindex]} {lcs.sz[sindex]}")
        else:
            print(f"site: {sindex} not occupied")

    # Write initial configuration to a file
    with open(f"mchains{fn}.csv", "w", newline='') as fptw:
        writer = csv.writer(fptw)
        writer.writerow(["x", "y", "z"])
        for i in range(nc):
            for j in range(bc):
                writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

    while totm < maxm:
        cnum = random.randint(0, nc - 1)
        # Move the first bead
        dEf = 0.0
        fcalc = beads[cnum].x[0] + lx * beads[cnum].y[0] + lx * ly * beads[cnum].z[0]
        findex = random.randint(0, 5)
        fmeval(cnum, findex, fcalc)
        if metrop(dEf) == 1:
            scalc = fpos.ex + lx * fpos.ey + lx * ly * fpos.ez
            accmov(cnum, 0, fcalc, scalc, fpos)
            sucm += 1
        # Move last bead
        dEl = 0.0
        lcalc = beads[cnum].x[nb] + lx * beads[cnum].y[nb] + lx * ly * beads[cnum].z[nb]
        lindex = random.randint(0, 5)
        lmeval(cnum, lindex, lcalc)
        if metrop(dEl) == 1:
            scalc = npos.ex + lx * npos.ey + lx * ly * npos.ez
            accmov(cnum, nb, lcalc, scalc, npos)
            sucm += 1
        # Move penultimate bead
        dE = 0.0
        pcalc = beads[cnum].x[pb] + lx * beads[cnum].y[pb] + lx * ly * beads[cnum].z[pb]
        kmoves(cnum, pb)
        scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
        dE = deltaE(pcalc, scalc)
        if metrop(dE) == 1:
            accmov(cnum, pb, pcalc, scalc, kpos)
            sucm += 1
        # Move second bead
        dE = 0.0
        pcalc = beads[cnum].x[1] + lx * beads[cnum].y[1] + lx * ly * beads[cnum].z[1]
        kmoves(cnum, 1)
        scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
        dE = deltaE(pcalc, scalc)
        if metrop(dE) == 1:
            accmov(cnum, 1, pcalc, scalc, kpos)
            sucm += 1
        # Move random internal bead of chain
        while True:
            bnum = random.randint(0, bc - 1)
            if bnum != 0 and bnum != 19:
                break
        dE = 0.0
        pcalc = beads[cnum].x[bnum] + lx * beads[cnum].y[bnum] + lx * ly * beads[cnum].z[bnum]
        kmoves(cnum, bnum)
        scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
        dE = deltaE(pcalc, scalc)
        if metrop(dE) == 1:
            accmov(cnum, bnum, pcalc, scalc, kpos)
            sucm += 1
        totm += 1
        if totm % 5000000 == 0:
            print("I am here")
            fn += 1
            with open(f"mchains{fn}.csv", "w", newline='') as fptw:
                writer = csv.writer(fptw)
                writer.writerow(["x", "y", "z"])
                for i in range(nc):
                    for j in range(bc):
                        writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

    print(f"Total number of moves {totm} and successful moves {sucm}")
    # Check if a site has more than one bead
    sstate()
    for sindex in range(ns):
        if lcs.sz[sindex] > 1:
            print(f"site: {sindex} {lcs.sx[sindex]} {lcs.sy[sindex]} {lcs.sz[sindex]}")
    with open("mchainsf.csv", "w", newline='') as fptw:
        writer = csv.writer(fptw)
        writer.writerow(["x", "y", "z"])
        for i in range(nc):
            for j in range(bc):
                writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

if __name__ == "__main__":
    main()
