import math
import random
import time
import csv

# Constants
nc = 36   # number of chains
bc = 20   # number of beads per chain
lx = 7    # box size in x-direction
ly = 7    # box size in y-direction
lz = 20   # box size in z-direction
ns = 980  # number of sites

# Energy parameters
Ea = -2.232  # Bead-Wall interaction energy
Eb = -0.457  # Bead-Bead interaction energy
Ex = 200.0   # Bead overlap energy

nb = bc - 1
pb = bc - 2
cb = bc - 3

# Data structures
class Pos:
    def __init__(self):
        self.x = [0] * bc
        self.y = [0] * bc
        self.z = [0] * bc

class Vec:
    def __init__(self, ex=0, ey=0, ez=0):
        self.ex = ex
        self.ey = ey
        self.ez = ez

class Site:
    def __init__(self):
        self.sx = [0] * ns
        self.sy = [0] * ns
        self.sz = [0] * ns

# Initialize
beads = [Pos() for _ in range(nc)]
lcs = Site()
dir = [Vec(1,0,0), Vec(-1,0,0), Vec(0,1,0), Vec(0,-1,0), Vec(0,0,1), Vec(0,0,-1)]
fpos = Vec()
npos = Vec()
kpos = Vec()

# Random seed
random.seed(int(time.time()))
seed = int(-1 - time.time()) if time.time() >= 0 else int(time.time())

def ran2(seed):
    # Placeholder for the ran2 random number generator
    return random.random()

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
            k_idx = xval + lx * yval + lx * ly * zval
            lcs.sx[k_idx] = i
            lcs.sy[k_idx] = j
            lcs.sz[k_idx] += 1

def fmoves(cnum, r):
    global fpos
    fvec = Vec(
        beads[cnum].x[1] - beads[cnum].x[0],
        beads[cnum].y[1] - beads[cnum].y[0],
        beads[cnum].z[1] - beads[cnum].z[0]
    )
    move = not (fvec.ex == dir[r].ex and fvec.ey == dir[r].ey and fvec.ez == dir[r].ez)
    if move:
        fpos.ex = beads[cnum].x[1] - dir[r].ex
        fpos.ey = beads[cnum].y[1] - dir[r].ey
        fpos.ez = beads[cnum].z[1] - dir[r].ez
    else:
        fpos.ex = beads[cnum].x[0]
        fpos.ey = beads[cnum].y[0]
        fpos.ez = beads[cnum].z[0]

    # Periodic boundary condition
    fpos.ex %= lx
    fpos.ey %= ly

    # Wall boundary condition in z
    if fpos.ez < 0 or fpos.ez > lz - 1:
        fpos.ez = beads[cnum].z[0]

def nmoves(cnum, r):
    global npos
    nvec = Vec(
        beads[cnum].x[nb] - beads[cnum].x[pb],
        beads[cnum].y[nb] - beads[cnum].y[pb],
        beads[cnum].z[nb] - beads[cnum].z[pb]
    )
    move = not (nvec.ex == dir[r].ex and nvec.ey == dir[r].ey and nvec.ez == dir[r].ez)
    if move:
        npos.ex = beads[cnum].x[pb] + dir[r].ex
        npos.ey = beads[cnum].y[pb] + dir[r].ey
        npos.ez = beads[cnum].z[pb] + dir[r].ez
    else:
        npos.ex = beads[cnum].x[nb]
        npos.ey = beads[cnum].y[nb]
        npos.ez = beads[cnum].z[nb]

    # Periodic boundary condition
    npos.ex %= lx
    npos.ey %= ly

    # Wall boundary condition in z
    if npos.ez < 0 or npos.ez > lz - 1:
        npos.ez = beads[cnum].z[nb]

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
    move = (dota == 0)
    if move:
        kpos.ex = beads[cnum].x[k-1] + vec1.ex
        kpos.ey = beads[cnum].y[k-1] + vec1.ey
        kpos.ez = beads[cnum].z[k-1] + vec1.ez
    else:
        kpos.ex = beads[cnum].x[k]
        kpos.ey = beads[cnum].y[k]
        kpos.ez = beads[cnum].z[k]

    # Periodic boundary condition
    kpos.ex %= lx
    kpos.ey %= ly

    # Wall boundary condition in z
    if kpos.ez < 0 or kpos.ez > lz - 1:
        kpos.ez = beads[cnum].z[k]

def deltaE(olds, news):
    Ediff = 0.0
    if lcs.sz[news] >= 1:
        return lcs.sz[news] * Ex
    # Neighbor calculations
    def get_neighbors(idx):
        neighbors = [0]*6
        neighbors[0] = idx + 1 if (idx + 1) % 7 != 0 else idx - (lx - 1)
        neighbors[1] = idx - 1 if idx % 7 != 0 else idx + (lx - 1)
        neighbors[2] = idx + lx if (idx + 1) % 49 != 0 else idx - lx * (ly - 1)
        neighbors[3] = idx - lx if idx % 49 != 0 else idx + lx * (ly - 1)
        neighbors[4] = idx + lx * ly if idx < 931 else -1
        neighbors[5] = idx - lx * ly if idx > 48 else -1
        return neighbors

    olist = get_neighbors(olds)
    nlist = get_neighbors(news)
    oldn = sum(lcs.sz[o] for o in olist if 0 <= o < ns)
    newn = sum(lcs.sz[n] for n in nlist if 0 <= n < ns)
    Ediff = (newn - oldn) * Eb
    return Ediff

def fmeval(cnum, dindex, fcalc):
    fmoves(cnum, dindex)
    scalc = fpos.ex + lx * fpos.ey + lx * ly * fpos.ez
    return deltaE(fcalc, scalc)

def lmeval(cnum, dindex, lcalc):
    nmoves(cnum, dindex)
    scalc = npos.ex + lx * npos.ey + lx * ly * npos.ez
    return deltaE(lcalc, scalc)

def metrop(delE):# delE is a variable that stores the value returned by deltaE
    if delE <= 0:
        return True
    pij = math.exp(-delE)
    rij = ran2(seed)
    return pij >= rij

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
    maxm = int(1e8)
    mcstep = int(5e6)
    sucm = 0
    totm = 0
    fn = 0

    # Read chain data from file
    with open("chains.csv", "r") as fptr:
        reader = csv.reader(fptr)
        next(reader)  # Skip header
        ccount = 0
        bcount = 0
        for row in reader:
            if len(row) != 3:
                continue
            beads[ccount].x[bcount] = int(row[0])
            beads[ccount].y[bcount] = int(row[1])
            beads[ccount].z[bcount] = int(row[2])
            bcount += 1
            if bcount == bc:
                ccount += 1
                bcount = 0
            if ccount == nc:
                break

    # Initial lattice state
    sstate()

    # Write initial configuration
    with open(f"mchains{fn}.csv", "w", newline='') as fptw:
        writer = csv.writer(fptw)
        writer.writerow(["x", "y", "z"])
        for i in range(nc):
            for j in range(bc):
                writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

    while totm < maxm:
        cnum = random.randint(0, nc-1)
        # Move first bead
        fcalc = beads[cnum].x[0] + lx * beads[cnum].y[0] + lx * ly * beads[cnum].z[0]
        findex = random.randint(0, 5)
        dEf = fmeval(cnum, findex, fcalc)
        if metrop(dEf):
            scalc = fpos.ex + lx * fpos.ey + lx * ly * fpos.ez
            accmov(cnum, 0, fcalc, scalc, fpos)
            sucm += 1

        # Move last bead
        lcalc = beads[cnum].x[nb] + lx * beads[cnum].y[nb] + lx * ly * beads[cnum].z[nb]
        lindex = random.randint(0, 5)
        dEl = lmeval(cnum, lindex, lcalc)
        if metrop(dEl):
            scalc = npos.ex + lx * npos.ey + lx * ly * npos.ez
            accmov(cnum, nb, lcalc, scalc, npos)
            sucm += 1

        # Move penultimate bead
        pcalc = beads[cnum].x[pb] + lx * beads[cnum].y[pb] + lx * ly * beads[cnum].z[pb]
        kmoves(cnum, pb)
        scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
        dE = deltaE(pcalc, scalc)
        if metrop(dE):
            accmov(cnum, pb, pcalc, scalc, kpos)
            sucm += 1

        # Move second bead
        pcalc = beads[cnum].x[1] + lx * beads[cnum].y[1] + lx * ly * beads[cnum].z[1]
        kmoves(cnum, 1)
        scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
        dE = deltaE(pcalc, scalc)
        if metrop(dE):
            accmov(cnum, 1, pcalc, scalc, kpos)
            sucm += 1

        # Move random internal bead
        while True:
            bnum = random.randint(0, bc-1)
            if bnum != 0 and bnum != 19:
                break
        pcalc = beads[cnum].x[bnum] + lx * beads[cnum].y[bnum] + lx * ly * beads[cnum].z[bnum]
        kmoves(cnum, bnum)
        scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
        dE = deltaE(pcalc, scalc)
        if metrop(dE):
            accmov(cnum, bnum, pcalc, scalc, kpos)
            sucm += 1

        totm += 1
        if totm % mcstep == 0:
            fn += 1
            with open(f"mchains{fn}.csv", "w", newline='') as fptw:
                writer = csv.writer(fptw)
                writer.writerow(["x", "y", "z"])
                for i in range(nc):
                    for j in range(bc):
                        writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

    print(f"Total number of moves {totm} and successful moves {sucm}")

    # Final lattice state
    sstate()
    with open("mchainsf.csv", "w", newline='') as fptw:
        writer = csv.writer(fptw)
        writer.writerow(["x", "y", "z"])
        for i in range(nc):
            for j in range(bc):
                writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

if __name__ == "__main__":
    main()
