# Program provides lattice moves of 36 chains in 7x7x20 lattice
# Creation Date:27Nov2024
# Modified Date:01Dec2024
# 01-12-2024: Introduced subroutine to determine state of a lattice site

import math 
import random
import time
import csv

# Constants
nc = 36  # number of chains
bc = 20  # number of beads per chain
lx = 7   # box size in x-direction
ly = 7   # box size in y-direction
lz = 20  # box size in z-direction
ns = 980 # number of sites

# Structures
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

# Chain beads position for nc chains
beads = [Pos() for _ in range(nc)]

# Lattice configuration and state for 7x7x20 sites
lcs = Site() # lcs is a site tracker object and be used to call a site (lcs.sx[k])

# Defining unit vectors for the directions in the cubic lattice
dir = [Vec(1,0,0), Vec(-1,0,0), Vec(0,1,0), Vec(0,-1,0), Vec(0,0,1), Vec(0,0,-1)]
fpos = Vec() # fpos is a tracker for the first bead position (fpos.ex, fpos.ey, fpos.ez)
npos = Vec() # npos is a tracker for the last bead position (npos.ex, npos.ey, npos.ez)
kpos = Vec() # kpos is a tracker for the kth bead position( kpos.ex, kpos.ey, kpos.ez)

nb = bc-1  # index of the last bead
pb = bc-2  # index of the penultimate bead
cb = bc-3  # index of bead connected to penultimate bead

# Determines lattice site details
def sstate():
    # Initialize lattice site status and details
    for k in range(ns):
        lcs.sx[k] = 0
        lcs.sy[k] = 0
        lcs.sz[k] = 0

    # Update occupation status and lattice site details
    for i in range(nc):
        for j in range(bc):
            xval = beads[i].x[j]
            yval = beads[i].y[j]
            zval = beads[i].z[j]
            kidx = xval + lx * yval + lx * ly * zval # Calculating the site index, (2,3,5) gives a value of
        # 268 for 0 based indexing(based on given input), and 269 for 1 based indexing(counting in reality)
            lcs.sx[kidx] = i
            lcs.sy[kidx] = j
            lcs.sz[kidx] = 1

# First bead of a chain is moved
def fmoves(cnum, r): #cnum is the chain number, r is the direction index
    global fpos
    # Determining the current orientation of the first segment
    fvec = Vec() #fvec points from the first bead to the second bead
    # gives the difference in x,y,z coordinates between the second and first bead of a chain
    fvec.ex = beads[cnum].x[1] - beads[cnum].x[0] 
    fvec.ey = beads[cnum].y[1] - beads[cnum].y[0]
    fvec.ez = beads[cnum].z[1] - beads[cnum].z[0]
# fvec.ex, fvec.ey, fvec.ez tells in which the chain is oriented, i.e. the direction of the first bead
    # Determining if the bead should be moved
    if (fvec.ex == dir[r].ex and fvec.ey == dir[r].ey and fvec.ez == dir[r].ez):#if the first bead is 
        # already in the direction of the move, then do not move it 
        move = 0
    else:
        move = 1

    # Preventing folding back
    fvec.ex = beads[cnum].x[2] - beads[cnum].x[1]
    fvec.ey = beads[cnum].y[2] - beads[cnum].y[1]
    fvec.ez = beads[cnum].z[2] - beads[cnum].z[1]
    if (fvec.ex == -dir[r].ex and fvec.ey == -dir[r].ey and fvec.ez == -dir[r].ez):
        move = 0
    else:
        move = 1
 
    # Moving the first bead
    if move == 1:
        fpos.ex = beads[cnum].x[1] - dir[r].ex
        fpos.ey = beads[cnum].y[1] - dir[r].ey
        fpos.ez = beads[cnum].z[1] - dir[r].ez
    else:# this is a temporary overlap later resolved by occupancy check
        fpos.ex = beads[cnum].x[1]
        fpos.ey = beads[cnum].y[1]
        fpos.ez = beads[cnum].z[1]

    # Constraining the bead to the simulation box,like if (-1,2,5) is moved then it resets to original position
    if (fpos.ex < 0 or fpos.ey < 0 or fpos.ez < 0 or
        fpos.ex > lx-1 or fpos.ey > ly-1 or fpos.ez > lz-1):
        reset = 1
    else:
        reset = 0

    if reset == 1:
        fpos.ex = beads[cnum].x[0]
        fpos.ey = beads[cnum].y[0]
        fpos.ez = beads[cnum].z[0]

# Last bead of a chain is moved
def nmoves(cnum, r):
    global npos
    nvec = Vec()# determining the current orientation of the last segment
    nvec.ex = beads[cnum].x[nb] - beads[cnum].x[pb]
    nvec.ey = beads[cnum].y[nb] - beads[cnum].y[pb]
    nvec.ez = beads[cnum].z[nb] - beads[cnum].z[pb]

    if (nvec.ex == dir[r].ex and nvec.ey == dir[r].ey and nvec.ez == dir[r].ez):
        move = 0
    else:
        move = 1

    nvec.ex = beads[cnum].x[pb] - beads[cnum].x[cb]
    nvec.ey = beads[cnum].y[pb] - beads[cnum].y[cb]
    nvec.ez = beads[cnum].z[pb] - beads[cnum].z[cb]
    
    if (nvec.ex == -dir[r].ex and nvec.ey == -dir[r].ey and nvec.ez == -dir[r].ez):
        move = 0
    else:
        move = 1

    if move == 1:
        npos.ex = beads[cnum].x[pb] + dir[r].ex
        npos.ey = beads[cnum].y[pb] + dir[r].ey
        npos.ez = beads[cnum].z[pb] + dir[r].ez
    else:
        npos.ex = beads[cnum].x[pb]
        npos.ey = beads[cnum].y[pb]
        npos.ez = beads[cnum].z[pb]

    if (npos.ex < 0 or npos.ey < 0 or npos.ez < 0 or
        npos.ex > lx-1 or npos.ey > ly-1 or npos.ez > lz-1):
        reset = 1
    else:
        reset = 0

    if reset == 1:
        npos.ex = beads[cnum].x[nb]
        npos.ey = beads[cnum].y[nb]
        npos.ez = beads[cnum].z[nb]

# kth non-terminal bead of a chain is moved
def kmoves(cnum, k):
    global kpos
    # Determining the current orientation of the kth segment
    vec1 = Vec()
    vec2 = Vec()
    vec1.ex = beads[cnum].x[k+1] - beads[cnum].x[k]
    vec1.ey = beads[cnum].y[k+1] - beads[cnum].y[k]
    vec1.ez = beads[cnum].z[k+1] - beads[cnum].z[k]

    vec2.ex = beads[cnum].x[k] - beads[cnum].x[k-1]
    vec2.ey = beads[cnum].y[k] - beads[cnum].y[k-1]
    vec2.ez = beads[cnum].z[k] - beads[cnum].z[k-1]
# to check if the bead is perpendicular to the previous and next bead
    print(f"Vectors associated with the {k} bead")
    print(f"{vec1.ex}, {vec1.ey}, {vec1.ez}")
    print(f"{vec2.ex}, {vec2.ey}, {vec2.ez}")
    dota = vec1.ex * vec2.ex + vec1.ey * vec2.ey + vec1.ez * vec2.ez
    if dota != 0:
        move = 0
    else:
        move = 1

    if move == 1:
        kpos.ex = beads[cnum].x[k-1] + vec1.ex
        kpos.ey = beads[cnum].y[k-1] + vec1.ey
        kpos.ez = beads[cnum].z[k-1] + vec1.ez
    else:
        kpos.ex = beads[cnum].x[k]
        kpos.ey = beads[cnum].y[k]
        kpos.ez = beads[cnum].z[k]

    if (kpos.ex < 0 or kpos.ey < 0 or kpos.ez < 0 or
        kpos.ex > lx-1 or kpos.ey > ly-1 or kpos.ez > lz-1):
        reset = 1
    else:
        reset = 0

    if reset == 1:
        kpos.ex = beads[cnum].x[k]
        kpos.ey = beads[cnum].y[k]
        kpos.ez = beads[cnum].z[k]

# Main simulation
def main():
    global beads, lcs
    ccount = 0
    bcount = 0
    read = 0
    dindex = 0
    sindex = 0
    scalc = 0
    pcalc = 0
    cnum = 0
    findex = 0
    lindex = 0
    fbmove = 0
    nbmove = 0

    # Reading data from the input file
    try:
        with open("chains.csv", "r") as fptr:
            reader = csv.reader(fptr)
            header = next(reader)
            for row in reader:
                if len(row) != 3:
                    print(f"read error {len(row)}")
                    continue
                beads[ccount].x[bcount] = int(row[0])
                beads[ccount].y[bcount] = int(row[1])
                beads[ccount].z[bcount] = int(row[2])
                bcount += 1
                if bcount == bc:
                    ccount += 1
                    bcount = 0
    except Exception as e:
        print("Error reading file", e)
        return 1

    random.seed(time.time())

    sstate()
    # Determining lattice state in the simulation box
    for sindex in range(ns):
        if lcs.sz[sindex] == 1:
            print(f"site: {sindex} {lcs.sx[sindex]} {lcs.sy[sindex]} {lcs.sz[sindex]}")
        else:
            print(f"site: {sindex} not occupied")

    # Randomly choose a chain
    cnum = random.randint(0, nc-1)
    print(f"Chosen chain {cnum}")

    print("Initial position of the second bead ")
    print(f"{beads[cnum].x[1]}, {beads[cnum].y[1]}, {beads[cnum].z[1]}")
    print("Initial position of the penultimate bead ")
    print(f"{beads[cnum].x[pb]}, {beads[cnum].y[pb]}, {beads[cnum].z[pb]}")

    fbmove = 0
    nbmove = 0
    # Moving beads of the chosen chain
    for dindex in range(6):
        print(f"------Iteration {dindex} starts------")
        fmoves(cnum, dindex)
        print("Move first bead to new position")
        print(f"{fpos.ex}, {fpos.ey}, {fpos.ez}")
        scalc = fpos.ex + lx * fpos.ey + lx * ly * fpos.ez
        if lcs.sz[scalc] == 1:
            print("First bead cannot move to this site")
            print(f"It is occupied by bead {lcs.sy[scalc]} of chain {lcs.sx[scalc]}")
        else:
            print("First bead can be moved to this site")
            findex = dindex
            fbmove = 1

        nmoves(cnum, dindex)
        print("Move last bead to new position")
        print(f"{npos.ex}, {npos.ey}, {npos.ez}")
        scalc = npos.ex + lx * npos.ey + lx * ly * npos.ez
        if lcs.sz[scalc] == 1:
            print("Last bead cannot move to this site")
            print(f"It is occupied by bead {lcs.sy[scalc]} of chain {lcs.sx[scalc]}")
        else:
            print("Last bead can be moved to this site")
            lindex = dindex
            nbmove = 1
        print(f"------Iteration {dindex} ends------")

    if fbmove == 1:
        fmoves(cnum, findex)
        pcalc = beads[cnum].x[0] + lx * beads[cnum].y[0] + lx * ly * beads[cnum].z[0]
        lcs.sx[pcalc] = 0
        lcs.sy[pcalc] = 0
        lcs.sz[pcalc] = 0
        beads[cnum].x[0] = fpos.ex
        beads[cnum].y[0] = fpos.ey
        beads[cnum].z[0] = fpos.ez
        lcs.sx[scalc] = cnum
        lcs.sy[scalc] = 0
        lcs.sz[scalc] = 1

    if nbmove == 1:
        nmoves(cnum, lindex)
        pcalc = beads[cnum].x[nb] + lx * beads[cnum].y[nb] + lx * ly * beads[cnum].z[nb]
        lcs.sx[pcalc] = 0
        lcs.sy[pcalc] = 0
        lcs.sz[pcalc] = 0
        beads[cnum].x[nb] = npos.ex
        beads[cnum].y[nb] = npos.ey
        beads[cnum].z[nb] = npos.ez
        lcs.sx[scalc] = cnum
        lcs.sy[scalc] = nb
        lcs.sz[scalc] = 1

    # Move penultimate bead
    kmoves(cnum, pb)
    print(f"Move bead {pb} to new position")
    print(f"{kpos.ex}, {kpos.ey}, {kpos.ez}")
    scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
    if lcs.sz[scalc] == 1:
        print("Penultimate bead cannot move to this site")
        print(f"It is occupied by bead {lcs.sy[scalc]} of chain {lcs.sx[scalc]}")
    else:
        print("Penultimate bead can be moved to this site")
        pcalc = beads[cnum].x[pb] + lx * beads[cnum].y[pb] + lx * ly * beads[cnum].z[pb]
        lcs.sx[pcalc] = 0
        lcs.sy[pcalc] = 0
        lcs.sz[pcalc] = 0
        beads[cnum].x[pb] = kpos.ex
        beads[cnum].y[pb] = kpos.ey
        beads[cnum].z[pb] = kpos.ez
        lcs.sx[scalc] = cnum
        lcs.sy[scalc] = pb
        lcs.sz[scalc] = 1

    # Move second bead
    kmoves(cnum, 1)
    print(f"Move bead {2} to new position")
    print(f"{kpos.ex}, {kpos.ey}, {kpos.ez}")
    scalc = kpos.ex + lx * kpos.ey + lx * ly * kpos.ez
    if lcs.sz[scalc] == 1:
        print("Second bead cannot move to this site")
        print(f"It is occupied by bead {lcs.sy[scalc]} of chain {lcs.sx[scalc]}")
    else:
        print("Second bead can be moved to this site")
        pcalc = beads[cnum].x[1] + lx * beads[cnum].y[1] + lx * ly * beads[cnum].z[1]
        lcs.sx[pcalc] = 0
        lcs.sy[pcalc] = 0
        lcs.sz[pcalc] = 0
        beads[cnum].x[1] = kpos.ex
        beads[cnum].y[1] = kpos.ey
        beads[cnum].z[1] = kpos.ez
        lcs.sx[scalc] = cnum
        lcs.sy[scalc] = 1
        lcs.sz[scalc] = 1

    with open("mchains.csv", "w", newline='') as fptw:
        writer = csv.writer(fptw)
        writer.writerow(["x", "y", "z"])
        for i in range(nc):
            for j in range(bc):
                writer.writerow([beads[i].x[j], beads[i].y[j], beads[i].z[j]])

    return 0

if __name__ == "__main__":
    main()
