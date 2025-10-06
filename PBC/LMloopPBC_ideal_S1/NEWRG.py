# Computes radius of gyration (Rg) for free chains
# Translation of +0.5 is applied only locally (Alternative 3)
# Box size: 7x7x20  (lattice sites 0–6, 0–6, 0–19)
# Requires: mchains{i}.csv (bead coords), com_{i}.csv (center of mass)
# Modified Date: 06 Oct 2025  (Alternative 3 version)

import math
import pandas as pd
import os
import array

cwd = os.getcwd()
print("Working directory:", cwd)

ipname = 'mchains'
icom = 'com'
opname = 'rg'

# lattice box parameters
lx, ly, lz = 7, 7, 20  # number of sites
xms, yms, zms = float(lx)-1, float(ly)-1, float(lz)-1

# periodic boundary condition helper (minimum image)
def apply_pbc(diff, box_length):
    return diff - box_length * round(diff / box_length)

for i in range(0, 41):
    # Initialize arrays
    sxs = array.array('d', (0.0 for _ in range(36)))
    sys = array.array('d', (0.0 for _ in range(36)))
    szs = array.array('d', (0.0 for _ in range(36)))
    crg = array.array('d', (0.0 for _ in range(36)))
    Rgx2 = array.array('d', (0.0 for _ in range(36)))
    Rgy2 = array.array('d', (0.0 for _ in range(36)))
    Rgz2 = array.array('d', (0.0 for _ in range(36)))

    # Read COM file
    infcom = f"{icom}_{i}.csv"
    com = pd.read_csv(infcom, dtype={"c": int, "x": float, "y": float, "z": float})
    comx = com['x']
    comy = com['y']
    comz = com['z']

    # Read bead positions
    infname = f"{ipname}{i}.csv"
    pos = pd.read_csv(infname, dtype={"x": float, "y": float, "z": float})

    # Output file
    opfname = f"{opname}{i}.csv"
    f = open(opfname, "w")
    f.write("Rg,Rgx2,Rgy2,Rgz2\n")

    # Compute Rg for each chain
    for c in range(36):
        for b in range(20):
            # Local +0.5 translation (to lattice center) applied here only
            bead_x = pos.iloc[c*20 + b, 0] + 0.5
            bead_y = pos.iloc[c*20 + b, 1] + 0.5
            bead_z = pos.iloc[c*20 + b, 2] + 0.5

            com_x = comx.iloc[c]
            com_y = comy.iloc[c]
            com_z = comz.iloc[c]

            # distance between bead and COM
            xdiff = com_x - bead_x
            ydiff = com_y - bead_y
            zdiff = com_z - bead_z

            # Apply periodic boundary conditions (minimum image convention)
            xdiff = apply_pbc(xdiff, xms)
            ydiff = apply_pbc(ydiff, yms)
            zdiff = apply_pbc(zdiff, zms)

            # accumulate squared distances
            sxs[c] += xdiff * xdiff
            sys[c] += ydiff * ydiff
            szs[c] += zdiff * zdiff

        # Debug print (optional)
        print(f"Chain {c}: sumX={sxs[c]:.4f}, sumY={sys[c]:.4f}, sumZ={szs[c]:.4f}")

    # Compute average and write results
    for c in range(36):
        Rgx2[c] = sxs[c] / 20
        Rgy2[c] = sys[c] / 20
        Rgz2[c] = szs[c] / 20
        crg[c] = math.sqrt(Rgx2[c] + Rgy2[c] + Rgz2[c])

        f.write(f"{crg[c]},{Rgx2[c]},{Rgy2[c]},{Rgz2[c]}\n")

    print(f"Frame {i} processed")
    f.close()
