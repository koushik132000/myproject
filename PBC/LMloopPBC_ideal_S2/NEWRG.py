# Computes radius of gyration after translating coordinates by +0.5
# Adaptive to wrapped/unwrapped configurations (consistent with COM logic)
# Box size: 7x7x20
# Modified: 08 Oct 2025
# Author: Koushik (M.Tech, IIT Hyderabad)

import math
import pandas as pd
import numpy as np
import os
import array

cwd = os.getcwd()
print("Working in", cwd)

ipname = 'mchains'
icom = 'com'
opname = 'rg'

# Box sizes (based on lattice sites)
xmax, ymax, zmax = 6.0, 6.0, 19.0  # for PBC difference (used in MIC)
xms, yms, zms = 7.0, 7.0, 20.0     # physical box lengths

nchains = 36
nbeads = 20
PI = math.pi


# ---------------------- Helper functions ----------------------

def apply_pbc(diff, box_length):
    """Applies minimum image convention to coordinate difference."""
    return diff - box_length * round(diff / box_length)

def circular_mean(values, box_length):
    """Computes circular mean for wrapped coordinates."""
    thetas = 2 * PI * (values / box_length)
    c_sum = np.mean(np.cos(thetas))
    s_sum = np.mean(np.sin(thetas))
    theta_avg = math.atan2(-s_sum, -c_sum) + PI
    return (box_length * theta_avg) / (2 * PI)

def is_wrapped(values, box_length):
    """Detects if data spans more than half the box — indicates wrapping."""
    return (max(values) - min(values)) > (box_length / 2.0)


# ---------------------- Main loop ----------------------

for i in range(0, 41):
    infcom = f"{icom}_{i}.csv"
    infname = f"{ipname}{i}.csv"

    # Read center of mass and chain positions
    com = pd.read_csv(infcom, dtype={"c": int, "x": float, "y": float, "z": float})
    pos = pd.read_csv(infname, dtype={"x": float, "y": float, "z": float})

    # Shift to cell centers
    xs = pos["x"] + 0.5
    ys = pos["y"] + 0.5
    zs = pos["z"] + 0.5

    # Output file
    opfname = f"{opname}{i}.csv"
    with open(opfname, "w") as f:
        f.write("Rg,Rgx2,Rgy2,Rgz2\n")

        # Initialize arrays
        sxs = np.zeros(nchains)
        sys = np.zeros(nchains)
        szs = np.zeros(nchains)

        # --- Per-chain calculation ---
        for c in range(nchains):
            # Extract block of 20 beads for each chain
            x_block = xs.iloc[c*nbeads:(c+1)*nbeads].to_numpy()
            y_block = ys.iloc[c*nbeads:(c+1)*nbeads].to_numpy()
            z_block = zs.iloc[c*nbeads:(c+1)*nbeads].to_numpy()

            cx = com.at[c, 'x']
            cy = com.at[c, 'y']
            cz = com.at[c, 'z']

            # For each bead in chain
            for b in range(nbeads):
                dx = x_block[b] - cx
                dy = y_block[b] - cy
                dz = z_block[b] - cz

                # Apply periodic wrapping (same as C code)
                dx = apply_pbc(dx, xms)
                dy = apply_pbc(dy, yms)

                # Decide if z is wrapped or not (like in COM logic)
                if is_wrapped(z_block, zms):
                    dz = apply_pbc(dz, zms)
                # else: keep dz as is (unwrapped)

                # Accumulate squared displacements
                sxs[c] += dx*dx
                sys[c] += dy*dy
                szs[c] += dz*dz

            # Mean square deviations → per-chain Rg
            Rgx2 = sxs[c] / nbeads
            Rgy2 = sys[c] / nbeads
            Rgz2 = szs[c] / nbeads
            Rg = math.sqrt(Rgx2 + Rgy2 + Rgz2)

            f.write(f"{Rg:.5f},{Rgx2:.5f},{Rgy2:.5f},{Rgz2:.5f}\n")

    print(f"Frame {i:02d} processed → {opfname}")
