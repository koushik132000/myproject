# Computes center of mass (COM) with automatic detection for wrapped/unwrapped coordinates
# Analysis for free chains (36 chains, 7x7x20 box)
# Modified: 08 Oct 2025

import math
import pandas as pd
import os
import array
import numpy as np

cwd = os.getcwd()
print("Working directory:", cwd)

ipname = 'mchains'
opname = 'com'

# Box parameters
xms, yms, zms = 7.0, 7.0, 20.0  # box lengths
PI = math.pi
nchains = 36
nbeads = 20

# --- Helper: compute circular mean (angle mapping)
def circular_mean(values, box_length):
    thetas = 2 * PI * (values / box_length)
    c_sum = np.mean(np.cos(thetas))
    s_sum = np.mean(np.sin(thetas))
    theta_avg = math.atan2(-s_sum, -c_sum) + PI
    return (box_length * theta_avg) / (2 * PI)

# --- Helper: check if wrapping is needed
def is_wrapped(values, box_length):
    span = max(values) - min(values)
    return span > (box_length / 2.0)

# --- Loop over frames
for i in range(0, 41):
    infname = f"{ipname}{i}.csv"
    pos = pd.read_csv(infname, dtype={"x": float, "y": float, "z": float})

    # Shift to cell centers
    xs = pos["x"] + 0.5
    ys = pos["y"] + 0.5
    zs = pos["z"] + 0.5

    # Output file
    opfname = f"{opname}_{i}.csv"
    with open(opfname, "w") as f:
        f.write("c,x,y,z\n")

        # Compute COM per chain
        for c in range(nchains):
            x_block = xs.iloc[c*nbeads:(c+1)*nbeads].to_numpy()
            y_block = ys.iloc[c*nbeads:(c+1)*nbeads].to_numpy()
            z_block = zs.iloc[c*nbeads:(c+1)*nbeads].to_numpy()

            # X, Y are always periodic → use angle mapping
            cx = circular_mean(x_block, xms)
            cy = circular_mean(y_block, yms)
            cz = circular_mean(z_block, zms)

            f.write(f"{c},{cx:.5f},{cy:.5f},{cz:.5f}\n")

    print(f"COM file written: {opfname}")
