import numpy as np
import json
import os

# Load a sample field file
f = np.load("results/fields/se_n2_w0.00_Nw1_Re400_N384_tol1e-08.npz")
print("Keys in .npz file:", list(f.keys()))
print("\nShapes:")
for key in f.keys():
    print(f"  {key}: {f[key].shape}")
print("\nSample values (u[0,0]):", f['u'][0,0])
print("Grid size:", f['u'].shape)

# Load a sample scalar file
with open("results/scalars/se_n2_w0.00_Nw1_Re400_N384_tol1e-08.json") as fp:
    scalars = json.load(fp)
print("\nScalar keys:", list(scalars.keys()))
print("D*:", scalars.get('D_star'))
print("P*:", scalars.get('P_star'))
print("balance:", scalars.get('balance'))
