import numpy as np
import h5py

amplitude = 0.002
efield_y_avg = 0

print(f'get data from 2d_144_144_efield_neq.h5')
with h5py.File(f"./2d_144_144_efield_neq.h5", "r") as f:
    efield_y_avg += f["effective_efield"][:, : ,1]**2
np.save(f"./2d_144_144_efield_neq_{amplitude}.npy", efield_y_avg)
print(f"shape {efield_y_avg.shape} in 2d_144_144_efield_neq_{amplitude}_neq.npy")