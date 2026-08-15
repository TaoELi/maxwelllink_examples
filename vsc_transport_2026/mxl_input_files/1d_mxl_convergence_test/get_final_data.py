import numpy as np
import h5py

repeat_list = [1, 3, 5, 10]
amp_list = ["0.001", "0.0016", "0.003", "0.005", "0.007"]
for y0 in range(5):
    for z0 in repeat_list:
        sp = 0
        for x0 in range(z0):
            with h5py.File(f"./data/multimode_cavmd_144_24_{amp_list[y0]}_{x0}_neq.h5", "r") as f:
                data = {key: f[key][:] for key in f.keys()}
                sp += np.sum(np.reshape(data["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2
        np.save(f"./scpdata/multimode_cavmd_144_24_{amp_list[y0]}_{z0}_neq_final.npy", sp / z0)
