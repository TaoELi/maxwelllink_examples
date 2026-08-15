import numpy as np
import sys, h5py
from maxwelllink.tools import ir_spectrum

sp_avg = 0
for i in range(0,10):
    with h5py.File(f"./data/multimode_cavmd_144_{i}_eq.h5", "r") as f:
        qc_y = f["qc"][:,:,1]
        sp = np.zeros((qc_y.shape[0]//2, qc_y.shape[1]))
        for j in range(qc_y.shape[1]):
            x, spj = ir_spectrum(qc_y[:,j], 2)
            sp[:,j] = spj
        sp_avg += sp
np.save("./qc_y.npy", {"freq" : x, "sp" : sp_avg/10})
