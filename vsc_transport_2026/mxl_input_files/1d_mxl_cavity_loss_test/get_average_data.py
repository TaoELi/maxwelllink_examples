import numpy as np
import h5py
from maxwelllink.tools import ir_spectrum
import sys

nframe = 9
amp = [5000, 2000, 1000, 800, 500, 100]

for i in range(6):

    efield_y_avg = 0
    amplitude = amp[i]

    for idx in range(nframe+1):
        print(f'get data from multimode_cavmd_{idx}_{amplitude}_neq.h5')
        with h5py.File(f"./data/multimode_cavmd_{idx}_{amplitude}_neq.h5", "r") as f:
            efield_y_avg += f["effective_efield"][:, : ,1]**2    
    
    np.save(f"./scpdata/effective_efield_{amplitude}.npy", efield_y_avg / (nframe+1))    
    print(f"shape {efield_y_avg.shape} in effective_efield_{amplitude}_neq.npy")

