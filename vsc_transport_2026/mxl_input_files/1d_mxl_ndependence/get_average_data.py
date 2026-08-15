import numpy as np
import h5py
from maxwelllink.tools import ir_spectrum
import sys

nframe = 9
gird_list = [144, 288, 576, 1152] 

for grid in gird_list:
    amplitude = 0.007
    up = 36
    efield_y_avg = 0

    for idx in range(nframe+1):
        print(f'get data from multimode_cavmd_{grid}_{idx}_{up}_{amplitude}_neq.h5')
        with h5py.File(f"./data/multimode_cavmd_{grid}_{idx}_{up}_{amplitude}_neq.h5", "r") as f:
            efield_y_avg += np.mean(np.reshape(f["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2
    
    np.save(f"./scpdata/effective_efield_{grid}_{up}_{amplitude}.npy", efield_y_avg / (nframe+1))    
    print(f"shape {efield_y_avg.shape} in effective_efield_{grid}_{up}_{amplitude}_neq.npy")


