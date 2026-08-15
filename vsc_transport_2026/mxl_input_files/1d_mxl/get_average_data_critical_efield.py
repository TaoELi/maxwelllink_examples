import numpy as np
import h5py
from maxwelllink.tools import ir_spectrum
import sys

nframe = 9
amp = ['0.0016', '0.0024']
up_list = ['24', '36'] 

for i in range(2):

    efield_y_avg = 0
    up = up_list[i]
    amplitude = amp[i]

    for idx in range(nframe+1):
        print(f'get data from multimode_cavmd_{idx}_{up}_{amplitude}_neq.h5')
        with h5py.File(f"./data/multimode_cavmd_{idx}_{up}_{amplitude}_neq.h5", "r") as f:
            efield_y_avg += f["effective_efield"][:, : ,1]**2    
    
    np.save(f"./scpdata/effective_efield_{up}_{amplitude}.npy", efield_y_avg / (nframe+1))    
    print(f"shape {efield_y_avg.shape} in effective_efield_{up}_{amplitude}_neq.npy")

