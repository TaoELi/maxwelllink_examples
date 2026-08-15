import numpy as np
import h5py
from maxwelllink.tools import ir_spectrum
import sys

nframe = 9
amp = ['0.001', '0.003', '0.005', '0.007']
up_list = ['24', '27', '30', '36', '42', '48', '54', '60', '66'] 

for amplitude in amp:
    for up in up_list:

        efield_y_avg = 0
        coenergy_bond_avg = 0

        for idx in range(nframe+1):
            print(f'get data from multimode_cavmd_{idx}_{up}_{amplitude}_neq.h5')
            with h5py.File(f"./data/multimode_cavmd_{idx}_{up}_{amplitude}_neq.h5", "r") as f:
                efield_y_avg += f["effective_efield"][:, : ,1]**2    
            coenergy_bond_avg += np.load(f"./data/coenergy_{idx}_{up}_{amplitude}.npy")[:, :]
        
        np.save(f"./scpdata/effective_efield_{up}_{amplitude}.npy", efield_y_avg / (nframe+1))    
        print(f"shape {efield_y_avg.shape} in effective_efield_{up}_{amplitude}_neq.npy")

        np.save(f"./scpdata/coenergy_{up}_{amplitude}.npy", coenergy_bond_avg / (nframe+1))
        print(f"shape {coenergy_bond_avg.shape} in coenergy_{up}_{amplitude}_neq.npy")

