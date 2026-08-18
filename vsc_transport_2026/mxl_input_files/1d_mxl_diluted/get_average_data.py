import numpy as np
import h5py
from maxwelllink.tools import ir_spectrum
import sys

nframe = 9
amp = [5000, 2000, 1000, 800, 500, 100]

for i in range(1):

    efield_y_avg = 0
    phe_avg = 0
    qc_sp = 0
    amplitude = amp[i]

    for idx in range(nframe+1):
        print(f'get data from multimode_cavmd_{idx}_{amplitude}_neq.h5')
        with h5py.File(f"./data/multimode_cavmd_{idx}_{amplitude}_neq.h5", "r") as f:
            efield_y_avg += f["effective_efield"][:, : ,1]**2    
            phe_avg += f["photonic_energy"][:, :]
            qc_y = f["qc"][:, :, 1]
            sp = np.zeros((qc_y.shape[0]//2, qc_y.shape[1]))
            for j in range(qc_y.shape[1]):
                x, spj = ir_spectrum(qc_y[:,j], 2)
                sp[:,j] = spj
            qc_sp += sp

    np.save(f"./scpdata/photonic_energy_36_{amplitude}.npy", phe_avg / (nframe+1))
    np.save(f"./scpdata/qc_spectra_36_{amplitude}.npy", {"freq": x, "sp": qc_sp / (nframe+1)})
    print(f"shape {phe_avg.shape} in photonic_energy_36_{amplitude}_neq.npy")
    print(f"shape {qc_sp.shape} in qc_spectra_36_{amplitude}_neq.npy")
    
    np.save(f"./scpdata/effective_efield_{amplitude}.npy", efield_y_avg / (nframe+1))    
    print(f"shape {efield_y_avg.shape} in effective_efield_{amplitude}_neq.npy")
