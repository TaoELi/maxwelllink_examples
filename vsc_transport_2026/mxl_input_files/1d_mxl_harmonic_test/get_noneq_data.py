import numpy as np
import h5py
from maxwelllink.tools import ir_spectrum

amp = [0.001, 0.007]
nframe = 9

for amplitude in amp:
    phe_avg = 0
    efield_avg = 0
    sp_avg = 0
    coenergy_avg = 0
    for i in range(nframe):
        print(f'get data from multimode_cavmd_{nframe}_96_{amplitude}_neq.h5')
        with h5py.File(f"./data/multimode_cavmd_{nframe}_96_{amplitude}_neq.h5", "r") as f:
            efield_avg += f["effective_efield"][:, : ,1]        
            phe_avg += f["photonic_energy"][:, :]
            qc_y = f["qc"][:, :, 1]
            sp = np.zeros((qc_y.shape[0]//2, qc_y.shape[1]))
            for j in range(qc_y.shape[1]):
                x, spj = ir_spectrum(qc_y[:,j], 2)
                sp[:,j] = spj
            sp_avg += sp
        coenergy_avg += np.load(f"./data/coenergy_{nframe}_96_{amplitude}.npy")
    
    np.save(f"./scpdata/photonic_energy_96_{amplitude}.npy", phe_avg / (nframe+1))
    print(f"shape {phe_avg.shape} in photonic_energy_96_{amplitude}_neq.npy")
    np.save(f"./scpdata/qc_spectra_96_{amplitude}.npy", {"freq": x, "sp": sp_avg / (nframe+1)})
    print(f"shape {sp_avg.shape} in qc_spectra_96_{amplitude}_neq.npy")
    np.save(f"./scpdata/effective_efield_96_{amplitude}.npy", efield_avg / (nframe+1))
    print(f"shape {efield_avg.shape} in effective_efield_96_{amplitude}_neq.npy")
    np.save(f"./scpdata/coenergy_{amplitude}_neq.npy", coenergy_avg / (nframe+1))
    print(f"shape {coenergy_avg.shape} in coenergy_{amplitude}_neq.npy")




