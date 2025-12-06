import numpy as np
import maxwelllink as mxl

dt_rttddft_au = 0.04

model = mxl.RTTDDFTModel(
    molecule_xyz="hcn.xyz",
    functional="b3lyp",
    basis="cc-pvdz",
    dt_rttddft_au=dt_rttddft_au,
    delta_kick_au=1.0e-3,
    delta_kick_direction="xyz",
    dft_grid_name="SG1"
)

model.initialize(dt_new=dt_rttddft_au, molecule_id=0)

# calculate LR-TDDFT spectrum
states = 60
poles, oscillator_strengths, res = model._get_lr_tddft_spectrum(states=states, tda=False, savefile=False)

# save LR-TDDFT results for plotting
np.savetxt(f"lrtddft_poles_au_{states}.txt", np.array(poles))
np.savetxt(f"lrtddft_oscillator_strengths_{states}.txt", np.array(oscillator_strengths))

# save LR-TDDFT results for multiscale benchmark
tdm_len = np.array(
                [r["ELECTRIC DIPOLE TRANSITION MOMENT (LEN)"] for r in res]
            )

print("Excitation Energies (a.u.):\n", poles)
print("Transition Dipole Moments (a.u.):\n", tdm_len)

# Identify the strongest absorption peak
idx = np.argmax(oscillator_strengths)
print(f"The strongest absorption peak is at {poles[idx]:.4f} a.u. ({poles[idx]*27.2114:.4f} eV) with oscillator strength {oscillator_strengths[idx]:.4f}.")
print(f"The transition dipole vector of this transition is {tdm_len[idx]} a.u.")

np.savetxt(f"lrtddft_tdm_len_{states}.txt", tdm_len)
