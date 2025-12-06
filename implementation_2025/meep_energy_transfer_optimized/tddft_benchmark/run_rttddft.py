import maxwelllink as mxl

dt_rttddft_au = 0.04

model = mxl.RTTDDFTModel(
    molecule_xyz="hcn.xyz",
    functional="b3lyp",
    basis="cc-pvdz",
    dt_rttddft_au=dt_rttddft_au,
    delta_kick_au=1.0e-2,
    delta_kick_direction="xyz",
    memory = "16GB",
    num_threads = 16,
    verbose = False,
    dft_grid_name="SG1"
)

model.initialize(dt_new=dt_rttddft_au, molecule_id=0)

# propagate standalone RT-TDDFT
times, energies, dipoles = model._propagate_full_rt_tddft(nsteps=40000, savefile=True)
