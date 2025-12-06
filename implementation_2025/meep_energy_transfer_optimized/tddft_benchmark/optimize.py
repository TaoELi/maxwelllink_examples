import psi4

psi4.set_memory("1 GB")
psi4.set_num_threads(1)

xyz = """
nocom
noreorient
0 1
    C            0.000000000000    -0.000000000000    -0.500327002867
    N            0.000000000000     0.000000000000     0.657588473262
    H            0.000000000000     0.000000000000    -1.577252899094
unit angstrom
"""

mol = psi4.geometry(xyz)

opts = {
            "basis": "cc-pvdz",
            "reference": "rks",
            "scf_type": "out_of_core",
            "e_convergence": 1e-12,
            "d_convergence": 1e-12,
            "guess": "sad",
            "puream": True,
            "g_convergence": "GAU_VERYTIGHT",
            "dft_grid_name": "SG1"
        }

psi4.set_options(opts)

energy, wfn = psi4.optimize("b3lyp/cc-pvdz", molecule=mol, return_wfn=True)
