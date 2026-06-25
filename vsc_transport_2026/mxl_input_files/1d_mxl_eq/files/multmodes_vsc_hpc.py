import numpy as np
import maxwelllink as mxl
from maxwelllink.tools import gaussian_pulse, k_parallel_pulse
from maxwelllink.tools.harmonic_oscillator_helper import MaxwellBoltzmannInitializer, LangevinThermostat
from maxwelllink import sockets as mxs
from maxwelllink.units import unit

rng = np.random.default_rng(seed=114514)
random_seed = rng.integers(1919810)
N_grid = 144
host, port = mxs.get_available_host_port(localhost=False)
hub = mxl.SocketHub(host=host, port=port, timeout=10.0, latency=1e-4)
print(f"SocketHub listening on {host}:{port}")
# save host and port number to a file so mxl_driver can read it
with open("tcp_host_port_info.txt", "w") as f:
    f.write(f"{host}\n{port}\n")

molecule_list  = []
for i in range(N_grid):
    molecule = mxl.Molecule(
        hub=hub,
        rescaling_factor=1.0,
        store_additional_data=False,
    )
    molecule_list.append(molecule)

coupling_strength = 5e-5 / np.sqrt(N_grid//36)
print(f"Coupling strength: {coupling_strength:.3e} au")
damping_au = 0e-4
dt_fs = 0.5

fbcavity = mxl.FabryPerotCavity(
        frequency_au=2320*unit("cm_inv"),
        coupling_strength=coupling_strength,
        coupling_axis="xy",
        n_grid_x=N_grid,
        y_grid_1d=[0.0],
        delta_omega_x_au=12.5*unit("cm_inv"),
        delta_omega_y_au=0.0*unit("cm_inv"),
        n_mode_x=144,
        n_mode_y=1,
        abc_cutoff=0.0,
        save_mode_functions=True,
)

sim = mxl.MultiModeSimulation(
    hub=hub,
    molecules=molecule_list,
    dt_au=dt_fs*unit("fs"),
    damping_au=damping_au,
    include_dse=True,
    cavity_geometry=fbcavity,
    initializer=MaxwellBoltzmannInitializer(temperature_au=300.0*unit("K"), random_seed=random_seed),
    thermostat=LangevinThermostat(temperature_au=300.0*unit("K"), dt_au=dt_fs*unit("fs"), tau_au=5000*unit("fs"), random_seed=random_seed)
)

sim.run(steps=int(20000/dt_fs),
        record_history=True,
        record_to_disk=True,
        disk_folder_address='./',
        h5_filename=f"multimode_cavmd_4t4_eq.h5",
        record_every_steps=4,
        record_list=['qc'],)

