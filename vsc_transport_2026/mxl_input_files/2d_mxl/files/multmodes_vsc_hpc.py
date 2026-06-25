import numpy as np
import maxwelllink as mxl
from maxwelllink.tools import gaussian_pulse, k_parallel_pulse
from maxwelllink.tools.harmonic_oscillator_helper import MaxwellBoltzmannInitializer, LangevinThermostat
from maxwelllink import sockets as mxs
from maxwelllink.units import unit

rng = np.random.default_rng(seed=114514)
random_seed = rng.integers(1919810)
N_grid_x = 144
N_grid_y = 144
N_grid = N_grid_x * N_grid_y
host, port = mxs.get_available_host_port(localhost=False)
#hub = mxl.SocketHub(host=host, port=port, timeout=1e6, latency=1e-4)
# change #1: use AggregatedSocketHub instead of SocketHub
hub = mxl.AggregatedSocketHub(host=host, port=port, timeout=6000.0, latency=1e-3)
print(f"SocketHub listening on {host}:{port}")
# save host and port number to a file for diagnostics / legacy direct-driver runs
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

# change #2: initialize the remote bridges and allocate 72 molecules per bridge/node
# the manifest "aggregation.json" will be written to the shared filesystem for the bridge nodes to read
hub.init_remote_bridges(
    molecule_list,
    molecules_per_bridge=576,
    unix_prefix="bridge_",
    save_file="aggregation.json",
)

coupling_strength = 5e-5 / np.sqrt(N_grid//36)
print(f"Coupling strength: {coupling_strength:.3e} au")
damping_au = 0e-4
dt_fs = 0.5

fbcavity = mxl.FabryPerotCavity(
        frequency_au=2320*unit("cm_inv"),
        coupling_strength=coupling_strength,
        coupling_axis="xy",
        n_grid_x=N_grid_x,
        n_grid_y=N_grid_y,
        delta_omega_x_au=12.5*unit("cm_inv"),
        delta_omega_y_au=12.5*unit("cm_inv"),
        n_mode_x=N_grid_x,
        n_mode_y=N_grid_y,
        abc_cutoff=[0.05,0.05],
        save_mode_functions=True,
)

pulse_envelope = gaussian_pulse(
    amplitude_au=1.0,
    t0_au=0*unit("fs"),
    sigma_au=5e2*unit("fs"),
    t_start_au=0.0*unit("fs"),
    t_end_au=1e3*unit("fs"),
)
# the pulse is basically envelope * cos(k_parallel * x - omega * t)
# with the envelope being a simple broad-band Gaussian pulse
molecule_source = k_parallel_pulse(
    cavity=fbcavity,
    envelope=pulse_envelope,
    omega_au=2438.84*unit("cm_inv"),
    k_parallel_au=[36*12.5*unit("cm_inv"),1*12.5*unit("cm_inv")],
    direction="xy",
    center=(0.2, 0.5),
    size=(0.3, 0.3),
    target="photon",
    amplitude_au=0.002,
)


sim = mxl.MultiModeSimulation(
    hub=hub,
    molecules=molecule_list,
    dt_au=dt_fs*unit("fs"),
    damping_au=damping_au,
    include_dse=True,
    cavity_geometry=fbcavity,
    excited_mode_list=molecule_source.excited_mode_list,
    photon_pulse_drive=molecule_source,
    photon_pulse_axis= "y",
    photon_partial_charge=0.43,
    initializer=MaxwellBoltzmannInitializer(temperature_au=300.0*unit("K"), random_seed=random_seed),
    thermostat=LangevinThermostat(temperature_au=300.0*unit("K"), dt_au=dt_fs*unit("fs"), tau_au=5000*unit("fs"), random_seed=random_seed)
)

sim.run(steps=int(8001/dt_fs),
        record_history=True,
        record_to_disk=True,
        disk_folder_address='./',
        h5_filename=f"2d_{N_grid_x}_{N_grid_y}_efield_neq.h5",
        record_every_steps=40,
        record_list=['effective_efield'],)
