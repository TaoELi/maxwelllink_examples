import numpy as np
import maxwelllink as mxl
from maxwelllink import sockets as mxs
import meep as mp

host, port = mxs.get_available_host_port(localhost=False)
hub = mxl.SocketHub(host=host, port=port, timeout=200.0, latency=1e-4)

print(f"SocketHub listening on {host}:{port}")

# save host and port number to a file so mxl_driver can read it
with open("tcp_host_port_info.txt", "w") as f:
    f.write(f"{host}\n{port}\n")


n_molecule=2

molecules = []

for i in range(n_molecule):
    molecule = mxl.Molecule(
        hub=hub,
        center=mp.Vector3(0, 0, 0),
        size=mp.Vector3(1, 1, 1),
        sigma=0.1,
        dimensions=2,
        # ensure the superradiance emission rate is independent of n_molecule
        rescaling_factor=1.0 / n_molecule**0.5,
        # store the whole trajectory for only the first molecule to save memory
        store_additional_data=True if i == 0 else False,
    )
    molecules.append(molecule)

sim = mxl.MeepSimulation(
    hub=hub,
    molecules=molecules,
    cell_size=mp.Vector3(8, 8, 0),
    boundary_layers=[mp.PML(3.0)],
    resolution=10,
    time_units_fs=0.1,
)

sim.run(until=400)

# 4. Obtain necessary data for post-processing
if mp.am_master():
    print("Simulation completed. Collecting data...")
    Pe = np.array([ad["Pe"] for ad in molecules[0].additional_data_history])
    time_au = np.array([ad["time_au"] for ad in molecules[0].additional_data_history])
    np.savez(
        "tls_se_data.npz",
        time_au=time_au,
        Pe=Pe,
    )
