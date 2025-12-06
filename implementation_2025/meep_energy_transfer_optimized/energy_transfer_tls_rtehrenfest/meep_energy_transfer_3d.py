import numpy as np
import maxwelllink as mxl
from maxwelllink import sockets as mxs
import meep as mp

host, port = mxs.get_available_host_port(localhost=False)
hub = mxl.SocketHub(host=host, port=port, timeout=200.0, latency=1e-4)

if mp.am_master():
    print(f"SocketHub listening on {host}:{port}")
    # save host and port number to a file so mxl_driver can read it
    with open("tcp_host_port_info.txt", "w") as f:
        f.write(f"{host}\n{port}\n")

# define the separation between the donor and accceptor molecules
separation_over_2 = 0.3

molecule_donor = mxl.Molecule(
    hub=hub,
    center=mp.Vector3(-separation_over_2, 0, 0),
    size=mp.Vector3(1, 1, 1),
    sigma=0.1,
    dimensions=3,
    rescaling_factor=1.0,
)

molecule_acceptor = mxl.Molecule(
    hub=hub,
    center=mp.Vector3(separation_over_2, 0, 0),
    size=mp.Vector3(1, 1, 1),
    sigma=0.1,
    dimensions=3,
    rescaling_factor=1.0,
)

sim = mxl.MeepSimulation(
    hub=hub,
    molecules=[molecule_donor, molecule_acceptor],
    cell_size=mp.Vector3(8, 8, 8),
    boundary_layers=[mp.PML(3.0)],
    resolution=10,
    time_units_fs=0.02,
    # converted to ~0.04 au per driver step
)

sim.run(steps=4000)

# 4. Obtain necessary data for post-processing
if mp.am_master():
    print("Simulation completed. Collecting data...")
    # donor should be a simple TLS
    time_au = np.array([ad["time_au"] for ad in molecule_donor.additional_data_history])
    Pe = np.array([ad["Pe"] for ad in molecule_donor.additional_data_history])
    energy_au = np.array([ad["energy_au"] for ad in molecule_donor.additional_data_history])
    mux_au = np.array([ad["mux_au"] for ad in molecule_donor.additional_data_history])
    muy_au = np.array([ad["muy_au"] for ad in molecule_donor.additional_data_history])
    muz_au = np.array([ad["muz_au"] for ad in molecule_donor.additional_data_history])
    np.savez(
        "donor_tls_data.npz",
        time_au=time_au,
        Pe=Pe,
        energy_au=energy_au,
        mux_au=mux_au,
        muy_au=muy_au,
        muz_au=muz_au,
    )

    # acceptor can be a simple TLS, QuTiP molecule, or TDDFT molecule
    time_au = np.array([ad["time_au"] for ad in molecule_acceptor.additional_data_history])
    energy_au = np.array([ad["energy_au"] for ad in molecule_acceptor.additional_data_history])
    mux_au = np.array([ad["mux_au"] for ad in molecule_acceptor.additional_data_history])
    muy_au = np.array([ad["muy_au"] for ad in molecule_acceptor.additional_data_history])
    muz_au = np.array([ad["muz_au"] for ad in molecule_acceptor.additional_data_history])
    # for TLS
    result_tls = molecule_acceptor.additional_data_history[0].get("Pe", None)
    if result_tls is not None:
        Pe = np.array([ad["Pe"] for ad in molecule_acceptor.additional_data_history])
        np.savez(
            "acceptor_tls_data.npz",
            time_au=time_au,
            Pe=Pe,
            energy_au=energy_au,
            mux_au=mux_au,
            muy_au=muy_au,
            muz_au=muz_au,
        )
    # for QuTiP
    result_qutip = molecule_acceptor.additional_data_history[0].get("rho_diag", None)
    if result_qutip is not None:
        rho_diag = np.array([ad["rho_diag"] for ad in molecule_acceptor.additional_data_history])
        np.savez(
            "acceptor_qutip_data.npz",
            time_au=time_au,
            rho_diag=rho_diag,
            energy_au=energy_au,
            mux_au=mux_au,
            muy_au=muy_au,
            muz_au=muz_au,
        )
    # for TDDFT and Ehrenfest results
    if result_tls is None and result_qutip is None:
        # only for Ehrenfest 
        energy_kin_au = np.array([ad["energy_kin_au"] for ad in molecule_acceptor.additional_data_history])
        np.savez(
            "acceptor_tddft_data.npz",
            time_au=time_au,
            energy_au=energy_au,
            energy_kin_au=energy_kin_au,
            mux_au=mux_au,
            muy_au=muy_au,
            muz_au=muz_au,
        )
    
