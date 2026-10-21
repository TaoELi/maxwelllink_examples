"""RT-TDDFTB methylene blue in an NPoM cavity using cylindrical (r,z), m=0 FDTD."""

from pathlib import Path
import os
import shutil
import subprocess
import time
import random

import numpy as np
import matplotlib.pyplot as plt
import meep as mp
import maxwelllink as mxl
from maxwelllink.cavity import NPoM
from maxwelllink.tools import rt_tddft_spectrum
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
is_master = mp.am_master()

ROOT = Path.cwd()
INPUT_DIR = ROOT / "dftbplus_input"
RUN_ROOT = ROOT / "dftbplus_runs"
RUN_ROOT.mkdir(exist_ok=True)

if is_master:
    print("Using MPI processes:", comm.Get_size(), flush=True)
    print("ROOT:", ROOT, flush=True)
    print("INPUT_DIR:", INPUT_DIR, flush=True)
    print("RUN_ROOT:", RUN_ROOT, flush=True)

AU_TO_FS = 0.024188843265857
LIGHT_SPEED = 299792458  # m/s

os.environ["DFTBPLUS_BIN"] = "/home/4207/.local/bin/dftb+"
DFTBPLUS_BIN = os.environ.get("DFTBPLUS_BIN") or shutil.which("dftb+")

if DFTBPLUS_BIN is None:
    raise RuntimeError("Could not find DFTB+. Set DFTBPLUS_BIN to a socket-enabled dftb+ executable.")

if is_master:
    print("Using DFTB+:", DFTBPLUS_BIN, flush=True)


def normalize_spectrum(sp):
    sp = np.asarray(sp, dtype=float)
    scale = np.nanmax(np.abs(sp)) if sp.size else 0.0
    return sp / scale if scale > 0.0 else sp


def prepare_molecular_dynamics(socket_name: str, dt_au_value: float, steps: int):
    run_dir = RUN_ROOT
    run_dir.mkdir(parents=True, exist_ok=True)

    text = (INPUT_DIR / "methylene-blue_socket.hsd").read_text()
    text = text.replace("__MXL_SOCKET__", socket_name)
    text = text.replace("__DT_AU__", f"{dt_au_value:.16e}")
    text = text.replace("__STEPS__", str(int(steps)))
    (run_dir / "dftb_in.hsd").write_text(text)

    return run_dir


def launch_dftb_molecular_dynamics(run_dir: Path, sleep_time: float = 0.5):
    stdout = open(run_dir / "dftbplus.out", "w")

    print("Launching DFTB+ in", run_dir, flush=True)

    proc = subprocess.Popen(
        [DFTBPLUS_BIN],
        cwd=run_dir,
        stdout=stdout,
        stderr=subprocess.STDOUT,
    )

    time.sleep(sleep_time)

    if proc.poll() is not None:
        stdout.close()
        tail = (run_dir / "dftbplus.out").read_text(errors="replace")[-4000:]
        raise RuntimeError(f"DFTB+ exited before socket coupling started. Tail of dftbplus.out:\n{tail}")

    print("DFTB+ PID:", proc.pid, flush=True)
    return proc, stdout


# NPoM cavity: cylindrical (r,z), m = 0
resolution = 4000
courant = 0.50
cavity = NPoM(
    radius_nm=20.0,
    gap_nm=1.0,
    spacer_index=1.4,
    film_nm=70.0,
    padding_nm=300.0,
    pml_nm=300.0,
    resolution=resolution,
    dimensions=mp.CYLINDRICAL,
)

time_units_fs = cavity.time_units_fs
dt_meep = courant / resolution
dt_fs = dt_meep * time_units_fs
dt_au = dt_fs / AU_TO_FS

# set simulation time
total_time_meep = 30.0
meep_steps = int(round(total_time_meep / dt_meep))
dftb_steps = meep_steps + 10

if is_master:
    print(cavity.summary() if hasattr(cavity, "summary") else cavity, flush=True)
    print("time_units_fs:", time_units_fs, flush=True)
    print("FDTD dt [Meep]:", dt_meep, flush=True)
    print("FDTD dt [fs]:", dt_fs, flush=True)
    print("DFTB+ dt [a.u.]:", dt_au, flush=True)
    print("Meep steps:", meep_steps, flush=True)
    print("DFTB+ steps:", dftb_steps, flush=True)

# Socket
if is_master:
    socket = f"dftbplus_methylene-blue_{random.randint(100000, 999999)}"
else:
    socket = None

socket = comm.bcast(socket, root=0)

if is_master:
    print("Socket name:", socket, flush=True)

hub = mxl.SocketHub(unixsocket=socket, timeout=20.0, latency=1e-5)

# Methylene blue molecule
# In cylindrical m=0, the self-consistent molecule-field coupling is restricted to the z molecular dipole channel.
# The target MB transition dipole should therefore be aligned along z.

molecule = cavity.place_molecule(
    hub=hub,
    size_nm=1.0,
    sigma_nm=0.5,
    #polarization_type="transverse",
    polarization_type="analytical-finite",
    remove_static_self_field=1.4**2,
    store_additional_data=True,
    rescaling_factor=1.0,
)


# External optical source
# The cylindrical NPoM optical setup uses an Ez-polarized ring source, which belongs to m=0 and drives 
# the rotationally symmetric gap mode.
# The source center frequency is resonant with the DFTB+ MB excitation:
# 2.269 eV -> approximately 546 nm.

resonance_ev = 1.88233907  #unit: eV
f0 = resonance_ev / 1.23984667 # the energy of wavelength of 1000 nm (in unit eV)
#f0 = 1.830065  # 2.269 eV, approximately 546 nm
df = 0.5
cutoff = 3
start_time = -cutoff / df
source_amplitude = 2e3

optical_setup = cavity.optical_setup()

pulse = mp.GaussianSource(
    frequency=f0,
    fwidth=df,
    cutoff=cutoff,
    start_time=start_time,
)

sources_non_molecule = [
    mp.Source(
        src=pulse,
        component=optical_setup["component"],
        center=optical_setup["excitation"]["center"],
        size=optical_setup["excitation"]["size"],
        amplitude=source_amplitude,
    )
]
# check source setting
if mp.am_master():
    print("source component:", optical_setup["component"])
    print("source center:", optical_setup["excitation"]["center"])
    print("source size:", optical_setup["excitation"]["size"])

# MaxwellLink + Meep simulation
sim = cavity.make_simulation(
    molecules=[molecule],
    hub=hub,
    sources=sources_non_molecule,
    Courant=courant,
)


# Record the field at the NPoM hotspot
p = cavity.hotspot_center
if mp.am_master():
    print("hotspot center:", cavity.hotspot_center)
field_time_meep = []
Er_history = []
Ep_history = []
Ez_history = []

def step_function(sim):
    field_time_meep.append(sim.meep_time())
    Er_history.append(np.real(sim.get_field_point(mp.Er, p)))
    Ep_history.append(np.real(sim.get_field_point(mp.Ep, p)))
    Ez_history.append(np.real(sim.get_field_point(mp.Ez, p)))

# Launch one DFTB+ process and run all Meep MPI ranks
proc = None
stdout = None
run = None
try:
    if is_master:
        run = prepare_molecular_dynamics(socket, dt_au, dftb_steps)
        proc, stdout = launch_dftb_molecular_dynamics(run)

    comm.Barrier()
    sim.run(mp.at_every(dt_meep, step_function), steps=meep_steps)

finally:
    if is_master:
        if proc is not None and proc.poll() is None:
            proc.terminate()

            try:
                proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=10)

        if stdout is not None and not stdout.closed:
            stdout.close()
comm.Barrier()

# Postprocessing on the master rank only
if is_master:
    print("DFTB+ run directory:", run, flush=True)
    history = molecule.additional_data_history
    if len(history) == 0:
        raise RuntimeError("No molecular history was received from DFTB+.")
    print("First molecular record:", history[0], flush=True)
    print("Last molecular record:", history[-1], flush=True)
    time_au = np.asarray([np.real(ad["time_au"]) for ad in history])
    time_fs = time_au * AU_TO_FS
    time_ps = time_fs / 1e3

    mux_au = np.asarray([np.real(ad["mux_au"]) for ad in history])
    muy_au = np.asarray([np.real(ad["muy_au"]) for ad in history])
    muz_au = np.asarray([np.real(ad["muz_au"]) for ad in history])
    energy_au = np.asarray([np.real(ad["energy_au"]) for ad in history])

    field_time_meep = np.asarray(field_time_meep)
    field_time_fs = field_time_meep * time_units_fs
    Er_history = np.asarray(Er_history)
    Ep_history = np.asarray(Ep_history)
    Ez_history = np.asarray(Ez_history)

    np.savetxt(
        "dipole_moment",
        np.column_stack([time_fs, mux_au, muy_au, muz_au]),
        header="time_fs mux_au muy_au muz_au",
    )

    np.savetxt(
        "electric_field",
        np.column_stack([field_time_fs, Er_history, Ep_history, Ez_history]),
        header="time_fs Er Ep Ez",
    )

    np.savetxt(
        "energy_evolution",
        np.column_stack([time_fs, energy_au]),
        header="time_fs energy_au",
    )


    # Fourier transform of the z dipole
    # Only the z molecular channel is self-consistently coupled in cylindrical m=0, so muz is the relevant component here.

    if len(time_au) >= 2:
        nskip = 5
        mu = muz_au[::nskip]
        mu = mu - mu[0]
        spectrum_dt_au = (time_au[1] - time_au[0]) * nskip

        freq, intensity, _, _ = rt_tddft_spectrum(
            mu,
            dt_au=spectrum_dt_au,
            sp_form="absolute",
            e_start_ev=0.5,
            e_cutoff_ev=25.0,
            sigma=1.0e4,
            w_step=2.0e-4,
        )
        np.savetxt(
            "excitation_spectrum",
            np.column_stack([freq, normalize_spectrum(intensity)]),
            header="energy_eV normalized_spectral_power",
        )

    # find peak position in 1.8-2.8 eV
    mask = (freq >= 1.8) & (freq <= 2.8)
    freq_window = freq[mask]
    intensity_window = intensity[mask]
    imax = np.argmax(intensity_window)
    peak_freq = freq_window[imax]
    peak_intensity = intensity_window[imax]
    print("Peak frequency =", peak_freq)
    print("Peak intensity =", peak_intensity)
    with open("peak_position", "w") as f:
        f.write(f"{peak_freq:.8f}\n")

    # Plot dipole moment
    plt.figure(figsize=(10, 4))
    plt.plot(time_fs, mux_au, label=r"$\mu_x$")
    plt.plot(time_fs, muy_au, label=r"$\mu_y$")
    plt.plot(time_fs, muz_au, label=r"$\mu_z$")
    plt.xlabel("time [fs]")
    plt.ylabel("dipole moment (a.u.)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("dipole_moment.png", dpi=300)
    plt.close()
    # Plot cylindrical electric field
    plt.figure(figsize=(10, 4))
    plt.plot(field_time_fs, Er_history, label=r"$E_r$")
    plt.plot(field_time_fs, Ep_history, linestyle="--", label=r"$E_\phi$")
    plt.plot(field_time_fs, Ez_history, label=r"$E_z$")
    plt.xlabel("time [fs]")
    plt.ylabel("electric field [m.u.]")
    plt.legend()
    plt.tight_layout()
    plt.savefig("electric_field.png", dpi=300)
    plt.close()
    # Plot molecular energy
    plt.figure(figsize=(10, 4))
    plt.plot(time_fs, energy_au - energy_au[0], label="molecular energy")
    plt.xlabel("time [fs]")
    plt.ylabel("energy (a.u.)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("energy_evolution.png", dpi=300)
    plt.close()
    # Plot excitation spectrum
    if len(time_au) >= 2:
        plt.figure(figsize=(10, 4))
        plt.plot(freq, normalize_spectrum(intensity), label="excitation spectrum")
        plt.xlabel("energy [eV]")
        plt.ylabel("Normalized spectral power [arb. u.]")
        plt.legend()
        plt.tight_layout()
        plt.savefig("excitation_spectrum.png", dpi=300)
        plt.close()
    print("Simulation and postprocessing completed.", flush=True)
