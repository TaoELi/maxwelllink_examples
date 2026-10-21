"""Calculate the NPoM + methylene-blue scattering spectrum using cylindrical FDTD + self-consistent RT-TDDFTB dynamics."""

from pathlib import Path
import os
import shutil
import subprocess
import time
import random

import numpy as np
import meep as mp
import maxwelllink as mxl
from maxwelllink.cavity import NPoM
from maxwelllink.measurements import MeepScatteringSpectroscopy
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
is_master = mp.am_master()

# DFTB+ paths
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

os.environ["DFTBPLUS_BIN"] = "/home/4207/.local/bin/dftb+"
DFTBPLUS_BIN = os.environ.get("DFTBPLUS_BIN") or shutil.which("dftb+")

if DFTBPLUS_BIN is None:
    raise RuntimeError("Could not find DFTB+. Set DFTBPLUS_BIN to a socket-enabled dftb+ executable.")

if is_master:
    print("Using DFTB+:", DFTBPLUS_BIN, flush=True)

# Prepare and launch DFTB+
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

# Cylindrical NPoM cavity
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

if is_master:
    cavity.summary()

# FDTD / RT-TDDFTB timestep
time_units_fs = cavity.time_units_fs
dt_meep = courant / resolution
dt_fs = dt_meep * time_units_fs
dt_au = dt_fs / AU_TO_FS

if is_master:
    print("time_units_fs:", time_units_fs, flush=True)
    print("FDTD dt [Meep]:", dt_meep, flush=True)
    print("FDTD dt [fs]:", dt_fs, flush=True)
    print("DFTB+ dt [a.u.]:", dt_au, flush=True)
# Scattering-spectrum measurement
# First construct the measurement without molecule or hub.
# The reference calculation contains only the Au film + spacer.

#source_amplitude = 1.0
measurement = MeepScatteringSpectroscopy(
    cavity,
    600.0,
    800.0,
    units="nm",
    nfreq=400,
    molecules=[],
    hub=None,
    max_time=50.0,
    #source_amplitude=source_amplitude,
    #Courant=courant,
)

if is_master:
    print("\nStarting reference simulation...", flush=True)

reference = measurement.reference()

if is_master:
    print("Reference simulation finished.", flush=True)


# ----------------------------------------------------------------------
# Create socket AFTER the reference simulation
# ----------------------------------------------------------------------

if is_master:
    socket = f"dftbplus_methylene-blue_{random.randint(100000, 999999)}"
else:
    socket = None

socket = comm.bcast(socket, root=0)

if is_master:
    print("Socket name:", socket, flush=True)

hub = mxl.SocketHub(unixsocket=socket, timeout=200.0, latency=1e-5)

# Methylene-blue molecule
# Cylindrical m=0 couples to the z-oriented molecular dipole.
# The MB geometry should therefore be oriented so that the transition dipole is along the z direction.

molecule = cavity.place_molecule(
    hub=hub,
    size_nm=1.0,
    sigma_nm=0.5,
    polarization_type="analytical-finite",
    remove_static_self_field=1.4**2,
    store_additional_data=False,
    rescaling_factor=1.0,
)

# Add the RT-TDDFTB molecule and socket hub to the same measurement object.
measurement.molecules = [molecule]
measurement.hub = hub

# DFTB+ step limit

dftb_time_cap_meep = 50.0
dftb_steps = int(np.ceil(dftb_time_cap_meep / dt_meep)) + 10

if is_master:
    print("Maximum DFTB+ steps:", dftb_steps, flush=True)

# Launch one DFTB+ process
proc = None
stdout = None
run = None

try:
    if is_master:
        run = prepare_molecular_dynamics(socket, dt_au, dftb_steps)
        proc, stdout = launch_dftb_molecular_dynamics(run)

    comm.Barrier()

    if is_master:
        print("\nStarting NPoM + methylene-blue signal simulation...", flush=True)

    signals = measurement.signal_run()

    if is_master:
        print("Signal simulation finished.", flush=True)

    spectrum = measurement.postprocess(reference, signals)

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

# Save scattering spectrum
if is_master:
    wavelength = spectrum["wavelength_nm"]
    scattering = spectrum["scattering"]
    absorption = spectrum["absorption"]
    extinction = spectrum["extinction"]

    order = np.argsort(wavelength)

    np.savetxt(
        "npom_methylene_blue_scattering_spectrum.csv",
        np.column_stack(
            [
                wavelength[order],
                scattering[order],
                absorption[order],
                extinction[order],
            ]
        ),
        delimiter=",",
        header="wavelength_nm,scattering,absorption,extinction",
        comments="",
    )

    print("\nSpectrum saved to npom_methylene_blue_scattering_spectrum.csv", flush=True)
    print("DFTB+ run directory:", run, flush=True)
