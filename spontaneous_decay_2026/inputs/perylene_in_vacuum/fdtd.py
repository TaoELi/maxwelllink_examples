from pathlib import Path
import os
import shutil
import subprocess
import time
import random

import matplotlib.pyplot as plt
import numpy as np

import meep as mp
import maxwelllink as mxl
from maxwelllink.tools import ir_spectrum, rt_tddft_spectrum

ROOT = Path.cwd()
INPUT_DIR = ROOT / "dftbplus_input"
RUN_ROOT = ROOT / "dftbplus_runs"
RUN_ROOT.mkdir(exist_ok=True)

print("ROOT:", ROOT)
print("INPUT_DIR:", INPUT_DIR)
print("RUN_ROOT:", RUN_ROOT)

AU_TO_CM = 219474.6313705
AU_TO_EV = 27.211386245988
AU_TO_FS = 0.024188843265857
H_TO_EV = 4.135667696e-15

os.environ["DFTBPLUS_BIN"] = "/home/4207/.local/bin/dftb+"
DFTBPLUS_BIN = os.environ.get("DFTBPLUS_BIN") or shutil.which("dftb+")
if DFTBPLUS_BIN is None:
    raise RuntimeError("Could not find DFTB+. Set DFTBPLUS_BIN to a socket-enabled dftb+ executable.")

print("Using DFTB+:", DFTBPLUS_BIN)

def normalize_spectrum(sp):
    sp = np.asarray(sp, dtype=float)
    scale = np.nanmax(np.abs(sp)) if sp.size else 0.0
    return sp / scale if scale > 0.0 else sp


def fft_spectrum_cm(signal, dt_au):
    signal = np.asarray(signal, dtype=float)
    centered = signal - np.mean(signal)
    window = np.hanning(centered.size)
    spectrum = np.abs(np.fft.rfft(centered * window)) ** 2
    omega_au = 2.0 * np.pi * np.fft.rfftfreq(centered.size, d=dt_au)
    return omega_au[1:] * AU_TO_CM, spectrum[1:]

# Prepare and launch DFTB+
def prepare_molecular_dynamics(socket_name: str, dt_au_value: float, steps: int):
    run_dir = RUN_ROOT
    run_dir.mkdir(parents=True, exist_ok=True)

    text = (INPUT_DIR / "perylene_socket.hsd").read_text()
    text = text.replace("__MXL_SOCKET__", socket_name)
    text = text.replace("__DT_AU__", f"{dt_au_value:.16e}")
    text = text.replace("__STEPS__", str(int(steps)))
    (run_dir / "dftb_in.hsd").write_text(text)
    return run_dir


def launch_dftb_molecular_dynamics(run_dir: Path, sleep_time: float = 0.5):
    stdout = open(run_dir / "dftbplus.out", "w")
    print("Launching DFTB+ in", run_dir)
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
    return proc, stdout

# Run the coupled electronic simulation
socket = f"dftbplus_perylene_{random.randint(100000, 999999)}"
resolution = 10
time_units_fs = 0.1
dt_au = time_units_fs / resolution / 2 * (1 / AU_TO_FS)
print('dt_au:', dt_au)
meep_steps = 200000
dftb_steps = 200010

hub = mxl.SocketHub(unixsocket=socket, timeout=20.0, latency=1e-5)
mol = mxl.Molecule(
                      hub=hub,
                      center=mp.Vector3(0, 0, 0),
                      size=mp.Vector3(1, 1, 1),
                      sigma=0.1,
                      resolution=resolution,
                      polarization_type="transverse",
                      dimensions=3,
                      rescaling_factor=5e1
                      )

# source setting
resonance_ev = 2.634  #unit: eV
frequency_meep = (resonance_ev / H_TO_EV) * (time_units_fs * 1e-15)
dE_eV = 0.5  #width of the excitation
fwidth = (dE_eV / H_TO_EV) * (time_units_fs * 1e-15)
cutoff = 3
#start_time = cutoff / fwidth * -1
start_time = 0e0
print('frequency, fwidth:', frequency_meep, fwidth)
sources_non_molecule = [
                        mp.Source(
                            src=mp.GaussianSource(
                                frequency=frequency_meep,
                                fwidth=fwidth,
                                cutoff=cutoff,
                                start_time=start_time),
                            component=mp.Ez,
                            center=mp.Vector3(0, 0, 0),   #source locates at molecular center
                            size=mp.Vector3(1.0, 1.0, 1.0),     #cover molecule box
                            amplitude=1e-2,)]

sim = mxl.MeepSimulation(
    hub=hub,
    molecules=[mol],
    sources=sources_non_molecule,
    cell_size=mp.Vector3(8, 8, 8),
    boundary_layers=[mp.PML(3.0)],
    resolution=resolution,
    time_units_fs=time_units_fs,  # fix a units system
)


# Run simulation
p = mp.Vector3(0,0,0)
Ez_history = []
def step_function(sim):
    Ez = sim.get_field_point(mp.Ez,p).real
    Ez_history.append(Ez)

run = prepare_molecular_dynamics(socket, dt_au, dftb_steps)
proc, stdout = launch_dftb_molecular_dynamics(run)
sim.run(
           mp.at_every(0.05, step_function),
           steps=meep_steps
           )
print("DFTB+ run directory:", run)
print(mol.additional_data_history[0])
print(mol.additional_data_history[-1])

# get and save data
history = mol.additional_data_history
t_au = np.array([ad["time_au"] for ad in history])
t_fs = t_au * AU_TO_FS
muz_au = np.array([ad["muz_au"] for ad in history])
energy_au = np.array([ad["energy_au"] for ad in history])
energy_kin_au = np.array([ad["energy_kin_au"] for ad in history])
energy_pot_au = np.array([ad["energy_pot_au"] for ad in history])

np.savetxt(
    "dipole_moment",
    np.column_stack([t_fs, muz_au]),
    #delimiter="",
    header="time_fs muz_au",
    #comments="",
)

np.savetxt(
    "electric_field",
    np.column_stack([t_fs[1:len(t_fs)], Ez_history]),
    #delimiter="",
    header="time_fs Ez",
    #comments="",
)

np.savetxt(
    "energy_evolution",
    np.column_stack([t_fs, energy_au-energy_au[0]]),
    #delimiter="",
    header="time_fs energy_au",
    #comments="",
)

#Fourier transform
nskip = 5
mu = muz_au[::nskip]
mu = mu - mu[0]          # remove static component
dt_au = (t_au[1] - t_au[0]) * nskip
freq, intensity, _, _ = rt_tddft_spectrum(
    mu,
    dt_au=dt_au,
    sp_form="absolute",
    e_start_ev=0.5,
    e_cutoff_ev=25.0,
    sigma=1.0e4,
    w_step=2.0e-4,
)

# Plot the dipole moment, electric field, energy and the excitation spectrum
plt.figure(figsize=(10, 4))
plt.plot(t_fs, muz_au, label=r"$\mu_z$")
plt.xlabel("Time (fs)")
plt.ylabel("Dipole moment (a.u.)")
#plt.title("H$_2$ electronic coupling energy trace")
plt.legend()
plt.tight_layout()
plt.savefig("dipole_moment.png", dpi=300)

plt.figure(figsize=(10, 4))
plt.plot(t_fs[1:len(t_fs)], Ez_history, label=r"$E_z$")
plt.xlabel("Time (fs)")
plt.ylabel("Electric field (m.u.)")
plt.legend()
plt.tight_layout()
plt.savefig("electric_field.png", dpi=300)

plt.figure(figsize=(10, 4))
#plt.plot(t_fs, energy_au, label="molecular energy")
plt.plot(t_fs, energy_au-energy_au[0], label="molecular energy")
#plt.plot(t_fs, energy_kin_au-energy_kin_au[0], label="kinetic energy")
#plt.plot(t_fs, energy_pot_au-energy_pot_au[0], label="potential energy")
plt.xlabel("Time (fs)")
plt.ylabel("Energy (a.u.)")
plt.legend()
plt.tight_layout()
plt.savefig("energy_evolution.png", dpi=300)

plt.figure(figsize=(10, 4))
plt.plot(freq, normalize_spectrum(intensity), label="excitation spectrum")
plt.xlabel("Energy (eV)")
plt.ylabel("Normalized spectral power")
plt.legend()
plt.tight_layout()
plt.savefig("excitation_spectrum.png", dpi=300)
plt.show()
