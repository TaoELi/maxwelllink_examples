"""Plot the scattering spectrum produced by run_spectrum.py."""

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

# =========================
# plot style
# =========================
plt.rcParams.update({
    # font
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 8,

    # axis
    'axes.linewidth': 0.8,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'axes.labelpad': 2,
    'axes.titlepad': 4,
    'axes.labelcolor': 'black',
    'axes.edgecolor': 'black',

    # axis tick marks
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,

    # line and marker
    'lines.linewidth': 0.8,
    'lines.markersize': 3,
    'lines.markeredgewidth': 0.3,
    'patch.linewidth': 0.6,

    # legend
    'legend.frameon': False,
    'legend.fontsize': 10,
    'legend.handlelength': 1.2,
    'legend.handletextpad': 0.3,
    'legend.borderaxespad': 0.2,
    'legend.labelspacing': 0.2,

    # grid and background
    'axes.grid': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',

    # figure size and save
    'figure.figsize': [2.36, 2.36],
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'savefig.transparent': True,

    # math font
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
})

data_npom = np.loadtxt("npom_scattering_spectrum.csv",delimiter=",",skiprows=1)
data_rescaling_1 = np.loadtxt("npom_methylene_blue_scattering_spectrum.csv",delimiter=",",skiprows=1)
    
wavelength_npom = data_npom[:, 0]
scattering_npom = data_npom[:, 1]
wavelength_rescaling_1 = data_rescaling_1[:, 0]
scattering_rescaling_1 = data_rescaling_1[:, 1]

plt.plot(wavelength_npom, scattering_npom, color='royalblue', label="NPoM")
plt.plot(wavelength_rescaling_1, scattering_rescaling_1, color='indianred', label="NPoM + molecule")
plt.xlabel("wavelength [nm]")
plt.ylabel("scattering intensity")
#plt.title("NPoM + mb molecule scattering spectra")
# =========================
# hide y-axis ticks, labels, and scientific offset, e.g. 1e-5
# =========================
ax = plt.gca()
ax.set_yticks([])
ax.tick_params(axis="y", left=False, labelleft=False)
ax.yaxis.get_offset_text().set_visible(False)
# =========================
# legend outside, just above the top frame, in one row
# =========================
plt.legend(
    frameon=False,
    ncol=2,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.00),
    borderaxespad=0.0,
    columnspacing=0.2,
    handlelength=1.0,
    handletextpad=0.2
)

plt.tight_layout()
plt.savefig("scattering_spectrum.png", dpi=400)
plt.show()
