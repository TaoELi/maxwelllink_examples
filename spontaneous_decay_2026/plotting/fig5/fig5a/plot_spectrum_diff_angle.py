"""Plot the scattering spectrum produced by run_spectrum.py."""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
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

cmap = LinearSegmentedColormap.from_list(
    "red_to_blue",
    ["indianred", "royalblue"]
)
colors = cmap(np.linspace(0, 1, 5))

colors = [
    "#d73027",  # red
    "#f46d43",  # orange-red
    "#ab5ab5",  # purple
    #"#5e67d8",  # blue-purple
    "#00A6D6",  # cyan-blue
    "#2b6cb0",  # blue
]

data_10 = np.loadtxt("npom_methylene_blue_scattering_spectrum_10.csv",delimiter=",",skiprows=1)
data_30 = np.loadtxt("npom_methylene_blue_scattering_spectrum_30.csv",delimiter=",",skiprows=1)
data_50 = np.loadtxt("npom_methylene_blue_scattering_spectrum_50.csv",delimiter=",",skiprows=1)
data_70 = np.loadtxt("npom_methylene_blue_scattering_spectrum_70.csv",delimiter=",",skiprows=1)
data_90 = np.loadtxt("npom_methylene_blue_scattering_spectrum_90.csv",delimiter=",",skiprows=1)
    
wavelength_10 = data_10[:, 0]
scattering_10 = data_10[:, 1]
wavelength_30 = data_30[:, 0]
scattering_30 = data_30[:, 1]
wavelength_50 = data_50[:, 0]
scattering_50 = data_50[:, 1]
wavelength_70 = data_70[:, 0]
scattering_70 = data_70[:, 1]
wavelength_90 = data_90[:, 0]
scattering_90 = data_90[:, 1]

plt.plot(wavelength_10, scattering_10, color=colors[0], label="10°")
plt.plot(wavelength_30, scattering_30, color=colors[1], label="30°")
plt.plot(wavelength_50, scattering_50, color=colors[2], label="50°")
plt.plot(wavelength_70, scattering_70, color=colors[3], label="70°")
plt.plot(wavelength_90, scattering_90, color=colors[4], label="90°")
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
    ncol=1,
    #loc="lower center",
    #bbox_to_anchor=(0.5, 1.00),
    #borderaxespad=0.0,
    columnspacing=0.2,
    handlelength=1.0,
    handletextpad=0.2
)

plt.tight_layout()
plt.savefig("scattering_spectrum_diff_angle.png", dpi=400)
plt.show()
