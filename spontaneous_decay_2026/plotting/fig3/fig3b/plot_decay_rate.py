import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
    'figure.figsize': [3.15, 2.36],
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


# =========================
# input files
# =========================
"""
datasets = [
    {
        "filename": "decay_rate-ethylene",
        "label": "ethylene",
        "color": "#0072B2",     # blue
    },
    {
        "filename": "decay_rate-anthracene",
        "label": "anthracene",
        "color": "#D55E00",     # vermillion/orange-red
    },
    {
        "filename": "decay_rate-perylene",
        "label": "perylene",
        "color": "#009E73",     # bluish green
    },
    {
        "filename": "decay_rate-DPA",
        "label": "DPA",
        "color": "#CC79A7",     # reddish purple
    },
]
"""
datasets = [
    {
        "filename": "decay_rate-ethylene",
        "label": "ethylene",
        "color": "#4C72B0",     # muted blue
    },
    {
        "filename": "decay_rate-anthracene",
        "label": "anthracene",
        "color": "#DD8452",     # muted orange
    },
    {
        "filename": "decay_rate-perylene",
        "label": "perylene",
        "color": "#55A868",     # muted green
    },
    {
        "filename": "decay_rate-DPA",
        "label": "DPA",
        "color": "#C44E52",     # muted red
    },
]

# =========================
# plot
# =========================
fig, ax = plt.subplots()

x_fit = np.linspace(0.0, 2500.0, 200)

for item in datasets:
    data = np.loadtxt(item["filename"])

    x = data[:, 0]   # gamma^2
    y = data[:, 1]   # decay rate

    # scatter points
    ax.scatter(
        x, y,
        color=item["color"],
        #edgecolor="black",
        #linewidth=0.3,
        s=25,
        label=item["label"],
        zorder=3
    )

    # fitting function: y = kx
    k = np.sum(x * y) / np.sum(x * x)
    y_fit = k * x_fit

    # fitting line
    ax.plot(
        x_fit, y_fit,
        color=item["color"],
        linestyle="-",
        linewidth=0.9,
        zorder=2
    )

    print(f"{item['label']}: decay_rate = {k:.6e} * gamma^2")


# =========================
# axis labels and range
# =========================
ax.set_xlabel(r"$\gamma^2$")
ax.set_ylabel(r"decay rate [ps$^{-1}$]")

ax.set_xlim(1, 3000)
ax.set_ylim(1e-5, 1e1)


# =========================
# legend
# =========================
handles, labels = ax.get_legend_handles_labels()

"""
fit_handle = Line2D(
    [0], [0],
    color="black",
    linestyle="-",
    linewidth=0.9,
    label="linear fitting"
)
"""

#handles.append(fit_handle)
#labels.append("linear fitting")

ax.legend(handles, labels, frameon=False)


# =========================
# save and show
# =========================
x_label = [1e0, 1.00001e0]
y_label = [7.01498e-5, 7.0149899001e-5]
plt.plot(x_label, y_label, "--",color='black', label='linear fitting')
plt.legend(frameon=False, ncol=2, columnspacing=0.4, loc=2, bbox_to_anchor=(-0.02, 1.02), borderaxespad=0)
plt.xscale("log")
plt.yscale("log")
plt.tight_layout()
plt.savefig("decay_rate_linear_fitting.png")
plt.show()
