import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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
# read data
data = np.genfromtxt('energy_evolution', skip_header=1)
time_fs = data[:, 0]
energy_au = data[:, 1]
energy_au = energy_au * 1e4
time_ps = np.array(time_fs) / 1e3

t0_step = 48000
# define function
def decay_func(t, E0, tau, E_infinite):
    #decay function
    return E0 * np.exp(-(t-time_ps[t0_step]) / tau) + E_infinite

# fitting parameters
popt, pcov = curve_fit(decay_func, time_ps[t0_step:len(time_ps)], energy_au[t0_step:len(energy_au)]-energy_au[0], p0=[energy_au[t0_step]-energy_au[0], 2e-2, 0e0])

# save data
E0 = popt[0]
tau = popt[1]
E_infinite = popt[2]
print('tau= ', tau, 'ps')
print(f"decay rate: 1/tau = {1/tau:.4f} ps^-1")
print('E0= ', E0, 'au')
print('E_infinite= ', E_infinite, 'au')

np.savetxt(
    "decay_rate",
    np.column_stack([tau, 1/tau, E0, E_infinite]),
        #delimiter="",
        header="lifetime[ps] decay_rate[ps^-1] E0[a.u.] E_infinite[a.u.]",
        #comments="",
    )
# plot

#plt.figure(figsize=(6,4))
#plt.plot(time_ps[t0_step:len(time_ps)], energy_au[t0_step:len(energy_au)]-energy_au[0], 'b-', label='Data')
plt.plot(time_ps, energy_au-energy_au[0], color='indianred', label='simulation')
plt.plot(time_ps[t0_step:len(time_ps)], decay_func(time_ps[t0_step:len(time_ps)], E0, tau, E_infinite), "--",color='black', label='fitting')
plt.xlabel('time [ps]')
#plt.ylabel(r"$E - E_0$ [a.u.]")
plt.ylabel(r"$E - E_0$ [10$^{-4}$ a.u.]")
plt.ylim(0.0, 1.5)
plt.yticks(np.arange(0.0, 1.50001, 0.5))
plt.text(0.05, 0.3, r'$\tau$ =11.1 fs', fontsize=10)
#plt.title(f"decay rate: {1/tau:.4f} ps^-1")
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
plt.savefig('fit_decay.png',dpi=400)
plt.show()
