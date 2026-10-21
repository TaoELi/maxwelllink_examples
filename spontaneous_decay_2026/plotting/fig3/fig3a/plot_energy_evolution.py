import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({
    # font
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],       # use Arial
    'font.size': 8,                     # fontsize: 8pt

    # axis
    'axes.linewidth': 0.8,              # 坐标轴线条宽度
    'axes.labelsize': 10,                # 轴标签字号
    'axes.titlesize': 10,                # 标题字号
    'axes.labelpad': 2,                 # 标签与轴之间的距离
    'axes.titlepad': 4,                 # 标题与图之间的距离
    'axes.labelcolor': 'black',
    'axes.edgecolor': 'black',

    # Axis tick narks
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 2.5,            # 主刻度线长度
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

    # figure size and save the figure
    'figure.figsize': [3.15, 2.36],     # 单列图大小（约 8 cm × 6 cm）
    'savefig.dpi': 600,                 # 高分辨率导出
    'savefig.bbox': 'tight',            # 紧凑布局
    'savefig.pad_inches': 0.02,
    'savefig.transparent': True,

    # math font
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
})

# read data
def read_data(filename):
    #data = np.genfromtxt('energy_evolution-rescaling-1e1', skip_header=1)
    data = np.genfromtxt(filename, skip_header=1)
    time_fs = data[:, 0]
    energy_au = data[:, 1]
    time_ps = np.array(time_fs) / 1e3
    return energy_au, time_ps

# define function
def decay_func(t, E0, tau):
    #decay function
    return E0 * np.exp(-(t-0.060) / tau)


# fitting and plot
"""
plt.figure(figsize=(3.15,2.36))
energy_au_1e1, time_ps_1e1 = read_data('energy_evolution-rescaling-1e1') 
popt, pcov = curve_fit(decay_func, time_ps_1e1[5000:len(time_ps_1e1)], energy_au_1e1[5000:len(energy_au_1e1)])
plt.plot(time_ps_1e1[5000:len(time_ps_1e1)], energy_au_1e1[5000:len(energy_au_1e1)], 'b-', label=r'$\gamma$=10')
plt.plot(time_ps_1e1[5000:len(time_ps_1e1)], decay_func(time_ps_1e1[5000:len(time_ps_1e1)], popt[0], popt[1]), 'r--')
energy_au_2e1, time_ps_2e1 = read_data('energy_evolution-rescaling-2e1') 
popt, pcov = curve_fit(decay_func, time_ps_2e1[5000:len(time_ps_2e1)], energy_au_2e1[5000:len(energy_au_2e1)])
plt.plot(time_ps_2e1[5000:len(time_ps_2e1)], energy_au_2e1[5000:len(energy_au_2e1)], 'b-', label=r'$\gamma$=20')
plt.plot(time_ps_2e1[5000:len(time_ps_2e1)], decay_func(time_ps_2e1[5000:len(time_ps_2e1)], popt[0], popt[1]), 'r--')
energy_au_3e1, time_ps_3e1 = read_data('energy_evolution-rescaling-3e1') 
popt, pcov = curve_fit(decay_func, time_ps_3e1[5000:len(time_ps_3e1)], energy_au_3e1[5000:len(energy_au_3e1)])
plt.plot(time_ps_3e1[5000:len(time_ps_3e1)], energy_au_3e1[5000:len(energy_au_3e1)], 'b-', label=r'$\gamma$=30')
plt.plot(time_ps_3e1[5000:len(time_ps_3e1)], decay_func(time_ps_3e1[5000:len(time_ps_3e1)], popt[0], popt[1]), 'r--')
energy_au_4e1, time_ps_4e1 = read_data('energy_evolution-rescaling-4e1') 
popt, pcov = curve_fit(decay_func, time_ps_4e1[5000:len(time_ps_4e1)], energy_au_4e1[5000:len(energy_au_4e1)])
plt.plot(time_ps_4e1[5000:len(time_ps_4e1)], energy_au_4e1[5000:len(energy_au_4e1)], 'b-', label=r'$\gamma$=40')
plt.plot(time_ps_4e1[5000:len(time_ps_4e1)], decay_func(time_ps_4e1[5000:len(time_ps_4e1)], popt[0], popt[1]), 'r--')
energy_au_5e1, time_ps_5e1 = read_data('energy_evolution-rescaling-5e1') 
popt, pcov = curve_fit(decay_func, time_ps_5e1[5000:len(time_ps_5e1)], energy_au_5e1[5000:len(energy_au_5e1)])
plt.plot(time_ps_5e1[5000:len(time_ps_5e1)], energy_au_5e1[5000:len(energy_au_5e1)], 'b-', label=r'$\gamma$=50')
plt.plot(time_ps_5e1[5000:len(time_ps_5e1)], decay_func(time_ps_5e1[5000:len(time_ps_5e1)], popt[0], popt[1]), 'r--')
plt.xlabel('time [ps]')
plt.ylabel('energy [a.u.]')
#plt.title(f"decay rate: {1/tau:.4f} ps^-1")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('energy_evolution.png',dpi=300)
plt.show()
"""
plt.figure(figsize=(3.15, 2.36))

gammas = [10, 20, 30, 40, 50]
files = [
    "energy_evolution-rescaling-1e1",
    "energy_evolution-rescaling-2e1",
    "energy_evolution-rescaling-3e1",
    "energy_evolution-rescaling-4e1",
    "energy_evolution-rescaling-5e1",
]

colors = plt.cm.coolwarm(np.linspace(0.1, 0.9, len(gammas)))

for gamma, fname, color in zip(gammas, files, colors):
    energy_au, time_ps = read_data(fname)

    x = time_ps[12000:]
    y = energy_au[12000:]

    popt, pcov = curve_fit(decay_func, x, y)

    plt.plot(
        x, y*1e5,
        color=color,
        #linewidth=1.5,
        label=rf"$\gamma$={gamma}"
    )

    plt.plot(
        x,
        decay_func(x, popt[0], popt[1])*1e5,
        "--",
        color='black',
        #linewidth=1.2
    )
x_label = [0.4325, 0.4325001]
y_label = [5.127e-1, 5.127001e-1]
plt.plot(x_label, y_label, "--",color='black', label='exponential fitting')

#plt.text('ethylene')
plt.text(0.5, 1.10, "ethylene", transform=plt.gca().transAxes, fontsize=10, ha="center", va="top")
plt.xlabel("time [ps]")
#plt.ylabel("energy [a.u.]")
#plt.ylabel(r"$E - E_0$ [a.u.]")
plt.ylabel(r"$E - E_0$ [10$^{-5}$ a.u.]")
plt.xlim(0, 1.0)
plt.ylim(0, 1.2e0)
plt.legend(frameon=False, ncol=2, columnspacing=0.4)
plt.tight_layout()
plt.savefig("energy_evolution.png", dpi=300)
plt.show()
