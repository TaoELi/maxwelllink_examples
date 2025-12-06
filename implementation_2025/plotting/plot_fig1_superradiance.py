import numpy as np
import matplotlib.pyplot as plt
import columnplots as clp

def get_superradiance_data(path="../meep_tls_socket_upperbound/"):
    filename_1tls = "ntls_1/tls_se_data.npz"
    data_1tls = np.load(path + filename_1tls)
    time_1tls = data_1tls["time_au"]
    Pe_1tls = data_1tls["Pe"]

    filename_4tls = "ntls_4_norescaling/tls_se_data.npz"
    data_4tls = np.load(path + filename_4tls)
    time_4tls = data_4tls["time_au"]
    Pe_4tls = data_4tls["Pe"]

    return time_1tls, Pe_1tls, time_4tls, Pe_4tls

def get_time_statistics(path="../meep_tls_socket_upperbound/"):
    ntls_lst = [16, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    time_lst = []
    for ntls in ntls_lst:
        filename = path + f"ntls_{ntls}/time_count.txt"
        data = np.loadtxt(filename)
        time_per_step = np.mean(data)
        print("time per step for", ntls, "TLSs:", time_per_step, "s")
        time_lst.append(time_per_step)
    ntls_lst = np.array(ntls_lst)
    time_lst = np.array(time_lst) / ntls_lst
    return ntls_lst, time_lst

# ---- data for Fig. a ----
t_1tls, Pe_1tls, t_4tls, Pe_4tls = get_superradiance_data(path="../meep_tls_socket_upperbound/")

dipole_moment = 0.1 # meep units of mu12 when time units of meep is 0.1 fs
frequency = 1.0 # meep units of omega when time units of meep is 0.1 fs
# analytical golden-rule decay rate
gamma = dipole_moment**2 * frequency**2 / 2.0
# convert time from au to fs
time_fs = t_1tls * 0.02418884254
time_meep = time_fs / 0.1
Pe_analytical_1 = Pe_1tls[0] * np.exp(-gamma * time_meep)  # convert time to a.u.
Pe_analytical_4 = Pe_1tls[0] * np.exp(-4 * gamma * time_meep)  # 4 TLSs

x1s, y1s = [time_fs]*4, [Pe_1tls, Pe_4tls, Pe_analytical_1, Pe_analytical_4]


# ---- data for Fig. b ----
ntls_lst, time_lst = get_time_statistics(path="../meep_tls_socket_upperbound/")
x2s, y2s = [ntls_lst]*2, [time_lst, time_lst[0] * np.ones(np.size(time_lst))]  # ideal scaling line

fig, axes = clp.initialize(2, 1, width=4.3, height=4.3*0.618*2, return_fig_args=True,
                           fontsize=12, labelthem=True, labelthemPosition=[0.95, 0.95],
                           labelsize=14,
                           LaTeX=True)

clp.plotone(x1s[0:3], y1s[0:3], axes[0], labels=['1 TLS', '4 TLSs', "analytical"], 
            linestyles=['-', '-', ":"], 
            colors=[clp.red, clp.navy_blue, clp.black],
            xlabel='time [fs]', ylabel=r'population $P_{\rm e}$',
            xlim=[0,40], ylim=[0, 1.e-4], alphaspacing=0.2,
            yscientificAtLabel=True)
clp.plotone(x1s[3:4], y1s[3:4], axes[0],
            colors=["k:"],
            showlegend=False, alpha=0.6,
            xlabel='time [fs]', ylabel=r'population $P_{\rm e}$',
            xlim=[0,40], ylim=[0, 1.e-4],
            yscientificAtLabel=True)


clp.plotone(x2s, y2s, axes[1],
            colors=[clp.cyan, clp.darkgray],
            linestyles=['-', '--'],
            markers=['o', ""],
            labels=['MaxwellLink runtime', 'baseline'],
            ylabel=r'stepping time / \# drivers [s]',
            xlabel=r'\# drivers (16 drivers per CPU core)',
            xlog=True, ylog=False,
            ylim=[0, 1e-3], xlim=[10, 1e5],
            alphaspacing=0.2,
            showlegend=True)

clp.adjust(tight_layout=True, savefile="fig1_superradiance.pdf")
