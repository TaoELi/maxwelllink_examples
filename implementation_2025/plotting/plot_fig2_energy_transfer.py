import numpy as np
import matplotlib.pyplot as plt
import columnplots as clp
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from maxwelllink.tools import rt_tddft_spectrum, lr_tddft_spectrum

au_to_fs = 0.02418884254

# load data
def get_energy_transfer_data(path="../meep_energy_transfer/energy_transfer_tls_tls/"):

    mode = path.split("/")[-2].split("_")[-1]
    mode_write = mode
    if mode == "rttddft" or mode == "rtehrenfest" or mode == "nonres" or mode == "rtehrenfestfine":
        mode_write = "tddft"
    if mode == "180" or mode == "100" or mode == "30" or mode == "60" or mode == "90" or mode == "182":
        mode_write = "qutip"
    print("Loading data for mode:", mode)

    filename_donor = "donor_tls_data.npz"
    filename_acceptor = f"acceptor_{mode_write}_data.npz"
    
    data_donor = np.load(path + filename_donor)
    time_donor = data_donor["time_au"] * au_to_fs
    energy_donor = data_donor["energy_au"]

    data_acceptor = np.load(path + filename_acceptor)
    time_acceptor = data_acceptor["time_au"] * au_to_fs
    energy_acceptor = data_acceptor["energy_au"]
    
    try:
        energy_kin_acceptor = data_acceptor["energy_kin_au"]
    except:
        energy_kin_acceptor = np.zeros_like(energy_acceptor)

    energy_acceptor -= energy_acceptor[0]  # normalize to initial energy

    return time_donor, energy_donor, time_acceptor, energy_acceptor, energy_kin_acceptor

# frequency domain spectrum data
dt_rttddft_au = 0.04 
rt_data = np.loadtxt("../meep_energy_transfer_optimized/tddft_benchmark/rt_tddft_energy_dipoles_0.txt")
time_au, energy_au, mux_au, muy_au, muz_au = rt_data[:, 0], rt_data[:, 1], rt_data[:, 2], rt_data[:, 3], rt_data[:, 4]

poles = np.loadtxt("../meep_energy_transfer_optimized/tddft_benchmark/lrtddft_poles_au_182.txt")
oscillator_strengths = np.loadtxt("../meep_energy_transfer_optimized/tddft_benchmark/lrtddft_oscillator_strengths_182.txt")

nskip = 5
nend = int(len(time_au)*0.8)
mux, muy, muz = mux_au[:nend:nskip], muy_au[:nend:nskip], muz_au[:nend:nskip]
sp_tot = 0.0
for mu in [mux, muy, muz]:
    freq_ev, sp, _, _ = rt_tddft_spectrum(mu, dt_au=dt_rttddft_au*nskip, sp_form="absorption", e_start_ev=8, e_cutoff_ev=40.0, sigma=1e4, w_step=1e-4)
    sp_tot += sp

freq_ev_lr, sp_lr = lr_tddft_spectrum(poles, oscillator_strengths, e_cutoff_ev=40.0, linewidth=1e-2, w_step=1e-5)

xs_freq = [freq_ev_lr, freq_ev]
ys_freq = [sp_lr / np.max(sp_lr), sp_tot / np.max(sp_tot)]  
labels_freq = ['LR-TDDFT', 'RT-TDDFT']

# time-domain energy transfer data
t_d_tls, e_d_tls, t_a_tls, e_a_tls, e_kin_a_tls = get_energy_transfer_data(path="../meep_energy_transfer_optimized/energy_transfer_tls_tls/")
t_d_qutip_30, e_d_qutip_30, t_a_qutip_30, e_a_qutip_30, e_kin_a_qutip_30 = get_energy_transfer_data(path="../meep_energy_transfer_optimized/energy_transfer_tls_qutip_30/")
t_d_qutip_180, e_d_qutip_180, t_a_qutip_180, e_a_qutip_180, e_kin_a_qutip_180 = get_energy_transfer_data(path="../meep_energy_transfer_optimized/energy_transfer_tls_qutip_182/")
t_d_rttddft, e_d_rttddft, t_a_rttddft, e_a_rttddft, e_kin_a_rttddft = get_energy_transfer_data(path="../meep_energy_transfer_optimized/energy_transfer_tls_rttddft/")
t_d_rtehrenfest, e_d_rtehrenfest, t_a_rtehrenfest, e_a_rtehrenfest, e_kin_a_rtehrenfest = get_energy_transfer_data(path="../meep_energy_transfer_optimized/energy_transfer_tls_rtehrenfest/")
xs = [t_a_rttddft, t_a_rtehrenfest, t_a_tls, t_a_qutip_30, t_a_qutip_180]
ys = [e_a_rttddft, e_a_rtehrenfest, e_a_tls, e_a_qutip_30, e_a_qutip_180]
ys = [y * 1.e6 for y in ys]  # scale to 1e-6 a.u.

labels = ['RT-TDDFT', 'RT-Eh', 'TLS', 'QuTiP (30 states)', 'QuTiP (182 states)']
fig, axes = clp.initialize(2, 1, width=5.3, height=5.3*0.618*2, return_fig_args=True,
                         fontsize=12, LaTeX=True, labelthem=True,
                         labelthemPosition=[0.07, 0.97], labelsize=14)

# plot linear-response spectrum as the first subplot
ax = axes[0]
img = mpimg.imread("hcn.png")
imagebox = OffsetImage(img, zoom=0.15)
ab = AnnotationBbox(imagebox, (0.18, 0.62), xycoords='axes fraction', frameon=False, zorder=1)
ax.add_artist(ab)

clp.plotone(xs_freq, ys_freq, ax,
            labels=labels_freq,
            xlabel='energy [eV]', ylabel='absorption strength [arb. units]',
            linestyles=['-', ':'],
            xlim=[8, 30], 
            colors=[clp.black, clp.red],
            legendFontSize=8,
            showlegend=True) 

au_to_ev = 27.211386245988
# add an arror to point at this peak
ax.annotate('strongest excitation', xy=(0.4932*au_to_ev, 0.7), xytext=(0.65*au_to_ev, 0.8),
            arrowprops=dict(facecolor='black', arrowstyle='->'),
            fontsize=10, ha='center')   

# add inset image to demonstrate energy transfer setup
ax = axes[1]
img = mpimg.imread("energy_transfer_demo.png")
imagebox = OffsetImage(img, zoom=0.15)
ab = AnnotationBbox(imagebox, (0.18, 0.82), xycoords='axes fraction', frameon=False, zorder=1)
ax.add_artist(ab)

clp.plotone(xs, ys, ax,
            labels=labels,
            xlabel='time [fs]', ylabel=r'acceptor energy gain [$\times 10^{-6}$ a.u.]',
            linestyles=['-', '-.', '-', '--', ":"],
            xlim=[0, 4], ylim=[0, 2.8],
            colors=[clp.red,  clp.darkgray, clp.black, clp.cyan, clp.dark_blue],
            yscientificAtLabelString="a.u.",
            legendFontSize=8, alphaspacing=0.01,
            legendloc="upper right",
            showlegend=True) 

# add a color shadow region between two curves to highlight the kinetic energy contribution
ax.fill_between(t_a_rtehrenfest, (e_a_rtehrenfest - e_kin_a_rtehrenfest)*1.e6, e_a_rtehrenfest*1.e6, color=clp.darkgray, alpha=0.3)
ax.text(3.52, 1.55, 'K.E.', color=clp.darkgray, fontsize=10)
# add a curved arrow between text and shadow region
ax.annotate('', xy=(3.6, 1.5), xytext=(3.5, 1.08),
            arrowprops=dict(facecolor=clp.darkgray, arrowstyle='->', connectionstyle="arc3,rad=-0.3", color=clp.darkgray),
            fontsize=10, ha='center')

clp.adjust(tight_layout=True, savefile="fig2_energy_transfer.pdf")