import h5py
from matplotlib.patches import Rectangle
import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
from maxwelllink.tools import ir_spectrum

def smooth(x,window_len=11,window='hamming'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
	x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
	the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[window_len//2-1:-window_len//2]

def plot_2d_background():

    fig, axes = clp.initialize(5, 3, width=3*4.3, height=4.3*5*0.86, LaTeX=True, fontsize=12, return_fig_args=True)
    efield_data_list = ["./data/mxl_144*144/2d_144_144_efield_neq_0.002.npy", "./data/mxl_144*144/2d_144_144_efield_neq_0.005.npy"]
    coenergy_data_list = ["./data/mxl_144*144/coenergy_36_0.002.npy", "./data/mxl_144*144/coenergy_36_0.005.npy"]
    time_list = [50, 150, 300]
    time_labels = [r"$t=%.1f$ ps" % (t/50) for t in time_list]
    label_list = {
        1: ["(c)","(d)","(e)"],
        2: ["(f)","(g)","(h)"],
        3: ["(i)","(j)","(k)"],
        4: ["(l)","(m)","(n)"],
    }
    amp_list = [r"$E_0=2\times10^{-3}$ a.u.", r"$E_0=5\times10^{-3}$ a.u."]
    row_specs = [
        ("efield", 0, "E-field intensity [a.u.]"),
        ("coenergy", 0, r"$E_{\rm C=O}$ per molecule [$10^2$ cm$^{-1}$]"),
        ("efield", 1, "E-field intensity [a.u.]"),
        ("coenergy", 1, r"$E_{\rm C=O}$ per molecule [$10^2$ cm$^{-1}$]"),
    ]

    axes[0, 0].axis("off")
    axes[0, 2].axis("off")

    data_cache = {}
    scale_cache = {}
    for kind, amp_index, _ in row_specs:
        key = (kind, amp_index)
        if key in data_cache:
            continue
        if kind == "efield":
            data = np.load(efield_data_list[amp_index])
            ref_sp = (np.abs(data[time_list[0], :]).reshape(144, 144))**2
            vmax = np.max(np.max(ref_sp))
            vmin = vmax * 0.005
        else:
            data = np.load(coenergy_data_list[amp_index]) * (219474.63 / 36) * 0.01
            vmin = np.percentile(data[time_list[0]], 98.5)
            vmax = np.percentile(data[time_list[0]], 99.9)
        data_cache[key] = data
        scale_cache[key] = (vmin, vmax)

    for row, (kind, amp_index, cbar_label) in enumerate(row_specs, start=1):
        data = data_cache[(kind, amp_index)]
        vmin, vmax = scale_cache[(kind, amp_index)]
        last_pos = None
        for col, time_index in enumerate(time_list):
            if kind == "efield":
                sp = (np.abs(data[time_index, :]).reshape(144, 144))**2
            else:
                sp = data[time_index]
            extent = [0 , 144, 0, 144]
            pos = axes[row, col].imshow(sp, aspect='equal', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            last_pos = pos
            axes[row, col].set_box_aspect(1)
            clp.plotone([], [], axes[row, col], colors=["c-","c-"], ylabel=r"$L_y$ \ position \ [$\mu\rm{m}$]" if col==0 else None, xlabel=r"$L_x$ \ position \ [$\mu\rm{m}$]" if row==4 else None, showlegend=False)
            axes[row, col].set_yticks([0,72, 144])
            axes[row, col].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[row, col].set_xticks([0,72, 144])
            axes[row, col].set_xticklabels([r"$0$", r"$200$", r"$400$"])
            if col != 0:
                axes[row, col].set_yticklabels([])
            if row != 4:
                axes[row, col].set_xticklabels([])
            axes[row, col].text(0.98, 0.18, time_labels[col], transform=axes[row, col].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[row, col].text(0.02, 0.98, label_list[row][col], transform=axes[row, col].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[row, col].text(0.98, 0.08, amp_list[amp_index], transform=axes[row, col].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            if col == 0:
                axes[row, col].text(0.55, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[row, col].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
                axes[row, col].text(0.55, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=axes[row, col].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")

            marker_center = (85, 72)
            marker_width = 20
            marker_height = 24
            marker_linewidth = 1.5
            center_x, center_y = marker_center
            center_marker = Rectangle(
                (center_x - marker_width / 2, center_y - marker_height / 2),
                marker_width,
                marker_height,
                fill=False,
                edgecolor="cyan",
                linewidth=marker_linewidth,
                zorder=3,
            )
            axes[row, col].add_patch(center_marker)
        cbar_ax = axes[row, 2].inset_axes([1.04, 0.0, 0.055, 1.0])
        cbar = fig.colorbar(last_pos, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(cbar_label, fontsize=12)
        if kind == "coenergy":
            if amp_index == 0:
                cbar.set_ticks([4, 5, 6, 7])
                cbar.set_ticklabels([r"$4$", r"$5$", r"$6$", r"$7$"])
            else:
                cbar.set_ticks([8, 9, 12, 16])
                cbar.set_ticklabels([r"$8$", r"$9$", r"$12$", r"$15$"])
        cbar.ax.yaxis.set_label_coords(4.2, 0.5)

    qc_file = np.load("./data/mxl_144*144/qc_y.npy", allow_pickle=True).item()
    sp = qc_file['sp']
    x  = qc_file['freq']
    dx = x[2] - x[1]
    nstart, nend = int(2200 / dx), int(3000 / dx)
    x = x[nstart:nend]
    sp = np.abs(sp[nstart:nend,:]) / 1e32
    sp = sp[::-1, :]
    freq_cav_inplane_min = 12.5
    freq_cav_inplane_max = 12.5 * 144
    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    vmax = np.max(np.max(sp))
    vmin = vmax * 0.001
    pos = axes[0,1].imshow(sp, aspect='auto', extent=extent,
            cmap=cm.inferno,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, axes[0,1], showlegend=False, colors=["orange", "c"], linestyles=["--", "--"], lw=1.2, xlim=[12.5,1800],
            xlabel=r"$\omega_{\parallel}$ [$\rm{cm}^{-1}$]",
            ylabel=r"IR frequency [$\rm{cm}^{-1}$]")
    axes[0,1].text(1190, 2550, "cavity photon", color='c', fontsize=12)
    axes[0,1].text(940, 2370, "C=O asym. stretch", color='orange', fontsize=12)
    axes[0,1].tick_params(color='c', labelsize='medium', width=2)

    axes[0,1].set_xticks([12.5,600,1200,1800])
    axes[0,1].set_yticks([2200,2400,2600,2800])
    axes[0,1].text(0.02, 0.98, "(b)", transform=axes[0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
    axes[0,1].annotate('', xy=(450, 2450), xytext=(450, 2550), arrowprops=dict(facecolor='w', edgecolor='w', arrowstyle='->', alpha=0.8, lw=2), fontsize=20)
    axes[0,1].text(0.38, 0.5, "UP excitation", transform=axes[0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")

    plt.rcParams["axes.axisbelow"] = False
    fig.subplots_adjust(left=0.10, right=0.90, bottom=0.06, top=0.97, wspace=0.032, hspace=0.112)
    fig_width, fig_height = fig.get_size_inches()
    first_row_pos = axes[0, 1].get_position()
    reference_pos = axes[1, 0].get_position()
    cbar_reference_pos = axes[1, 2].get_position()
    first_second_gap_boost = reference_pos.height * 0.04
    for row in range(1, 5):
        for col in range(3):
            pos_to_shift = axes[row, col].get_position()
            axes[row, col].set_position([
                pos_to_shift.x0,
                pos_to_shift.y0 - first_second_gap_boost,
                pos_to_shift.width,
                pos_to_shift.height,
            ])
    cbar_gap = cbar_reference_pos.width * 0.04
    cbar_width = cbar_reference_pos.width * 0.055
    dispersion_height = reference_pos.height
    scaled_square_width = dispersion_height * fig_height / fig_width
    dispersion_width = dispersion_height * fig_height / fig_width / (0.618 * 1.1)
    inset_width = scaled_square_width * 2400 / (2400 - 397)
    first_row_image_gap = scaled_square_width * 0.24
    first_row_left = axes[1, 0].get_position().x0
    first_row_right = axes[1, 2].get_position().x1 + cbar_gap + cbar_width
    first_row_group_width = inset_width + first_row_image_gap + dispersion_width + cbar_gap + cbar_width
    first_row_group_left = first_row_left + (first_row_right - first_row_left - first_row_group_width) / 2
    b_left = first_row_group_left + inset_width + first_row_image_gap
    axes[0, 1].set_position([
        b_left,
        first_row_pos.y0 + (first_row_pos.height - dispersion_height) / 2 + first_second_gap_boost * 0.5,
        dispersion_width,
        dispersion_height,
    ])
    b_pos = axes[0, 1].get_position()
    cbar_left = b_pos.x1 + cbar_gap
    cbar_ax = fig.add_axes([
        cbar_left,
        b_pos.y0,
        cbar_width,
        b_pos.height,
    ])
    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("spectral intensity [arb. units]", fontsize=12)
    clp.adjust(savefile=f"./figures/fig1_2d_efield_coenergy_background.png")

def plot_2d_final():

    plot_2d_background()
    from pathlib import Path
    from PIL import Image

    base_file = Path("./figures/fig1_2d_efield_coenergy_background.png")
    inset_file = Path("./figures/fig1a.png")
    output_file = Path("./figures/fig1_2d_final.png")

    base = Image.open(base_file).convert("RGBA")
    inset = Image.open(inset_file).convert("RGBA")
    base_rgb = np.asarray(base.convert("RGB"))
    height, width = base_rgb.shape[:2]

    b_x0 = int(width * 0.35)
    b_x1 = int(width * 0.90)
    b_y1 = int(height * 0.20)
    b_region = base_rgb[:b_y1, b_x0:b_x1]
    b_dark = b_region.mean(axis=2) < 22
    b_rows = np.where(b_dark.sum(axis=1) > b_dark.shape[1] * 0.35)[0]
    if b_rows.size == 0:
        raise RuntimeError("Could not locate panel (b) height in the background image.")
    paste_top = int(b_rows[0])
    paste_height = int(b_rows[-1] - b_rows[0] + 1)

    b_cols = np.where(b_dark.sum(axis=0) > b_dark.shape[0] * 0.65)[0]
    if b_cols.size == 0:
        raise RuntimeError("Could not locate panel (b) left edge in the background image.")
    b_left = int(b_x0 + b_cols[0])

    resampling = getattr(Image, "Resampling", Image).LANCZOS
    top_crop = 397
    inset = inset.crop((0, top_crop, inset.width, inset.height))
    paste_width = int(round(inset.width * paste_height / inset.height))
    inset = inset.resize((paste_width, paste_height), resampling)
    paste_left = int(round(b_left - paste_height * 0.24 - paste_width))

    final = base.copy()
    final.alpha_composite(inset, dest=(paste_left, paste_top))
    output_file.parent.mkdir(parents=True, exist_ok=True)

    dpi = 300
    fig = plt.figure(figsize=(base.width / dpi, base.height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(final)
    ax.axis("off")
    fig.text((paste_left + paste_width * 0.98) / base.width, 1 - (paste_top + paste_height * 0.02) / base.height,
             "(a)", fontsize=12, fontweight='bold', color="w", ha="right", va="top")
    fig.savefig(output_file, dpi=dpi, pad_inches=0)
    plt.close(fig)
    print(f"Saved final figure to {output_file}")

def plot_1d_efield():

    axes = clp.initialize(3, 2, width=2*4.3, height=4.3*0.618*1.1*3, LaTeX=True, fontsize=12)
    uplist = [36, 24]
    k_parallel_list = [450, 300]
    coeff = np.array([0.11638023, 0.07069759])
    wc_up = np.array([0.61, 0.55])
    length_list = [2600]*2
    vg = np.array([r"$%.3f\,c$"%(coeff[i]) for i in range(2)])
    label_list = {0 : ["(a)","(b)"], 1: ["(c)","(d)"], 2: ["(e)","(f)"]}
    amp_dict = {0 : ["0.001", "0.0024", "0.003", "0.005", "0.007"], 1: ["0.001", "0.0016", "0.003", "0.005", "0.007"]}
    for j in range(2):
        amp_list = [r"$%.1f\times10^{-3}$" % (float(i)*1000) for i in amp_dict[j]]
        sp_list = [np.load(f"./data/mxl_144_efield/effective_efield_{uplist[j]}_{amp_dict[j][i]}.npy")[50:2551].T for i in range(len(amp_dict[j]))]
        mt_list, mmsd_list = get_adaptive_emsd(uplist[j], coeff[j], amp_list=amp_dict[j], sp_list=sp_list)
        clp.plotone(mt_list, mmsd_list, axes[2,j], ylabel=r"$\mathrm{MSD}$ \ [$10^{4} \ \mu\rm{m}^2$]" if j == 0 else None, xlabel="time [ps]",
                    showlegend=True, legendloc=(0.05,0.3), legendFontSize=8, 
                    xlim=[0,mt_list[0][-1]], ylim=[0, 0.05] if j==1 else [0, 0.2],
                    labels=amp_list, colorMap=plt.cm.hot, colorMap_endpoint=0.6, alpha=0.5)
        quad_t = mt_list[0][(mt_list[0] >= 1.0) & (mt_list[0] <= 4.0)]
        linear_t = mt_list[0][(mt_list[0] >= 2.0) & (mt_list[0] <= 5.0)]
        quad_fit_params = [[0.010412225466983087, -0.009692823652307786, 0.002930743601640022],
                           [0.0034124054608374583, -0.003780895807991756, 0.001182095496525906]]
        linear_fit_params = [[0.0659900795219766, -0.13787260151910052],
                             [0.01722295767678947, -0.03394572375462868]]
        quad_sampled_t = quad_fit_params[j][0] * quad_t**2 + quad_fit_params[j][1] * quad_t + quad_fit_params[j][2]
        linear_sampled_t = linear_fit_params[j][0] * linear_t + linear_fit_params[j][1]
        clp.plotone([quad_t], [quad_sampled_t], axes[2,j], colors=["g--"], showlegend=False, lw=1.5)
        clp.plotone([linear_t], [linear_sampled_t], axes[2,j], colors=["b--"], showlegend=False, lw=1.5)
        axes[2,j].text(0.02, 0.98, label_list[2][j], transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,j].text(0.68, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[j]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[2,j].text(0.32, 0.85, r"$W_{\rm ph}=$"+rf"$ \ {wc_up[j]} \ $", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,j].axvline(x=1, linestyle='-.', alpha=0.3, color="c")
        axes[2,j].text(0.1, 0.9, "pulse \n on", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', alpha=0.3, color="c")
        axes[2,j].text(0.22, 0.18, r"$\propto t^2$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="g")
        axes[2,j].text(0.85, 0.65, r"$\propto t$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="b")
    
    axes[2,1].set_yticks([0, 0.025, 0.05])
    axes[2,1].set_yticklabels([r"$0$", r"$0.4$", r"$0.8$"])
    axes[2,0].get_legend().set_title(r"$E_0$ [a.u.]")
    axes[2,1].get_legend().set_title(r"$E_0$ [a.u.]")
    axes[2,0].set_yticks([0, 0.05, 0.1, 0.15, 0.2])
    axes[2,0].set_yticklabels([r"$0$", r"$0.8$", r"$1.6$", r"$2.4$", r"$3.2$"])

    amp = [0.001, 0.007]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    for y0 in range(2):
        ref = np.load(f"./data/mxl_144_efield/effective_efield_24_{amp[y0]}.npy")
        vmax = np.max(np.max(ref))
        vmin = vmax * 0.01 if y0 == 0 else vmax * 0.005
        for x0 in range(2):
            sp_avg = np.load(f"./data/mxl_144_efield/effective_efield_{uplist[x0]}_{amp[y0]}.npy")
            ntimes = sp_avg.shape[0]
            sp_avg = np.abs(sp_avg).T
            extent = [0 , 20, 0, 144]
            pos = axes[y0, x0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            t = np.linspace(0,(ntimes-1)*2,ntimes)[:length_list[x0]]/1000
            xs = [t]
            ys = [144*0.2+144*0.75*coeff[x0]*t]
            clp.plotone(xs, ys, axes[y0, x0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if x0==0 else None, xlabel="time [ps]" if y0==1 else None, showlegend=False, lw=1.5)
            axes[y0, x0].set_yticks([0,72, 144])
            axes[y0, x0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[y0, x0].text(0.6, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[x0]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.6, 0.85, r"$W_{\rm ph}=$"+rf"$ \ {wc_up[x0]} \ $", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.98, 0.1, amp_list[y0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[x0, y0].text(0.02, 0.98, label_list[x0][y0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[x0, y0].text(0.02, 0.1, r"$\tilde{v}_g \ = \ $"+f"{vg[y0]}", transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="c")
            cbar_ax = axes[y0, x0].inset_axes([1.04, 0.0, 0.055, 1.0])
            cbar = axes[y0, x0].figure.colorbar(pos, cax=cbar_ax)
            cbar.ax.tick_params(labelsize=12)
            cbar.set_label("E-field intensity [a.u.]", fontsize=12)
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/fig2_1d_efield.png")

def get_adaptive_emsd(up, vg, amp_list, folder="./data/mxl_144_scanamp/", ngrid=144,
                      reference_ngrid=144, scale_emsd_by_size=True, sp_list=None):

    def smooth_nan(x, window_len=151, window="hamming"):
        x = np.asarray(x, dtype=float)
        valid = np.isfinite(x)
        if window_len < 3:
            return x
        if window_len % 2 == 0:
            window_len += 1
        if window == "flat":
            w = np.ones(window_len, dtype=float)
        else:
            w = getattr(np, window)(window_len)
        w /= np.sum(w)
        values = np.where(valid, x, 0.0)
        norm = np.convolve(valid.astype(float), w, mode="same")
        y = np.convolve(values, w, mode="same") / np.maximum(norm, 1e-12)
        y[norm < 0.1] = np.nan
        return y

    def ballistic_guided_centroid(sp, t, x_grid, x0, vtilde, sigma=0.08, search_radius=0.24):
        x_mean = np.full_like(t, np.nan, dtype=float)
        prev = x0
        for n, tn in enumerate(t):
            col = sp[:, n]
            ref = x0 + 0.75 * vtilde * tn
            mask = np.abs(x_grid - ref) <= search_radius
            if not mask.any() or ref > 1.05:
                continue
            idxs = np.where(mask)[0]
            peak = x_grid[idxs[np.nanargmax(col[mask])]]
            soft_window = np.exp(-0.5 * ((x_grid - peak) / sigma) ** 2)
            floor = 0.003 * np.nanmax(col[mask])
            weight = np.maximum(col - floor, 0.0) * soft_window
            den = np.sum(weight)
            if den <= 0:
                x_mean[n] = prev
                continue
            value = np.sum(x_grid * weight) / den
            if n > 0 and np.isfinite(prev):
                value = prev + np.clip(value - prev, -0.018, 0.018)
            x_mean[n] = value
            prev = value
        return smooth_nan(x_mean, window_len=181)

    def localization_aware_centroid(sp, t, x_grid, x0, vtilde, sigma=0.038, search_radius=0.14, t_ballistic=1.1, max_step=0.006, forward_bias=0.15):
        x_mean = np.full_like(t, np.nan, dtype=float)
        prev = x0
        min_signal = np.nanmax(sp) * 2e-4
        for n, tn in enumerate(t):
            col = sp[:, n]
            if np.nanmax(col) < min_signal:
                x_mean[n] = prev
                continue
            ballistic_ref = x0 + 0.75 * vtilde * tn
            ref = ballistic_ref if n == 0 or tn <= t_ballistic else prev
            mask = np.abs(x_grid - ref) <= search_radius
            if not mask.any():
                x_mean[n] = prev
                continue
            idxs = np.where(mask)[0]
            score = col[idxs] * (1.0 + forward_bias * np.maximum(x_grid[idxs] - prev, 0.0) / search_radius)
            peak = x_grid[idxs[np.nanargmax(score)]]
            soft_window = np.exp(-0.5 * ((x_grid - peak) / sigma) ** 2)
            floor = 0.015 * np.nanmax(col[mask])
            weight = np.maximum(col - floor, 0.0) * soft_window
            den = np.sum(weight)
            if den <= 0:
                x_mean[n] = prev
                continue
            value = np.sum(x_grid * weight) / den
            if n > 0 and np.isfinite(prev):
                value = prev + np.clip(value - prev, -max_step, max_step)
            x_mean[n] = value
            prev = value
        return smooth_nan(x_mean, window_len=181)

    def keep_before_reflection(t, displacement, min_time=2.0, min_drop=0.03):
        valid = np.isfinite(displacement)
        if np.sum(valid) < 3:
            return valid
        idx = np.where(valid)[0]
        d = smooth_nan(displacement, window_len=101)
        d_valid = d[idx]
        t_valid = t[idx]
        peak_so_far = np.maximum.accumulate(d_valid)
        drop_from_peak = peak_so_far - d_valid
        reflected = (t_valid > min_time) & (drop_from_peak > min_drop)
        if not np.any(reflected):
            return valid
        first_reflected = np.where(reflected)[0][0]
        peak_before_reflection = np.argmax(d_valid[:first_reflected + 1])
        cutoff = idx[peak_before_reflection]
        keep = valid.copy()
        keep[cutoff + 1:] = False
        return keep

    mt_list, mmsd_list = [], []
    if sp_list is not None:
        if len(sp_list) != len(amp_list):
            raise ValueError("sp_list must have the same length as amp_list.")
    for j in range(len(amp_list)):
        mtimes = None
        if folder == "./data/mxl_n_dependence/" :
            if sp_list is not None:
                sp = sp_list[j]
            else:
                row_data = np.load(f"{folder}/multimode_cavmd_{ngrid}_{up}_{amp_list[j]}_neq.npy")
                sp = np.mean(np.reshape(row_data, (10000, 144, -1)), axis=2)**2
                sp = sp[50:2551, :].T
        else :
            if sp_list is not None:
                sp = sp_list[j]
            else:
                with h5py.File(f"{folder}/multimode_cavmd_{ngrid}_{up}_{amp_list[j]}_neq.h5", "r") as f:
                    #sp = np.abs(f["effective_efield"][50:2551, :, 1]).T**2
                    sp = np.mean(np.reshape(f["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2
                    sp = sp[50:2551, :].T
        actual_ngrid = sp.shape[0]
        x_grid_1d = np.arange(1, actual_ngrid + 1) / (actual_ngrid + 1)
        x0 = x_grid_1d[np.argmin(np.abs(x_grid_1d - 0.2))]
        mtimes = 2 * (sp.shape[1] - 1) / 1000
        t = np.linspace(0, mtimes, sp.shape[1])
        if amp_list[j] in ["0.001"]:
            x_mean = ballistic_guided_centroid(sp, t, x_grid_1d, x0, vg)
        else:
            x_mean = localization_aware_centroid(sp, t, x_grid_1d, x0, vg)
        finite = np.isfinite(x_mean)
        if not np.any(finite):
            continue
        first = np.where(finite)[0][0]
        displacement = x_mean - x_mean[first]
        if amp_list[j] in ["0.001"]:
            finite &= keep_before_reflection(t, displacement)
        size_scale = actual_ngrid / reference_ngrid if scale_emsd_by_size else 1.0
        MMSD = smooth_nan(displacement**2, window_len=101) * size_scale**2
        positive = finite & (t > 0) & np.isfinite(MMSD) & (MMSD > 0)
        mt_list.append(t[positive])
        mmsd_list.append(MMSD[positive])
    return mt_list, mmsd_list

def plot_mmsd():

    axes = clp.initialize(2, 2, width=2*4.3, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12)
    up_list = np.array([36, 42, 48, 54, 60, 66])
    vg_up = np.array([0.11638023, 0.14234382, 0.17112003, 0.20126065, 0.23186782, 0.26334189])
    amp = ["0.001", "0.007"]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    label_list = {0 : ["(a)","(b)"], 1: ["(c)","(d)"]}
    x_grid_1d = np.arange(1, 145) / 145
    evaluate_range = np.where((x_grid_1d >= 0.05) & (x_grid_1d <= 0.95))[0]
    for y0 in range(2):
        k_list, mean_list, std_list = [], [], []
        for idx, up in enumerate(up_list):
            sp = np.load(f"./data/mxl_144_coenergy/coenergy_{up}_{amp[y0]}.npy")
            sp = np.clip(sp - sp[0, None], 0.0, None)
            t_ps = np.arange(sp.shape[0]) * 2 / 1000
            time_window = (t_ps >= 10) & (t_ps <= 20)
            weight = sp[np.ix_(time_window, evaluate_range)]
            denominator = np.sum(weight, axis=1)
            valid = denominator > 0
            if not np.any(valid):
                continue
            weight = weight[valid]
            denominator = denominator[valid]
            x_eval = x_grid_1d[evaluate_range]
            x_mean_t = np.sum(x_eval[None, :] * weight, axis=1) / denominator
            x_var_t = np.sum((x_eval[None, :] - x_mean_t[:, None]) ** 2 * weight, axis=1) / denominator
            x_mean = np.nanmean(x_mean_t)
            x_var = np.nanmean(x_var_t)
            k_list.append(up * 12.5)
            mean_list.append(x_mean * 144) # match the number of grid points
            std_list.append(np.sqrt(x_var) * 144) # match the number of grid points

        k_parallel = np.array(k_list)
        msd_array = np.array(mean_list)
        variance_array = np.array(std_list)
        lower = np.maximum(msd_array - variance_array, 0.0)
        upper = msd_array + variance_array
        axes[y0,1].fill_between(k_parallel, lower, upper, color="0.75", alpha=0.5, linewidth=0)
        vg_dot_tf = (0.2 + vg_up * 0.75 * 2) * 144
        clp.plotone([k_parallel]*2, [msd_array, vg_dot_tf], axes[y0,1], colors=["ro-", "c--"], showlegend=True, legendloc=(0.02,0.7), legendFontSize=8, labels=["MPL", r"$\tilde{v}_{\rm g} \cdot t_{\rm f}$"],
                    xlim=[425, 850],
                    ylabel=r"most probable location [$\mu$m]", xlabel=r"$k_{\parallel}$ [$\rm{cm}^{-1}$]" if y0==1 else None)
        axes[y0,1].text(0.08, 0.97, label_list[y0][1], transform=axes[y0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[y0,1].text(0.05, 0.15, "average over time 10-20 ps", transform=axes[y0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[y0,1].set_yticks([0,72, 144])
        axes[y0,1].set_yticklabels([r"$0$", r"$200$", r"$400$"])

    for y0 in range(2):
        sp_avg = np.load(f"./data/mxl_144_coenergy/coenergy_36_{amp[y0]}.npy")
        ntimes = sp_avg.shape[0]
        sp_avg = np.abs(sp_avg).T * (219474.63 / 36) * 10**(-(y0+2)) # convert to cm^-1 and scale for better visualization
        vmin = np.percentile(sp_avg, 60)
        vmax = np.percentile(sp_avg, 99)
        extent = [0 , 20, 0, 144]
        pos = axes[y0,0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                cmap=cm.hot, interpolation='nearest',
                norm=LogNorm(vmin=vmin, vmax=vmax))
        clp.plotone([], [], axes[y0,0], colors=["c-","c-"], ylabel=r"$L_x$ position [$\mu$m]", xlabel="time [ps]" if y0==1 else None, showlegend=False)
        axes[y0,0].set_yticks([0,72, 144])
        axes[y0,0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
        axes[y0,0].text(0.1, 0.97, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[y0,0].transAxes, fontsize=13, fontweight='bold', va='top', ha='left', color="#00e5ff")
        axes[y0,0].text(0.98, 0.97, amp_list[y0], transform=axes[y0,0].transAxes, fontsize=13, fontweight='bold', va='top', ha='right', color="#00e5ff")
        axes[y0,0].text(0.08, 0.97, label_list[y0][0], transform=axes[y0,0].transAxes, fontsize=13, fontweight='bold', va='top', ha='right', color="#00e5ff")
        cbar_ax = axes[y0,0].inset_axes([1.04, 0.0, 0.055, 1.0])
        cbar = axes[y0,0].figure.colorbar(pos, cax=cbar_ax)
        if y0 == 0:
            cbar.set_ticks([2.3, 2.4, 2.5, 2.6, 2.7, 2.8])
            cbar.set_ticklabels([r"$2.3$", r"$2.4$", r"$2.5$", r"$2.6$", r"$2.7$", r"$2.8$"])
        else:
            cbar.set_ticks([1, 2, 3, 4])
            cbar.set_ticklabels([r"$1.0$", r"$2.0$", r"$3.0$", r"$4.0$"])
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(r"$E_{\mathrm{C=O}}$ per molecule [$10^{%d} \ $cm$^{-1}$]"%(y0+2), fontsize=12)
        axes[y0,0].hlines(144*0.05, 0, 20, colors="c", linestyles="--", lw=1.5)
        axes[y0,0].hlines(144*0.35, 0, 20, colors="c", linestyles="--", lw=1.5)

    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/fig3_mmsd.png")

def plot_neq_spectrum_with_model():
    axes = clp.initialize(3, 2, width=4.3*2, height=4.3*0.618*1.1*3, LaTeX=True, fontsize=12)
    molecule_color = "#ffb000"
    photon_color = "#00d5ff"
    annotation_color = "w"
    inset_color = "w"
    label_list = {0 : ["(a)","(c)","(e)"], 1 : ["(b)","(d)","(f)"]}
    amp_list = [r"$E_0=1\times 10^{-3}$ a.u.", r"$E_0=7\times 10^{-3}$ a.u."]
    amp = [0.001, 0.007]
    dt_fs = 2
    N = 144
    domega = 50*36/N
    common_xlim = [12.5, 1800]
    common_ylim = [2200, 2980]
    xticks = [12.5, 600, 1200, 1800]
    yticks = [2200, 2400, 2600, 2800]

    def crop_spectrum_to_window(spectrum, extent):
        x_axis = np.linspace(extent[0], extent[1], spectrum.shape[1])
        y_axis = np.linspace(extent[3], extent[2], spectrum.shape[0])
        x_keep = (x_axis >= common_xlim[0]) & (x_axis <= common_xlim[1])
        y_keep = (y_axis >= common_ylim[0]) & (y_axis <= common_ylim[1])
        cropped = spectrum[np.ix_(y_keep, x_keep)]
        cropped_extent = [x_axis[x_keep][0], x_axis[x_keep][-1], y_axis[y_keep][-1], y_axis[y_keep][0]]
        return cropped, cropped_extent

    def crop_spectrum_to_window_origin_lower(spectrum, extent):
        x_axis = np.linspace(extent[0], extent[1], spectrum.shape[1])
        y_axis = np.linspace(extent[2], extent[3], spectrum.shape[0])
        x_keep = (x_axis >= common_xlim[0]) & (x_axis <= common_xlim[1])
        y_keep = (y_axis >= common_ylim[0]) & (y_axis <= common_ylim[1])
        cropped = spectrum[np.ix_(y_keep, x_keep)]
        cropped_extent = [x_axis[x_keep][0], x_axis[x_keep][-1], y_axis[y_keep][0], y_axis[y_keep][-1]]
        return cropped, cropped_extent

    for x0 in range(2):
        sim_ax = axes[0,x0]
        energy_ax = axes[1,x0]
        qc_file = np.load(f"./data/mxl_144_qc/qc_spectra_36_{amp[x0]}.npy", allow_pickle=True).item()
        x = qc_file["freq"]
        sp = qc_file["sp"]
        dx = x[1] - x[0]
        nstart, nmid, nend = int(2200 / dx), int(2350 / dx), int(3000 / dx)
        x = x[nstart:nend]
        sp = np.abs(sp[nstart:nend,:]) / 1e34
        sp = sp[::-1, :]
        freq_cav_inplane_min = domega
        freq_cav_inplane_max = domega * N
        extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]
        sp, extent = crop_spectrum_to_window(sp, extent)

        vmax = np.max(np.max(sp))
        vmin = vmax * 0.0001
        pos = sim_ax.imshow(sp, aspect='auto', extent=extent,
                cmap=cm.inferno,
                interpolation='nearest',
                norm=LogNorm(vmin=vmin, vmax=vmax)
                )
        cbar_ax = sim_ax.inset_axes([1.04, 0.0, 0.055, 1.0])
        cbar = sim_ax.figure.colorbar(pos, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label("spectral intensity [arb. units]", fontsize=12)

        freq_cav_inplane = np.linspace(common_xlim[0], common_xlim[1], 400)

        xs = [freq_cav_inplane]*2
        ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2)**0.5]
        clp.plotone(xs, ys, sim_ax, showlegend=False, colors=[molecule_color, photon_color], linestyles=["--", "--"], lw=1.2, xlim=common_xlim,
                xlabel=r"$\omega_{\parallel}$ $[$cm$^{-1}]$",
                ylabel="IR frequency [cm$^{-1}$]" if x0==0 else None)
        sim_ax.text(1100, 2500, "cavity photon", color=photon_color, fontsize=12)
        sim_ax.text(850, 2370, "C=O asym. stretch", color=molecule_color, fontsize=12)
        sim_ax.tick_params(color=annotation_color, labelsize='medium', width=2)
        sim_ax.set_xlim(common_xlim)
        sim_ax.set_ylim(common_ylim)
        sim_ax.set_xticks(xticks)
        sim_ax.set_yticks(yticks)
        sim_ax.text(0.99, 0.12, amp_list[x0], transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        sim_ax.text(0.09, 0.97, label_list[x0][0], transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        sim_ax.annotate('', xy=(450, 2450), xytext=(450, 2550), arrowprops=dict(facecolor=annotation_color, edgecolor=annotation_color, arrowstyle='->', alpha=0.8, lw=2), fontsize=20)
        sim_ax.text(0.45, 0.52, "UP excitation", transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color=annotation_color)

        phe = np.load(f"./data/mxl_144_qc/photonic_energy_36_{amp[x0]}.npy")[:3001, :]
        t = np.arange(phe.shape[0]) * 2 / 1000
        phes = [smooth(np.sum(phe[:,:5], axis=1)), smooth(np.sum(phe[:,5:10], axis=1)), smooth(np.sum(phe[:,25:30], axis=1)), smooth(np.sum(phe[:,30:35], axis=1)), smooth(phe[:,35])]
        labels = [r"$0$-$62.5$", r"$62.5$-$125$", r"$312.5$-$375$", r"$375$-$450$", r"$450$"]
        clp.plotone([t]*5, phes, energy_ax, 
                    showlegend=True if x0==0 else False, legendloc=(0.5,0.4), legendFontSize=9, labels=labels,
                    colorMap=plt.cm.hot, colorMap_endpoint=0.6, lw=1.2, xlim=[0,6],
                    xlabel="time [ps]",
                    ylabel="photon energy [a.u.]" if x0==0 else None)
        if x0 == 0:
            legend = energy_ax.get_legend()
            legend.set_title(r"$k_{\parallel}$ [cm$^{-1}$]")
        energy_ax.text(0.09, 0.97, label_list[x0][1], transform=energy_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    axes[1,0].set_yticks([0, 0.25, 0.5])
    axes[1,1].set_yticks([0, 5, 10, 15])
    axes[1,1].text(0.65, 0.85, "delayed excitation \n"+r"of $k_{\parallel}\approx 0$", transform=axes[1,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', color="r")
    axes[1,1].annotate('', xy=(3, 7), xytext=(4, 10), arrowprops=dict(facecolor=annotation_color, edgecolor="#ff4d6d", arrowstyle='->', alpha=0.8, lw=2), fontsize=20)

    model_ax = axes[2,0]
    model_data = np.load("./data/mxl_144_qc/model_1d_inverted_gaussian_spectrum.npz", allow_pickle=True)
    model_spectrum = model_data["spectrum"]
    model_extent = model_data["extent"]
    omega_m = float(model_data["omega_m"])
    omega_perp = float(model_data["omega_perp"])
    model_spectrum, model_extent = crop_spectrum_to_window_origin_lower(model_spectrum, model_extent)

    vmax = np.max(np.max(model_spectrum))
    vmin = vmax * 0.0001
    model_pos = model_ax.imshow(model_spectrum, aspect='auto', extent=model_extent,
            cmap=cm.inferno,
            origin="lower",
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    cbar_ax = model_ax.inset_axes([1.04, 0.0, 0.055, 1.0])
    cbar = model_ax.figure.colorbar(model_pos, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("spectral intensity [arb. units]", fontsize=12)
    freq_cav_inplane = np.linspace(common_xlim[0], common_xlim[1], 400)
    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * omega_m, (omega_perp**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, model_ax, showlegend=False, colors=[molecule_color, photon_color], linestyles=["--", "--"], lw=1.2, xlim=common_xlim,
            xlabel=r"$\omega_{\parallel}$ $[$cm$^{-1}]$",
            ylabel="IR frequency [cm$^{-1}$]")
    model_ax.text(1100, 2500, "cavity photon", color=photon_color, fontsize=12)
    model_ax.text(850, 2370, "C=O asym. stretch", color=molecule_color, fontsize=12)
    model_ax.tick_params(color=annotation_color, labelsize='medium', width=2)
    model_ax.set_xlim(common_xlim)
    model_ax.set_ylim(common_ylim)
    model_ax.set_xticks(xticks)
    model_ax.set_yticks(yticks)
    model_ax.text(0.09, 0.97, "(e)", transform=model_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
    density_ax = model_ax.inset_axes([0.15, 0.55, 0.36, 0.28])
    density_ax.set_facecolor("none")
    density_ax.patch.set_alpha(0.0)
    density_ax.plot(model_data["x_grid"], model_data["density"], color=inset_color, lw=1.5)
    density_ax.axvline(0.2, color=molecule_color, lw=1.0, ls="--", alpha=0.9)
    density_ax.set_xlim(float(np.min(model_data["x_grid"])), float(np.max(model_data["x_grid"])))
    density_ax.set_xticks([0, 0.5, 1])
    density_ax.set_ylim(0.0, 1.1 * float(np.max(model_data["density"])))
    density_ax.set_xlabel(r"$x/L_x$", fontsize=8)
    density_ax.set_ylabel(r"$n$", fontsize=8)
    density_ax.tick_params(labelsize=7, width=1.0, direction="in", bottom=True, top=True, left=True, right=True, colors=inset_color)
    density_ax.xaxis.label.set_color(inset_color)
    density_ax.yaxis.label.set_color(inset_color)
    for spine in density_ax.spines.values():
        spine.set_color(inset_color)
    model_ax.text(0.35, 0.95, "inverted gaussian", transform=model_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='center', color=annotation_color)

    axes[2,1].axis("off")
    axes[2,1].text(0.09, 0.97, "(f)", transform=axes[2,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/fig4_neq_spectrum_with_model.png")

def plot_efield_coenergy_animation():

    from pathlib import Path
    from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter

    shape = (144, 144)
    frame_stride = 1
    fps = 25
    dpi = 160
    max_frames = None
    dt_fs = 20
    time_index_for_scale = 50
    efield_data_list = ["./data/mxl_144*144/2d_144_144_efield_neq_0.002.npy", "./data/mxl_144*144/2d_144_144_efield_neq_0.005.npy"]
    coenergy_data_list = ["./data/mxl_144*144/coenergy_36_0.002.npy", "./data/mxl_144*144/coenergy_36_0.005.npy"]
    amp_values = ["0.002", "0.005"]
    amp_list = [r"$E_0=2\times10^{-3}$ a.u.", r"$E_0=5\times10^{-3}$ a.u."]
    cbar_labels = ["E-field intensity [a.u.]", r"$E_{\rm C=O}$ per molecule [$10^2$ cm$^{-1}$]"]
    label_list = ["(a)", "(b)"]

    for amp_index, amp_value in enumerate(amp_values):
        efield_data = np.load(efield_data_list[amp_index])
        coenergy_data = np.load(coenergy_data_list[amp_index]) * (219474.63 / 36) * 0.01
        ntimes = min(efield_data.shape[0], coenergy_data.shape[0])
        frame_ids = np.arange(0, ntimes, frame_stride)
        if max_frames is not None:
            frame_ids = frame_ids[:max_frames]
        t_ps = np.arange(ntimes) * dt_fs / 1000

        efield_ref = (np.abs(efield_data[time_index_for_scale, :]).reshape(shape))**2
        efield_vmax = np.max(np.max(efield_ref))
        efield_vmin = efield_vmax * 0.005
        coenergy_vmin = np.percentile(coenergy_data[time_index_for_scale], 98.5)
        coenergy_vmax = np.percentile(coenergy_data[time_index_for_scale], 99.9)
        norms = [
            LogNorm(vmin=efield_vmin, vmax=efield_vmax),
            LogNorm(vmin=coenergy_vmin, vmax=coenergy_vmax),
        ]

        fig, axes = clp.initialize(1, 2, width=10.6, height=4.8, LaTeX=True, fontsize=12, return_fig_args=True)
        images = []
        time_texts = []
        markers = []
        frame0_list = [
            (np.abs(efield_data[0, :]).reshape(shape))**2,
            coenergy_data[0],
        ]

        for col, frame0 in enumerate(frame0_list):
            pos = axes[col].imshow(frame0, aspect='equal', extent=[0, 144, 0, 144], origin="lower",
                    cmap=cm.hot, interpolation='nearest', norm=norms[col])
            images.append(pos)
            axes[col].set_box_aspect(1)
            clp.plotone([], [], axes[col], colors=["c-","c-"], ylabel=r"$L_y$ \ position \ [$\mu\rm{m}$]" if col==0 else None, xlabel=r"$L_x$ \ position \ [$\mu\rm{m}$]", showlegend=False)
            axes[col].set_yticks([0,72, 144])
            axes[col].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[col].set_xticks([0,72, 144])
            axes[col].set_xticklabels([r"$0$", r"$200$", r"$400$"])
            if col != 0:
                axes[col].set_yticklabels([])
            time_texts.append(axes[col].text(0.98, 0.18, "", transform=axes[col].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w"))
            axes[col].text(0.02, 0.98, label_list[col], transform=axes[col].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[col].text(0.98, 0.08, amp_list[amp_index], transform=axes[col].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            if col == 0:
                axes[col].text(0.55, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[col].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
                axes[col].text(0.55, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=axes[col].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")

            marker_center = (85, 72)
            marker_width = 20
            marker_height = 24
            marker_linewidth = 1.5
            center_x, center_y = marker_center
            center_marker = Rectangle(
                (center_x - marker_width / 2, center_y - marker_height / 2),
                marker_width,
                marker_height,
                fill=False,
                edgecolor="cyan",
                linewidth=marker_linewidth,
                zorder=3,
            )
            axes[col].add_patch(center_marker)
            markers.append(center_marker)

            cbar_ax = axes[col].inset_axes([1.04, 0.0, 0.055, 1.0])
            cbar = fig.colorbar(pos, cax=cbar_ax)
            cbar.ax.tick_params(labelsize=12)
            cbar.set_label(cbar_labels[col], fontsize=12)
            if col == 1:
                if amp_index == 0:
                    cbar.set_ticks([4, 5, 6, 7])
                    cbar.set_ticklabels([r"$4$", r"$5$", r"$6$", r"$7$"])
                else:
                    cbar.set_ticks([8, 9, 12, 16])
                    cbar.set_ticklabels([r"$8$", r"$9$", r"$12$", r"$15$"])
            cbar.ax.yaxis.set_label_coords(4.2, 0.5)

        fig.subplots_adjust(left=0.08, right=0.88, bottom=0.15, top=0.94, wspace=0.48)

        def update(frame_id):
            images[0].set_data((np.abs(efield_data[frame_id, :]).reshape(shape))**2)
            images[1].set_data(coenergy_data[frame_id])
            for time_text in time_texts:
                time_text.set_text(r"$t=%.1f$ ps" % t_ps[frame_id])
            return images + time_texts + markers

        ani = FuncAnimation(fig, update, frames=frame_ids, interval=1000 / fps, blit=False)
        output_file = Path(f"figures/efield_coenergy_2d_144_144_{amp_value}.mp4")
        output_file.parent.mkdir(parents=True, exist_ok=True)
        if output_file.suffix.lower() == ".gif":
            writer = PillowWriter(fps=fps)
        else:
            writer = FFMpegWriter(fps=fps, bitrate=2400)
        ani.save(output_file, writer=writer, dpi=dpi)
        plt.close(fig)
        print(f"Saved animation to {output_file}")

def plot_vg():
    fig, axes = clp.initialize(1, 2, width=4.3*2, height=4.3*0.618*1.1, LaTeX=True, fontsize=12, return_fig_args=True)
    qc_file = np.load("./data/mxl_144*144/qc_y.npy", allow_pickle=True).item()
    sp_full = qc_file['sp']
    x = qc_file['freq']
    dx = x[2] - x[1]
    nstart, nmid, nend = int(2200 / dx), int(2350 / dx), int(3000 / dx)
    x = x[nstart:nend]
    sp = np.abs(sp_full[nstart:nend,:]) / 1e32
    sp = sp[::-1, :]
    freq_cav_inplane_min = 12.5
    freq_cav_inplane_max = 12.5 * 144
    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]

    vmax = np.max(np.max(sp))
    vmin = vmax * 0.001
    pos = axes[0].imshow(sp, aspect='auto', extent=extent,
            cmap=cm.inferno,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    cbar_ax = axes[0].inset_axes([1.04, 0.0, 0.055, 1.0])
    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("spectral intensity [arb. units]", fontsize=12)

    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, axes[0], showlegend=False, colors=["orange", "c"], linestyles=["--", "--"], lw=1.2, xlim=[12.5,1800],
            xlabel=r"$\omega_{\parallel}$ [$\rm{cm}^{-1}$]",
            ylabel=r"IR frequency [$\rm{cm}^{-1}$]")
    axes[0].text(1190, 2550, "cavity photon", color='c', fontsize=12)
    axes[0].text(940, 2370, "C=O asym. stretch", color='orange', fontsize=12)
    axes[0].tick_params(color='c', labelsize='medium', width=2)

    axes[0].set_xticks([12.5,600,1200,1800])
    axes[0].set_yticks([2200,2400,2600,2800])
    axes[0].annotate('', xy=(450, 2450), xytext=(450, 2550), arrowprops=dict(facecolor='w', edgecolor='w', arrowstyle='->', alpha=0.8, lw=2), fontsize=20)

    lp = np.array([x[np.argmax(sp_full[nstart:nmid, i])] for i in range(144)])
    up = np.array([x[nmid - nstart + np.argmax(sp_full[nmid:nend, i])] for i in range(144)])

    omega_0 = 2320
    lp, up = lp.reshape(-1), up.reshape(-1)
    wp = np.linspace(12.5, 12.5*144, 144)
    wk = (omega_0**2 + wp**2)**0.5
    vg_up = (0.5 + (wk - omega_0) / (4 * up - 2 * wk - 2 * omega_0)) * (wp / wk)
    vg_lp = (0.5 + (wk - omega_0) / (4 * lp - 2 * wk - 2 * omega_0)) * (wp / wk)
    wc_p = (up - omega_0) / (2 * up - omega_0 - wk)
    xs = [wp]*2
    ys = [vg_lp, vg_up]
    clp.plotone(xs, ys, axes[1], showlegend=True, legendloc=(0.7,0.3), legendFontSize=9,
                labels=["LP","UP"], colors=["b-", "r-"], lw=1.2, xlim=[12.5,1800], ylim=[-0.01,0.6],
                xlabel=r"$\omega_{\parallel}$ $[$cm$^{-1}]$",
                ylabel=r"$v_g$ $[c]$")
    plt.rcParams["axes.axisbelow"] = False
    axes[0].text(0.08, 0.97, "(a)", transform=axes[0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
    axes[1].text(0.08, 0.97, "(b)", transform=axes[1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    axes[1].set_xticks([12.5,600,1200,1800])
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS1_vg.png")

def plot_1d_efield_else():

    axes = clp.initialize(3, 3, width=3*4.3, height=4.3*0.618*1.1*3, LaTeX=True, fontsize=12)
    uplist = [27, 30, 48]
    k_parallel_list = [337.5, 375, 600]
    coeff = np.array([0.08123523, 0.09236882, 0.17112003])
    wc_up = np.array([0.56, 0.58, 0.68])
    length_list = [2600]*3
    vg = np.array([r"$%.3f\,c$"%(coeff[i]) for i in range(3)])
    label_list = {0 : ["(a)","(b)", "(c)"], 1: ["(d)","(e)","(f)"], 2: ["(g)","(h)","(i)"]}
    amp_list = [r"$%d\times10^{-3}$" % i for i in range(1,9,2)]
    for j in range(3):
        sp_list = [np.load(f"./data/mxl_144_efield/effective_efield_{uplist[j]}_{i}.npy")[50:2551].T for i in ["0.001", "0.003", "0.005", "0.007"]]
        mt_list, mmsd_list = get_adaptive_emsd(uplist[j], coeff[j], amp_list=["0.001", "0.003", "0.005", "0.007"], sp_list=sp_list)
        clp.plotone(mt_list, mmsd_list, axes[2,j], ylabel=r"$\mathrm{MSD}$ \ [$10^{4} \ \mu\rm{m}^2$]" if j == 0 else None, xlabel="time [ps]",
                    showlegend=True if j==0 else False, legendloc=(0.25,0.3), legendFontSize=8, 
                    xlim=[0,mt_list[0][-1]], ylim=[0,0.05*(j+1)],
                    labels=amp_list, colorMap=plt.cm.hot, colorMap_endpoint=0.6, alpha=0.4)
        quad_t = mt_list[0][(mt_list[0] >= 1.0) & (mt_list[0] <= 4.5)]
        linear_t = mt_list[0][(mt_list[0] >= 2.0) & (mt_list[0] <= 5.0)]
        quad_fit_params = [[0.004531527975714095, -0.003958012370730053, 0.0010863142741930323],
                           [0.005996671732206156, -0.005798312959837168, 0.002028280858851599],
                           [0.021316608498041262, -0.017371908799149804, 0.004642095526640439]]
        linear_fit_params = [[0.024528881948415075, -0.04702675946500917],
                             [0.03571352305202849, -0.07123119558557142],
                             [0.13980444375345502, -0.29105386505480835]]
        quad_sampled_t = quad_fit_params[j][0] * quad_t**2 + quad_fit_params[j][1] * quad_t + quad_fit_params[j][2]
        linear_sampled_t = linear_fit_params[j][0] * linear_t + linear_fit_params[j][1]
        clp.plotone([quad_t], [quad_sampled_t], axes[2,j], colors=["g--"], showlegend=False, lw=1.5)
        clp.plotone([linear_t], [linear_sampled_t], axes[2,j], colors=["b--"], showlegend=False, lw=1.5)
        axes[2,j].text(0.08, 0.97, label_list[2][j], transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[2,j].text(0.55, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[j]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,j].text(0.55, 0.85, r"$W_{\rm ph}=$"+rf"$ \ {wc_up[j]} \ $", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,j].axvline(x=1, linestyle='-.', alpha=0.3, color="c")
        axes[2,j].text(0.1, 0.85, "pulse \n on", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', alpha=0.3, color="c")
        axes[2,j].text(0.25, 0.18, r"$\propto t^2$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="g")
        axes[2,j].text(0.85, 0.4 if j == 2 else 0.3, r"$\propto t$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="b")
    
    axes[2,0].get_legend().set_title(r"$E_0$ [a.u.]")
    axes[2,0].set_yticks([0, 0.05, 0.1])
    axes[2,0].set_yticklabels([r"$0$", r"$0.8$", r"$1.6$"])
    axes[2,1].set_yticks([0, 0.05, 0.1, 0.15])
    axes[2,1].set_yticklabels([r"$0$", r"$0.8$", r"$1.6$", r"$2.4$"])
    axes[2,2].set_yticks([0, 0.2, 0.4])
    axes[2,2].set_yticklabels([r"$0$", r"$3.2$", r"$6.4$"])

    amp = [0.001, 0.007]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    for y0 in range(2):
        ref = np.load(f"./data/mxl_144_efield/effective_efield_24_{amp[y0]}.npy")
        vmax = np.max(np.max(ref))
        vmin = vmax * 0.01 if y0 == 0 else vmax * 0.005
        for x0 in range(3):
            sp_avg = np.load(f"./data/mxl_144_efield/effective_efield_{uplist[x0]}_{amp[y0]}.npy")
            ntimes = sp_avg.shape[0]
            sp_avg = np.abs(sp_avg).T
            extent = [0 , 20, 0, 144]
            pos = axes[y0, x0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            t = np.linspace(0,(ntimes-1)*2,ntimes)[:length_list[x0]]/1000
            xs = [t]
            ys = [144*0.2+144*0.75*coeff[x0]*t]
            clp.plotone(xs, ys, axes[y0, x0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if x0==0 else None, xlabel="time [ps]" if y0==1 else None, showlegend=False, lw=1.5)
            axes[y0, x0].set_yticks([0,72, 144])
            axes[y0, x0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[y0, x0].text(0.55, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[x0]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.55, 0.85, r"$W_{\rm ph}=$"+rf"$ \ {wc_up[x0]} \ $", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.98, 0.1, amp_list[y0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[y0, x0].text(0.08, 0.98, label_list[y0][x0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[y0, x0].text(0.02, 0.1, r"$\tilde{v}_g \ = \ $"+f"{vg[x0]}", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="c")
            if x0 == 2 : 
                cbar_ax = axes[y0, x0].inset_axes([1.04, 0.0, 0.055, 1.0])
                cbar = axes[y0, x0].figure.colorbar(pos, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=12)
                cbar.set_label("E-field intensity [a.u.]", fontsize=12)
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS2_1d_efield_else.png")

def plot_1d_efield_ndependence():

    axes = clp.initialize(2, 4, width=4.3*4, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12)
    label_list = {0 : ["(a)","(b)","(c)","(d)"], 1: ["(e)","(f)","(g)","(h)"]}
    amp_list = [r"$7\times10^{-3}$"]
    n_list = [144, 288, 576, 1152]

    for j in range(4):
        sp_list = [np.load(f"./data/mxl_n_dependence/effective_efield_{n_list[j]}_36_0.007.npy")[50:2551].T]
        mt_list, mmsd_list = get_adaptive_emsd(36, 0.11638023, amp_list=["0.007"], folder="./data/mxl_n_dependence/", sp_list=sp_list)
        clp.plotone(mt_list, mmsd_list, axes[1,j], ylabel=r"$\mathrm{MSD}$ \ [$10^{4} \ \mu\rm{m}^2$]" if j == 0 else None, xlabel="time [ps]",
                    showlegend=True if j==0 else False, legendloc=(0.72,0.7), legendFontSize=8, 
                    xlim=[0,mt_list[0][-1]], ylim=[0,0.05],
                    labels=amp_list, colorMap=plt.cm.hot, colorMap_endpoint=0.6, alpha=0.4)
        axes[1,j].text(0.08, 0.98, label_list[1][j], transform=axes[1,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[1,j].text(0.6, 0.35, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[1,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[1,j].text(0.6, 0.25, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=axes[1,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[1,j].text(0.6, 0.15, r"$N_{\rm grid}= \ %d$" % n_list[j], transform=axes[1,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[1,j].set_yticks([0, 0.025, 0.05])
        axes[1,j].set_yticklabels([r"$0$", r"$0.4$", r"$0.8$"])
        axes[1,j].axvline(x=1, linestyle='-.', alpha=0.3, color="c")
        axes[1,j].text(0.1, 0.85, "pulse \n on", transform=axes[1,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', alpha=0.3, color="c")
    axes[1,0].get_legend().set_title(r"$E_0$ [a.u.]")

    refsp = np.load(f"./data/mxl_n_dependence/effective_efield_1152_36_0.007.npy")
    vmax = np.max(np.max(refsp))
    vmin = vmax * 0.01
    for x0 in range(4):
        sp_avg = np.load(f"./data/mxl_n_dependence/effective_efield_{n_list[x0]}_36_0.007.npy")
        ntimes = sp_avg.shape[0]
        sp_avg = np.abs(sp_avg).T
        extent = [0 , 20, 0, 144]
        pos = axes[0, x0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                cmap=cm.hot, interpolation='nearest',
                norm=LogNorm(vmin=vmin, vmax=vmax))
        t = np.linspace(0,(ntimes-1)*2,ntimes)[:2600]/1000
        xs = [t]
        ys = [144*0.2+144*0.75*0.11638023*t]
        clp.plotone(xs, ys, axes[0, x0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if x0==0 else None, xlabel="time [ps]", showlegend=False, lw=1.5)
        axes[0, x0].set_yticks([0,72, 144])
        axes[0, x0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
        axes[0, x0].text(0.6, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        axes[0, x0].text(0.6, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=axes[0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        axes[0, x0].text(0.6, 0.75, r"$N_{\rm grid}= \ %d$" % n_list[x0], transform=axes[0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        axes[0, x0].text(0.98, 0.1, r"$E_0=7\times10^{-3}$ a.u.", transform=axes[0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        axes[0, x0].text(0.08, 0.98, label_list[0][x0], transform=axes[0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        axes[0, x0].text(0.02, 0.1, r"$\tilde{v}_g \ = \ 0.116\,c$", transform=axes[0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="c")
        if x0 == 3: 
            cbar_ax = axes[0, x0].inset_axes([1.04, 0.0, 0.055, 1.0])
            cbar = axes[0, x0].figure.colorbar(pos, cax=cbar_ax)
            cbar.ax.tick_params(labelsize=12)
            cbar.set_label("E-field intensity [a.u.]", fontsize=12)
            cbar.ax.minorticks_off()
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS3_1d_efield_numberdependence.png")

def plot_addition_mmsd():

    axes = clp.initialize(2, 4, width=4*4.3, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12, sharex=True, sharey=True)
    up_list = np.array([24, 27, 30, 48])
    k_parallel_list = [300, 337.5, 375, 600]
    amp = ["0.001", "0.007"]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    label_list = {0 : ["(a)","(b)","(c)","(d)"], 1 : ["(e)","(f)","(g)","(h)"]}

    for y0 in range(2):
        for x0 in range(4):
            sp_avg = np.load(f"./data/mxl_144_coenergy/coenergy_{up_list[x0]}_{amp[y0]}.npy")
            ntimes = sp_avg.shape[0]
            sp_avg = np.abs(sp_avg).T * (219474.63 / 36) * 10**(-(y0+2)) # convert to cm^-1 and scale for better visualization
            vmin = np.percentile(sp_avg, 60)
            vmax = np.percentile(sp_avg, 99)
            extent = [0 , 20, 0, 144]
            pos = axes[y0,x0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            clp.plotone([], [], axes[y0,x0], colors=["c-","c-"], ylabel=r"$L_x$ position [$\mu$m]" if x0==0 else None, xlabel="time [ps]" if y0==1 else None, showlegend=False)
            axes[y0,x0].set_yticks([0,72, 144])
            axes[y0,x0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[y0,x0].text(0.1, 0.97, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[x0]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[y0,x0].transAxes, fontsize=13, fontweight='bold', va='top', ha='left', color="#00e5ff")
            axes[y0,x0].text(0.98, 0.97, amp_list[y0], transform=axes[y0,x0].transAxes, fontsize=13, fontweight='bold', va='top', ha='right', color="#00e5ff")
            axes[y0,x0].text(0.08, 0.97, label_list[y0][x0], transform=axes[y0,x0].transAxes, fontsize=13, fontweight='bold', va='top', ha='right', color="#00e5ff")
            cbar_ax = axes[y0,x0].inset_axes([1.04, 0.0, 0.055, 1.0])
            cbar = axes[y0,x0].figure.colorbar(pos, cax=cbar_ax)
            if y0 == 0:
                cbar.set_ticks([2.2, 2.5, 2.8, 3])
                cbar.set_ticklabels([r"$2.2$", r"$2.5$", r"$2.8$", r"$3.0$"])
                cbar.ax.minorticks_off()
            else:
                cbar.set_ticks([1, 2, 3, 4])
                cbar.set_ticklabels([r"$1.0$", r"$2.0$", r"$3.0$", r"$4.0$"])
                cbar.ax.minorticks_off()
            cbar.ax.tick_params(labelsize=12)
            if x0 == 3: cbar.set_label(r"C=O energy [$10^{%d} \ $cm$^{-1}$]"%(y0+2), fontsize=12)
            axes[y0,x0].hlines(144*0.05, 0, 20, colors="c", linestyles="--", lw=1.5)
            axes[y0,x0].hlines(144*0.35, 0, 20, colors="c", linestyles="--", lw=1.5)

    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS4_addition_coenergy.png")

def plot_convergence():

    axes = clp.initialize(3, 4, width=4*4.3, height=4.3*0.618*1.1*3, LaTeX=True, fontsize=12)
    label_list = {0 : ["(a)","(b)","(c)","(d)"], 1: ["(e)","(f)","(g)","(h)"], 2: ["(i)","(j)","(k)","(l)"]}
    repeat_list = [1, 3, 5, 10]
    amp = [0.001, 0.007]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    for y0 in range(2):
        eqref = np.load(f"./data/mxl_convergence_test/multimode_cavmd_144_24_{amp[y0]}_1_neq_final.npy")
        vmax = np.max(np.max(eqref))
        vmin = vmax * 0.01 if y0 == 0 else vmax * 0.005
        for x0 in range(4):
            sp_avg = np.load(f"./data/mxl_convergence_test/multimode_cavmd_144_24_{amp[y0]}_{repeat_list[x0]}_neq_final.npy")
            ntimes = sp_avg.shape[0]
            sp_avg = np.abs(sp_avg).T
            extent = [0 , 20, 0, 144]
            pos = axes[y0, x0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            t = np.linspace(0,(ntimes-1)*2,ntimes)[:2600]/1000
            xs = [t]
            ys = [144*0.2+144*0.75*0.07069759*t]
            clp.plotone(xs, ys, axes[y0, x0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if x0==0 else None, xlabel="time [ps]" if y0==1 else None, showlegend=False, lw=1.5)
            axes[y0, x0].set_yticks([0,72, 144])
            axes[y0, x0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[y0, x0].text(0.6, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 300 \ $"+r"$\rm{cm}^{-1}$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.6, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.55 \ $", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.98, 0.1, amp_list[y0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[y0, x0].text(0.02, 0.98, label_list[y0][x0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.02, 0.1, r"$\tilde{v}_g \ = \ $"+r"$0.071 \ c$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="c")
            axes[y0, x0].text(0.98, 0.2, f"averaged over {repeat_list[x0]} times", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            cbar_ax = axes[y0, x0].inset_axes([1.04, 0.0, 0.055, 1.0])
            cbar = axes[y0, x0].figure.colorbar(pos, cax=cbar_ax)
            cbar.ax.tick_params(labelsize=12)
            cbar.set_label("E-field intensity [a.u.]", fontsize=12)

    amp_name_list = ["0.001", "0.0016", "0.003", "0.005", "0.007"]
    amp_list = [r"$%.1f\times10^{-3}$" % (float(i)*1000) for i in amp_name_list]
    for x0 in range(4):
        sp_list = []
        for y0 in range(len(amp_name_list)):
            sp_avg = np.load(f"./data/mxl_convergence_test/multimode_cavmd_144_24_{amp_name_list[y0]}_{repeat_list[x0]}_neq_final.npy")
            sp_list.append(sp_avg[50:2551, :].T)
        
        mt_list, mmsd_list = get_adaptive_emsd(24, 0.07069759, 
                                               amp_list=amp_name_list, 
                                               folder="./data/mxl_convergence_test/",
                                               sp_list=sp_list)
        clp.plotone(mt_list, mmsd_list, axes[2,x0], ylabel=r"$\mathrm{MSD}$ \ [$10^{4} \ \mu\rm{m}^2$]" if x0 == 0 else None, xlabel="time [ps]",
                    showlegend=True if x0 == 0 else False, legendloc=(0.05,0.3), legendFontSize=8, 
                    xlim=[0,mt_list[0][-1]], ylim=[0, 0.05],
                    labels=amp_list, colorMap=plt.cm.hot, colorMap_endpoint=0.6, alpha=0.5)

        axes[2,x0].text(0.02, 0.98, label_list[2][x0], transform=axes[2,x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,x0].text(0.68, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 300 \ $"+r"$\rm{cm}^{-1}$", transform=axes[2,x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[2,x0].text(0.32, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.55 \ $", transform=axes[2,x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,x0].axvline(x=1, linestyle='-.', alpha=0.3, color="c")
        axes[2,x0].text(0.1, 0.9, "pulse \n on", transform=axes[2,x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', alpha=0.3, color="c")
        axes[2,x0].text(0.22, 0.18, r"$\propto t^2$", transform=axes[2,x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="g")
        axes[2,x0].text(0.85, 0.65, r"$\propto t$", transform=axes[2,x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="b")
        axes[2,x0].set_yticks([0, 0.025, 0.05])
        axes[2,x0].set_yticklabels([r"$0$", r"$0.4$", r"$0.8$"])
    
    axes[2,0].get_legend().set_title(r"$E_0$ [a.u.]")
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS5_convergence.png")

def plot_harmonic():
    axes = clp.initialize(4, 2, width=4.3*2, height=4.3*0.618*1.1*4, LaTeX=True, fontsize=12)
    molecule_color = "#ffb000"
    photon_color = "#00d5ff"
    annotation_color = "w"
    inset_color = "w"
    label_list = {0 : ["(a)","(c)","(e)","(g)"], 1 : ["(b)","(d)","(f)","(h)"]}
    amp_list = [r"$E_0=1\times 10^{-3}$ a.u.", r"$E_0=7\times 10^{-3}$ a.u."]
    amp = [0.001, 0.007]
    N = 144
    domega = 5
    common_xlim = [5, 720]
    common_ylim = [400, 1100]
    xticks = [5, 240, 480, 720]
    yticks = [400, 600, 800, 1000]

    def crop_spectrum_to_window(spectrum, extent):
        x_axis = np.linspace(extent[0], extent[1], spectrum.shape[1])
        y_axis = np.linspace(extent[3], extent[2], spectrum.shape[0])
        x_keep = (x_axis >= common_xlim[0]) & (x_axis <= common_xlim[1])
        y_keep = (y_axis >= common_ylim[0]) & (y_axis <= common_ylim[1])
        cropped = spectrum[np.ix_(y_keep, x_keep)]
        cropped_extent = [x_axis[x_keep][0], x_axis[x_keep][-1], y_axis[y_keep][-1], y_axis[y_keep][0]]
        return cropped, cropped_extent

    for x0 in range(2):
        sim_ax = axes[0,x0]
        energy_ax = axes[1,x0]
        qc_file = np.load(f"./data/mxl_144_harmonic/qc_spectra_96_{amp[x0]}.npy", allow_pickle=True).item()
        x = qc_file["freq"]
        sp = qc_file["sp"]
        dx = x[1] - x[0]
        nstart, nmid, nend = int(400 / dx), int(650 / dx), int(1100 / dx)
        x = x[nstart:nend]
        sp = np.abs(sp[nstart:nend,:]) / 1e34
        sp = sp[::-1, :]
        freq_cav_inplane_min = domega
        freq_cav_inplane_max = domega * N
        extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]
        sp, extent = crop_spectrum_to_window(sp, extent)

        vmax = np.max(np.max(sp))
        vmin = vmax * 0.0001
        pos = sim_ax.imshow(sp, aspect='auto', extent=extent,
                cmap=cm.inferno,
                interpolation='nearest',
                norm=LogNorm(vmin=vmin, vmax=vmax)
                )
        cbar_ax = sim_ax.inset_axes([1.04, 0.0, 0.055, 1.0])
        cbar = sim_ax.figure.colorbar(pos, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label("spectral intensity [arb. units]", fontsize=12)

        freq_cav_inplane = np.linspace(common_xlim[0], common_xlim[1], 400)

        xs = [freq_cav_inplane]*2
        ys = [np.ones(len(freq_cav_inplane)) * 655, (655**2 + freq_cav_inplane**2)**0.5]
        clp.plotone(xs, ys, sim_ax, showlegend=False, colors=[molecule_color, photon_color], linestyles=["--", "--"], lw=1.2, xlim=common_xlim,
                xlabel=r"$\omega_{\parallel}$ $[$cm$^{-1}]$",
                ylabel="IR frequency [cm$^{-1}$]" if x0==0 else None)
        sim_ax.text(450, 750, "cavity photon", color=photon_color, fontsize=12)
        sim_ax.text(400, 670, "O=C=O bending", color=molecule_color, fontsize=12)
        sim_ax.tick_params(color=annotation_color, labelsize='medium', width=2)
        sim_ax.set_xlim(common_xlim)
        sim_ax.set_ylim(common_ylim)
        sim_ax.set_xticks(xticks)
        sim_ax.set_yticks(yticks)
        sim_ax.text(0.99, 0.12, amp_list[x0], transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        sim_ax.text(0.09, 0.97, label_list[x0][0], transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
        sim_ax.annotate('', xy=(480, 950), xytext=(480, 1050), arrowprops=dict(facecolor=annotation_color, edgecolor=annotation_color, arrowstyle='->', alpha=0.8, lw=2), fontsize=20)
        sim_ax.text(0.55, 0.9, "UP excitation", transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color=annotation_color)

        phe = np.load(f"./data/mxl_144_harmonic/photonic_energy_96_{amp[x0]}.npy")[:3001, :]
        t = np.arange(phe.shape[0]) * 2 / 1000
        phes = [smooth(np.sum(phe[:,:24], axis=1)), smooth(np.sum(phe[:,24:48], axis=1)), smooth(np.sum(phe[:,48:72], axis=1)), smooth(np.sum(phe[:,72:96], axis=1)), smooth(phe[:,96])]
        labels = [r"$0$-$115$", r"$120$-$235$", r"$240$-$355$", r"$360$-$475$", r"$480$"]
        clp.plotone([t]*5, phes, energy_ax, 
                    showlegend=True if x0==1 else False, legendloc=(0.5,0.4), legendFontSize=9, labels=labels,
                    colorMap=plt.cm.hot, colorMap_endpoint=0.6, lw=1.2, xlim=[0,6],
                    xlabel="time [ps]",
                    ylabel="photon energy [a.u.]" if x0==0 else None)
        if x0 == 1:
            legend = energy_ax.get_legend()
            legend.set_title(r"$k_{\parallel}$ [cm$^{-1}$]")
        energy_ax.text(0.09, 0.97, label_list[x0][1], transform=energy_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    axes[1,0].set_yticks([0, 0.25, 0.5])
    axes[1,1].set_yticks([0, 5, 10, 15])

    amp = [0.001, 0.007]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    for y0 in range(2):
        sp_avg = np.load(f"./data/mxl_144_harmonic/effective_efield_96_{amp[y0]}.npy")
        vmax = np.max(np.max(sp_avg))
        vmin = vmax * 0.08
        ntimes = sp_avg.shape[0]
        sp_avg = np.abs(sp_avg).T
        extent = [0 , 20, 0, 144]
        pos = axes[2,y0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                cmap=cm.hot, interpolation='nearest',
                norm=LogNorm(vmin=vmin, vmax=vmax))
        clp.plotone([], [], axes[2,y0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if y0==0 else None, xlabel="time [ps]", showlegend=False, lw=1.5)
        axes[2,y0].set_yticks([0,72, 144])
        axes[2,y0].set_yticklabels([r"$0$", r"$500$", r"$1000$"])
        axes[2,y0].text(0.02, 0.98, label_list[y0][2], transform=axes[2,y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        cbar_ax = axes[2,y0].inset_axes([1.04, 0.0, 0.055, 1.0])
        cbar = axes[2,y0].figure.colorbar(pos, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label("E-field intensity [a.u.]", fontsize=12)

    amp = [0.001, 0.007]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    for y0 in range(2):
        sp_avg = np.load(f"./data/mxl_144_harmonic/coenergy_{amp[y0]}_neq.npy") * (219474.63 / 36) * 10**(-(y0+2))
        vmin = np.percentile(sp_avg, 60)
        vmax = np.percentile(sp_avg, 99)
        ntimes = sp_avg.shape[0]
        sp_avg = np.abs(sp_avg).T
        extent = [0 , 20, 0, 144]
        pos = axes[3,y0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                cmap=cm.hot, interpolation='nearest',
                norm=LogNorm(vmin=vmin, vmax=vmax))
        clp.plotone([], [], axes[3,y0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if y0==0 else None, xlabel="time [ps]", showlegend=False, lw=1.5)
        axes[3,y0].set_yticks([0,72, 144])
        axes[3,y0].set_yticklabels([r"$0$", r"$500$", r"$1000$"])
        axes[3,y0].text(0.02, 0.98, label_list[y0][3], transform=axes[3,y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        cbar_ax = axes[3,y0].inset_axes([1.04, 0.0, 0.055, 1.0])
        cbar = axes[3,y0].figure.colorbar(pos, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(r"$E_{\mathrm{C=O}}$ per molecule [$10^{%d} \ $cm$^{-1}$]"%(y0+2), fontsize=12)
        if y0 == 0:
            cbar.set_ticks([2.2, 2.5])
            cbar.set_ticklabels([r"$2.2$", r"$2.5$"])
            cbar.ax.minorticks_off()
        else:
            cbar.set_ticks([1, 2])
            cbar.set_ticklabels([r"$1.0$", r"$2.0$"])
            cbar.ax.minorticks_off()
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS6_harmonic_example.png")

def plot_diluted():
    axes = clp.initialize(1, 3, width=4.3*3, height=4.3*0.618*1.1, LaTeX=True, fontsize=12)
    molecule_color = "#ffb000"
    photon_color = "#00d5ff"
    annotation_color = "w"
    inset_color = "w"
    label_list = ["(a)","(b)","(c)"]
    amp_list = [r"$E_0=7\times 10^{-3}$ a.u."]
    amp = [0.001, 0.007]
    N = 144
    domega = 12.5
    common_xlim = [12.5, 1800]
    common_ylim = [2200, 2980]
    xticks = [12.5, 600, 1200, 1800]
    yticks = [2200, 2400, 2600, 2800]

    def crop_spectrum_to_window(spectrum, extent):
        x_axis = np.linspace(extent[0], extent[1], spectrum.shape[1])
        y_axis = np.linspace(extent[3], extent[2], spectrum.shape[0])
        x_keep = (x_axis >= common_xlim[0]) & (x_axis <= common_xlim[1])
        y_keep = (y_axis >= common_ylim[0]) & (y_axis <= common_ylim[1])
        cropped = spectrum[np.ix_(y_keep, x_keep)]
        cropped_extent = [x_axis[x_keep][0], x_axis[x_keep][-1], y_axis[y_keep][-1], y_axis[y_keep][0]]
        return cropped, cropped_extent

    sim_ax = axes[0]
    energy_ax = axes[1]
    qc_file = np.load(f"./data/mxl_144_diluted/qc_spectra_36_5000.npy", allow_pickle=True).item()
    x = qc_file["freq"]
    sp = qc_file["sp"]
    dx = x[1] - x[0]
    nstart, nmid, nend = int(2200 / dx), int(2350 / dx), int(3000 / dx)
    x = x[nstart:nend]
    sp = np.abs(sp[nstart:nend,:]) / 1e34
    sp = sp[::-1, :]
    freq_cav_inplane_min = domega
    freq_cav_inplane_max = domega * N
    extent = [freq_cav_inplane_min, freq_cav_inplane_max, x[0] , x[-1]]
    sp, extent = crop_spectrum_to_window(sp, extent)

    vmax = np.max(np.max(sp))
    vmin = vmax * 0.0001
    pos = sim_ax.imshow(sp, aspect='auto', extent=extent,
            cmap=cm.inferno,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    cbar_ax = sim_ax.inset_axes([1.04, 0.0, 0.055, 1.0])
    cbar = sim_ax.figure.colorbar(pos, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("spectral intensity [arb. units]", fontsize=12)

    freq_cav_inplane = np.linspace(common_xlim[0], common_xlim[1], 400)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 655, (655**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, sim_ax, showlegend=False, colors=[molecule_color, photon_color], linestyles=["--", "--"], lw=1.2, xlim=common_xlim,
            xlabel=r"$\omega_{\parallel}$ $[$cm$^{-1}]$",
            ylabel="IR frequency [cm$^{-1}$]")
    sim_ax.text(1100, 2500, "cavity photon", color=photon_color, fontsize=12)
    sim_ax.text(850, 2370, "C=O asym. stretch", color=molecule_color, fontsize=12)
    sim_ax.tick_params(color=annotation_color, labelsize='medium', width=2)
    sim_ax.set_xlim(common_xlim)
    sim_ax.set_ylim(common_ylim)
    sim_ax.set_xticks(xticks)
    sim_ax.set_yticks(yticks)
    sim_ax.text(0.99, 0.12, amp_list[0], transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
    sim_ax.text(0.09, 0.97, label_list[0], transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
    sim_ax.annotate('', xy=(450, 2450), xytext=(450, 2550), arrowprops=dict(facecolor=annotation_color, edgecolor=annotation_color, arrowstyle='->', alpha=0.8, lw=2), fontsize=20)
    sim_ax.text(0.45, 0.52, "UP excitation", transform=sim_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color=annotation_color)

    phe = np.load(f"./data/mxl_144_diluted/photonic_energy_36_5000.npy")[:3001, :]
    t = np.arange(phe.shape[0]) * 2 / 1000
    phes = [smooth(np.sum(phe[:,:5], axis=1)), smooth(np.sum(phe[:,5:10], axis=1)), smooth(np.sum(phe[:,25:30], axis=1)), smooth(np.sum(phe[:,30:35], axis=1)), smooth(phe[:,35])]
    labels = [r"$0$-$62.5$", r"$62.5$-$125$", r"$312.5$-$375$", r"$375$-$450$", r"$450$"]
    clp.plotone([t]*5, phes, energy_ax, 
                showlegend=True, legendloc=(0.5,0.4), legendFontSize=9, labels=labels,
                colorMap=plt.cm.hot, colorMap_endpoint=0.6, lw=1.2, xlim=[0,6],
                xlabel="time [ps]",
                ylabel="photon energy [a.u.]")
    legend = energy_ax.get_legend()
    legend.set_title(r"$k_{\parallel}$ [cm$^{-1}$]")
    energy_ax.text(0.09, 0.97, label_list[1], transform=energy_ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
    energy_ax.set_yticks([0, 5, 10, 15])

    sp_avg = np.load(f"./data/mxl_144_diluted/effective_efield_5000.npy")
    vmax = np.max(np.max(sp_avg))
    vmin = vmax * 0.01
    ntimes = sp_avg.shape[0]
    sp_avg = np.abs(sp_avg).T
    extent = [0 , 20, 0, 144]
    pos = axes[2].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
            cmap=cm.hot, interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax))
    clp.plotone([], [], axes[2], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]", xlabel="time [ps]", showlegend=False, lw=1.5)
    axes[2].set_yticks([0,72, 144])
    axes[2].set_yticklabels([r"$0$", r"$500$", r"$1000$"])
    axes[2].text(0.02, 0.98, label_list[2], transform=axes[2].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
    cbar_ax = axes[2].inset_axes([1.04, 0.0, 0.055, 1.0])
    cbar = axes[2].figure.colorbar(pos, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("E-field intensity [a.u.]", fontsize=12)

    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS7_diluted.png")

def plot_cavity_loss_scan():

    axes = clp.initialize(2, 3, width=3*4.3, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12)
    uplist = [36]
    k_parallel_list = [450]
    coeff = np.array([0.11638023])
    wc_up = np.array([0.61])
    length_list = [2600]
    vg = np.array([r"$%.3f\,c$"%(coeff[i]) for i in range(1)])
    label_list = {0 : ["(a)","(b)", "(c)"], 1: ["(d)","(e)","(f)"]}
    amp_list = [r"$E_0=7\times10^{-3}$ a.u."]
    tau_dict = {0 : [5000, 2000, 1000], 1 : [800, 500, 100]}
    qc_dict = {0 : [0.066, 0.105, 0.151], 1 : [0.171, 0.210, 0.480]}
    ref = np.load(f"./data/mxl_144_cavity_loss_test/effective_efield_5000.npy")
    vmax = np.max(np.max(ref))
    vmin = vmax * 0.005
    for y0 in range(2):
        for x0 in range(3):
            sp_avg = np.load(f"./data/mxl_144_cavity_loss_test/effective_efield_{tau_dict[y0][x0]}.npy")
            ntimes = sp_avg.shape[0]
            sp_avg = np.abs(sp_avg).T
            extent = [0 , 20, 0, 144]
            pos = axes[y0, x0].imshow(sp_avg, aspect='auto', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            t = np.linspace(0,(ntimes-1)*2,ntimes)[:length_list[0]]/1000
            xs = [t]
            ys = [144*0.2+144*0.75*coeff[0]*t]
            clp.plotone(xs, ys, axes[y0, x0], colors=["c--"], ylabel=r"$L_x$ position \ [$\mu\rm{m}$]" if x0==0 else None, xlabel="time [ps]" if y0==1 else None, showlegend=False, lw=1.5)
            axes[y0, x0].set_yticks([0,72, 144])
            axes[y0, x0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[y0, x0].text(0.55, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[0]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.55, 0.85, r"$W_{\rm ph}=$"+rf"$ \ {wc_up[0]} \ $", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.55, 0.75, r"$\tau=$"+rf"$ \ {tau_dict[y0][x0]} \ $" + r"$\rm{fs}$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.55, 0.65, r"$Q_{\rm c}=$"+rf"$ \ {qc_dict[y0][x0]} \ $" + r"$\rm{a.u.}$", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
            axes[y0, x0].text(0.98, 0.1, amp_list[0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[y0, x0].text(0.08, 0.98, label_list[y0][x0], transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[y0, x0].text(0.02, 0.1, r"$\tilde{v}_g \ = \ $"+f"{vg[0]}", transform=axes[y0, x0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="c")
            if x0 == 2 : 
                cbar_ax = axes[y0, x0].inset_axes([1.04, 0.0, 0.055, 1.0])
                cbar = axes[y0, x0].figure.colorbar(pos, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=12)
                cbar.set_label("E-field intensity [a.u.]", fontsize=12)
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS8_cavity_loss_scan.png")



if __name__ == "__main__":
    plot_2d_final()
    plot_1d_efield()
    plot_mmsd()
    plot_neq_spectrum_with_model()
    plot_efield_coenergy_animation()
    plot_vg()
    plot_1d_efield_else()
    plot_1d_efield_ndependence()
    plot_addition_mmsd()
    plot_convergence()
    plot_harmonic()
    plot_diluted()
    plot_cavity_loss_scan()
