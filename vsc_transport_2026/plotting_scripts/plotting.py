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

def plot_2d_efield():

    fig, axes = clp.initialize(2, 4, width=4*4.3, height=4.3*2*0.86, LaTeX=True, fontsize=12, return_fig_args=True)
    data_list = ["./data/mxl_144*144/2d_144_144_efield_neq_0.002.npy", "./data/mxl_144*144/2d_144_144_efield_neq_0.005.npy"]
    time_list = [50, 150, 300]
    time_labels = [r"$t=%.1f$ ps" % (t/50) for t in time_list]
    label_list = {0 : ["(b)","(c)","(d)"], 1: ["(e)","(f)", "(g)"]}
    amp_list = [r"$E_0=2\times10^{-3}$ a.u.", r"$E_0=5\times10^{-3}$ a.u."]

    for x0 in range(2):
        data = np.load(data_list[x0])
        ref_sp = (np.abs(data[time_list[0], :]).reshape(144, 144))**2
        vmax = np.max(np.max(ref_sp))
        vmin = vmax * 0.005
        for y0 in range(1,4):
            sp = (np.abs(data[time_list[y0-1], :]).reshape(144, 144))**2
            extent = [0 , 144, 0, 144]
            pos = axes[x0, y0].imshow(sp, aspect='equal', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            axes[x0, y0].set_box_aspect(1)
            clp.plotone([], [], axes[x0, y0], colors=["c-","c-"], ylabel=r"$L_y$ \ position \ [$\mu\rm{m}$]" if y0==1 else None, xlabel=r"$L_x$ \ position \ [$\mu\rm{m}$]" if x0==1 else None, showlegend=False)
            axes[x0, y0].set_yticks([0,72, 144])
            axes[x0, y0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[x0, y0].set_xticks([0,72, 144])
            axes[x0, y0].set_xticklabels([r"$0$", r"$200$", r"$400$"])
            axes[x0, y0].text(0.98, 0.18, time_labels[y0-1], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[x0, y0].text(0.08, 0.98, label_list[x0][y0-1], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[x0, y0].text(0.98, 0.08, amp_list[x0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")

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
            axes[x0, y0].add_patch(center_marker)
            if y0 == 3:
                cbar_ax = axes[x0, y0].inset_axes([1.04, 0.0, 0.055, 1.0])
                cbar = fig.colorbar(pos, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=12)
                cbar.set_label("E-field intensity [a.u.]", fontsize=12)

    for j in range(2):
        axes[j, 1].text(0.6, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[j, 1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        axes[j, 1].text(0.6, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=axes[j, 1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")

    with h5py.File(f"./data/mxl_144*144/multimode_cavmd_144_eq.h5", "r") as f:
        data = {key: f[key][:] for key in f.keys()}
    sp = np.zeros((5000,144))
    for i in range(144):
        x, sp[:,i] = ir_spectrum(data['qc'][:,i,1], 2)
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
    pos = axes[1,0].imshow(sp, aspect='auto', extent=extent,
            cmap=cm.inferno,
            interpolation='nearest',
            norm=LogNorm(vmin=vmin, vmax=vmax)
            )
    cbar_ax = axes[1,0].inset_axes([1.04, 0.0, 0.055, 1.0])
    cbar = fig.colorbar(pos, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label("spectral intensity [arb. units]", fontsize=12)

    freq_cav_inplane = np.linspace(freq_cav_inplane_min, freq_cav_inplane_max, x.size)

    xs = [freq_cav_inplane]*2
    ys = [np.ones(len(freq_cav_inplane)) * 2327, (2320.0**2 + freq_cav_inplane**2)**0.5]
    clp.plotone(xs, ys, axes[1,0], showlegend=False, colors=["orange", "c"], linestyles=["--", "--"], lw=1.2, xlim=[12.5,1800],
            xlabel=r"$\omega_{\parallel}$ [$\rm{cm}^{-1}$]",
            ylabel=r"IR frequency [$\rm{cm}^{-1}$]")
    axes[1,0].text(1190, 2550, "cavity photon", color='c', fontsize=12)
    axes[1,0].text(940, 2370, "C=O asym. stretch", color='orange', fontsize=12)
    axes[1,0].tick_params(color='c', labelsize='medium', width=2)

    axes[1,0].set_xticks([12.5,600,1200,1800])
    axes[1,0].set_yticks([2200,2400,2600,2800])
    axes[1,0].annotate('', xy=(450, 2450), xytext=(450, 2550), arrowprops=dict(facecolor='w', edgecolor='w', arrowstyle='->', alpha=0.8, lw=2), fontsize=20)
    axes[1,0].text(0.38, 0.5, "UP excitation", transform=axes[1,0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")

    axes[0,0].axis("off")

    plt.rcParams["axes.axisbelow"] = False
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.10, top=0.94, wspace=0.18, hspace=0.20)
    fig_width, fig_height = fig.get_size_inches()
    image_group_shift = 0.055
    image_col_gap = 0.027
    for row in range(2):
        base_pos = axes[row, 1].get_position()
        square_width = base_pos.height * fig_height / fig_width
        for col in range(1, 4):
            axes[row, col].set_position([
                base_pos.x0 + image_group_shift + (col - 1) * (square_width + image_col_gap),
                base_pos.y0,
                square_width,
                base_pos.height,
            ])
    clp.adjust(savefile=f"./figures/fig1_2d_efield.png")

def plot_2d_efield_final():

    plot_2d_efield()
    from pathlib import Path
    from PIL import Image

    base_file = Path("./figures/fig1_2d_efield.png")
    inset_file = Path("./figures/fig1a.png")
    output_file = Path("./figures/fig1_2d_efield_final.png")
    paste_left = 190
    paste_top = 35
    paste_right = 1294
    paste_height = 860
    label_pos = (paste_right - 25, paste_top + 32)

    base = Image.open(base_file).convert("RGBA")
    inset = Image.open(inset_file).convert("RGBA")

    width = paste_right - paste_left
    resampling = getattr(Image, "Resampling", Image).LANCZOS
    inset = inset.resize((width, width), resampling)
    top_crop = 185
    inset = inset.crop((0, top_crop, width, top_crop + paste_height))

    final = base.copy()
    final.alpha_composite(inset, dest=(paste_left, paste_top))
    output_file.parent.mkdir(parents=True, exist_ok=True)

    dpi = 300
    fig = plt.figure(figsize=(base.width / dpi, base.height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(final)
    ax.axis("off")
    fig.text(label_pos[0] / base.width, 1 - label_pos[1] / base.height, "(a)", fontsize=12, fontweight=500, color="w", ha="right", va="top")
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
    amp_list = [r"$%d\times10^{-3}$" % i for i in range(1,9,2)]
    for j in range(2):
        mt_list, mmsd_list = get_adaptive_emsd(uplist[j], coeff[j], amp_list=["0.001", "0.003", "0.005", "0.007"])
        clp.plotone(mt_list, mmsd_list, axes[2,j], ylabel=r"$\mathrm{MSD}$ \ [$10^{4} \ \mu\rm{m}^2$]" if j == 0 else None, xlabel="time [ps]",
                    showlegend=True if j==0 else False, legendloc=(0.22,0.35), legendFontSize=8, 
                    xlim=[0,mt_list[0][-1]], ylim=[0, 0.05] if j==1 else [0, 0.15],
                    labels=amp_list, colorMap=plt.cm.hot, colorMap_endpoint=0.6, alpha=0.4)
        quad_t = mt_list[0][(mt_list[0] >= 1.0) & (mt_list[0] <= 4.0)]
        linear_t = mt_list[0][(mt_list[0] >= 2.0) & (mt_list[0] <= 5.0)]
        quad_fit_params = [[0.008, -0.008608133962878038, 0.003689480540663409],
                           [0.003316594835064373, -0.001936205081767636, 0.00029550969845358165]]
        linear_fit_params = [[0.044870112820745794, -0.08911820204841633],
                             [0.013325053829792976, -0.019916468480877206]]
        quad_sampled_t = quad_fit_params[j][0] * quad_t**2 + quad_fit_params[j][1] * quad_t + quad_fit_params[j][2]
        linear_sampled_t = linear_fit_params[j][0] * linear_t + linear_fit_params[j][1]
        clp.plotone([quad_t], [quad_sampled_t], axes[2,j], colors=["g--"], showlegend=False, lw=1.5)
        clp.plotone([linear_t], [linear_sampled_t], axes[2,j], colors=["b--"], showlegend=False, lw=1.5)
        axes[2,j].text(0.02, 0.98, label_list[2][j], transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,j].text(0.68, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ {k_parallel_list[j]} \ $"+r"$\rm{cm}^{-1}$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[2,j].text(0.32, 0.85, r"$W_{\rm ph}=$"+rf"$ \ {wc_up[j]} \ $", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[2,j].axvline(x=1, linestyle='-.', alpha=0.3, color="c")
        axes[2,j].text(0.1, 0.85, "pulse \n on", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='center', alpha=0.3, color="c")
        axes[2,j].text(0.22, 0.18, r"$\propto t^2$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="g")
        axes[2,j].text(0.85, 0.65, r"$\propto t$", transform=axes[2,j].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="b")
    
    axes[2,1].set_yticks([0, 0.025, 0.05])
    axes[2,1].set_yticklabels([r"$0$", r"$0.4$", r"$0.8$"])
    axes[2,0].get_legend().set_title(r"$E_0$ [a.u.]")
    axes[2,0].set_yticks([0, 0.05, 0.1, 0.15])
    axes[2,0].set_yticklabels([r"$0$", r"$0.8$", r"$1.6$", r"$2.4$"])

    amp = [0.001, 0.007]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    for y0 in range(2):
        with h5py.File(f"./data/mxl_144_scanamp/multimode_cavmd_144_36_{amp[y0]}_neq.h5", "r") as feq:
            eqref = {key: feq[key][:] for key in feq.keys()}
        eqref = np.sum(np.reshape(eqref["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2    
        vmax = np.max(np.max(eqref))
        vmin = vmax * 0.01 if y0 == 0 else vmax * 0.005
        for x0 in range(2):
            with h5py.File(f"./data/mxl_144_scanamp/multimode_cavmd_144_{uplist[x0]}_{amp[y0]}_neq.h5", "r") as f:
                data = {key: f[key][:] for key in f.keys()}
            sp_avg = np.sum(np.reshape(data["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2
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
                      reference_ngrid=144, scale_emsd_by_size=True):

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
    for j in range(len(amp_list)):
        mtimes = None
        if folder == "./data/mxl_n_dependence/" :
            row_data = np.load(f"{folder}/multimode_cavmd_{ngrid}_{up}_{amp_list[j]}_neq.npy")
            sp = np.mean(np.reshape(row_data, (10000, 144, -1)), axis=2)**2
            sp = sp[50:2551, :].T
        else :
            with h5py.File(f"{folder}/multimode_cavmd_{ngrid}_{up}_{amp_list[j]}_neq.h5", "r") as f:
                #sp = np.abs(f["effective_efield"][50:2551, :, 1]).T**2
                sp = np.mean(np.reshape(f["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2
                sp = sp[50:2551, :].T
        actual_ngrid = sp.shape[0]
        x_grid_1d = np.arange(1, actual_ngrid + 1) / (actual_ngrid + 1)
        x0 = x_grid_1d[np.argmin(np.abs(x_grid_1d - 0.2))]
        mtimes = 2 * (sp.shape[1] - 1) / 1000
        t = np.linspace(0, mtimes, sp.shape[1])
        if amp_list[j] in ["0.001", "0.002"]:
            x_mean = ballistic_guided_centroid(sp, t, x_grid_1d, x0, vg)
        else:
            x_mean = localization_aware_centroid(sp, t, x_grid_1d, x0, vg)
        finite = np.isfinite(x_mean)
        if not np.any(finite):
            continue
        first = np.where(finite)[0][0]
        displacement = x_mean - x_mean[first]
        if amp_list[j] in ["0.001", "0.002"]:
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
    amp = ["0.001", "0.007"]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    label_list = {0 : ["(a)","(b)"], 1: ["(c)","(d)"]}
    x_grid_1d = np.arange(1, 145) / 145
    evaluate_range = np.where((x_grid_1d >= 0.05) & (x_grid_1d <= 0.95))[0]
    for y0 in range(2):
        k_list, mean_list, std_list = [], [], []
        for idx, up in enumerate(up_list):
            sp = np.load(f"./data/mxl_144_co/coenergy_144_{up}_{amp[y0]}.npy")
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
        clp.plotone([k_parallel], [msd_array], axes[y0,1], colors=["ro-"], showlegend=False, 
                    xlim=[425, 850],
                    ylabel=r"most probable location [$\mu$m]", xlabel=r"$k_{\parallel}$ [$\rm{cm}^{-1}$]" if y0==1 else None)
        axes[y0,1].text(0.08, 0.97, label_list[y0][1], transform=axes[y0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="k")
        axes[y0,1].text(0.05, 0.15, "average over time 10-20 ps", transform=axes[y0,1].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="k")
        axes[y0,1].set_yticks([0,72, 144])
        axes[y0,1].set_yticklabels([r"$0$", r"$200$", r"$400$"])

    for y0 in range(2):
        sp_avg = np.load(f"./data/mxl_144_co/coenergy_144_36_{amp[y0]}.npy")
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
            cbar.set_ticks([2.4, 2.6, 2.8, 3, 3.2])
            cbar.set_ticklabels([r"$2.4$", r"$2.6$", r"$2.8$", r"$3.0$", r"$3.2$"])
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
        with h5py.File(f"./data/mxl_144_qc/multimode_cavmd_144_36_{amp[x0]}_neq.h5", "r") as f:
            data = {key: f[key][:] for key in f.keys()}
        sp = np.zeros((5000,144))
        for i in range(144):
            x, sp[:,i] = ir_spectrum(data['qc'][:,i,1], dt_fs)
        dx = x[2] - x[1]
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

        phe = data["photonic_energy"][:3001, :]
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

def plot_efieldy_animation():

    from pathlib import Path
    from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
    from matplotlib.patches import Rectangle

    axis_fontsize = 16
    efieldidx = 2 # 2 corresponds to E0=0.002 a.u., 5 corresponds to E0=0.005 a.u.
    input_file = Path(f"data/mxl_144*144/2d_144_144_efield_neq_0.00{efieldidx}.npy")
    output_file = Path(f"figures/efieldy_2d_144_144_0.00{efieldidx}.mp4")
    shape = (144, 144)
    frame_stride = 1
    fps = 25
    dpi = 160
    max_frames = None
    dt_fs = 20

    input_file = Path(input_file)
    output_file = Path(output_file)
    amp = np.load(input_file)**2

    ntimes, ngrid = amp.shape
    frame_ids = np.arange(0, ntimes, frame_stride)
    if max_frames is not None:
        frame_ids = frame_ids[:max_frames]

    t_ps = np.arange(ntimes) * dt_fs / 1000
    vmax = np.max(amp)
    norm = LogNorm(vmin=0.001*vmax, vmax=vmax)
    frame0 = amp[0].reshape(shape)

    fig, ax = clp.initialize(1, 1, width=8.0, height=7.0, LaTeX=True, fontsize=12, return_fig_args=True)

    image = ax.imshow(
        frame0,
        aspect="equal",
        origin="lower",
        cmap=cm.hot,
        interpolation="nearest",
        norm=norm,
        extent=(0.5, shape[1] + 0.5, 0.5, shape[0] + 0.5),
    )

    marker_center = (85, 72)
    marker_width = 30
    marker_height = 24
    marker_linewidth = 1.5
    if marker_center is None:
        marker_center = ((shape[1] + 1) / 2, (shape[0] + 1) / 2)
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
    ax.add_patch(center_marker)

    time_text = ax.text(0.98, 0.95, "", transform=ax.transAxes, fontsize=axis_fontsize, fontweight="bold", va="top", ha="right", color="w")
    ax.text(0.98, 0.08, r"$E_0=%d\times10^{-3} \ \rm{a.u.}$" % efieldidx, transform=ax.transAxes, fontsize=axis_fontsize, fontweight="bold", va="top", ha="right", color="w")
    ax.text(0.05, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=ax.transAxes, fontsize=axis_fontsize, fontweight='bold', va='top', ha='left', color="w")
    ax.text(0.05, 0.9, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=ax.transAxes, fontsize=axis_fontsize, fontweight='bold', va='top', ha='left', color="w")
    ax.set_xlabel(r"$L_x$ \ position \ [$\mu\rm{m}$]", fontsize=axis_fontsize)
    ax.set_ylabel(r"$L_y$ \ position \ [$\mu\rm{m}$]", fontsize=axis_fontsize)
    ax.set_box_aspect(shape[0] / shape[1])
    ax.set_xlim(0.5, shape[1] + 0.5)
    ax.set_ylim(0.5, shape[0] + 0.5)
    ax.set_xticks([1, max(1, shape[1] // 2), shape[1]])
    ax.set_yticks([1, max(1, shape[0] // 2), shape[0]])
    ax.set_xticklabels([r"$0$", r"$200$", r"$400$"])
    ax.set_yticklabels([r"$0$", r"$200$", r"$400$"])
    ax.tick_params(axis="both", which="major", labelsize=axis_fontsize)
    cbar_ax = fig.add_axes([0.84, 0.12, 0.025, 0.79])
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label(r"$\rm{a.u.}$", fontsize=axis_fontsize)
    cbar.ax.tick_params(labelsize=axis_fontsize)
    cbar.set_label("E-field intensity [a.u.]", fontsize=12)
    fig.subplots_adjust(left=0.13, right=0.82, bottom=0.12, top=0.91)

    def update(frame_id):
        image.set_data(amp[frame_id].reshape(shape))
        time_text.set_text(f"t = {t_ps[frame_id]:.3f} ps")
        return image, time_text, center_marker

    ani = FuncAnimation(fig, update, frames=frame_ids, interval=1000 / fps, blit=False)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.suffix.lower() == ".gif":
        writer = PillowWriter(fps=fps)
    else:
        writer = FFMpegWriter(fps=fps, bitrate=2400)
    ani.save(output_file, writer=writer, dpi=dpi)
    plt.close(fig)
    print(f"Saved animation to {output_file}")

def plot_coenergy_animation():

    from pathlib import Path
    from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
    from matplotlib.patches import Rectangle

    axis_fontsize = 16
    efieldidx = 5 # 2 corresponds to E0=0.002 a.u., 5 corresponds to E0=0.005 a.u.
    input_file = Path(f"data/mxl_144*144/coenergy_36_0.00{efieldidx}.npy")
    output_file = Path(f"figures/coenergy_2d_144_144_0.00{efieldidx}.mp4")
    shape = (144, 144)
    frame_stride = 1
    fps = 25
    dpi = 160
    max_frames = None
    dt_fs = 20

    input_file = Path(input_file)
    output_file = Path(output_file)
    amp = np.abs(np.load(input_file)) * (219474.63 / 36) * 0.01

    ntimes = amp.shape[0]
    frame_ids = np.arange(0, ntimes, frame_stride)
    if max_frames is not None:
        frame_ids = frame_ids[:max_frames]

    t_ps = np.arange(ntimes) * dt_fs / 1000
    vmin = np.percentile(amp, 60)
    vmax = np.percentile(amp, 99)
    norm = LogNorm(vmin=vmin, vmax=vmax)
    frame0 = amp[0].reshape(shape)

    fig, ax = clp.initialize(1, 1, width=8.0, height=7.0, LaTeX=True, fontsize=12, return_fig_args=True)

    image = ax.imshow(
        frame0,
        aspect="equal",
        origin="lower",
        cmap=cm.hot,
        interpolation="nearest",
        norm=norm,
        extent=(0.5, shape[1] + 0.5, 0.5, shape[0] + 0.5),
    )

    marker_center = (85, 72)
    marker_width = 30
    marker_height = 24
    marker_linewidth = 1.5
    if marker_center is None:
        marker_center = ((shape[1] + 1) / 2, (shape[0] + 1) / 2)
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
    ax.add_patch(center_marker)

    time_text = ax.text(0.98, 0.95, "", transform=ax.transAxes, fontsize=axis_fontsize, fontweight="bold", va="top", ha="right", color="w")
    ax.text(0.98, 0.08, r"$E_0=%d\times10^{-3} \ \rm{a.u.}$" % efieldidx, transform=ax.transAxes, fontsize=axis_fontsize, fontweight="bold", va="top", ha="right", color="w")
    ax.text(0.05, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=ax.transAxes, fontsize=axis_fontsize, fontweight='bold', va='top', ha='left', color="w")
    ax.text(0.05, 0.9, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=ax.transAxes, fontsize=axis_fontsize, fontweight='bold', va='top', ha='left', color="w")
    ax.set_xlabel(r"$L_x$ \ position \ [$\mu\rm{m}$]", fontsize=axis_fontsize)
    ax.set_ylabel(r"$L_y$ \ position \ [$\mu\rm{m}$]", fontsize=axis_fontsize)
    ax.set_box_aspect(shape[0] / shape[1])
    ax.set_xlim(0.5, shape[1] + 0.5)
    ax.set_ylim(0.5, shape[0] + 0.5)
    ax.set_xticks([1, max(1, shape[1] // 2), shape[1]])
    ax.set_yticks([1, max(1, shape[0] // 2), shape[0]])
    ax.set_xticklabels([r"$0$", r"$200$", r"$400$"])
    ax.set_yticklabels([r"$0$", r"$200$", r"$400$"])
    ax.tick_params(axis="both", which="major", labelsize=axis_fontsize)
    cbar_ax = fig.add_axes([0.84, 0.12, 0.025, 0.79])
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label(r"$\rm{a.u.}$", fontsize=axis_fontsize)
    cbar.ax.tick_params(labelsize=axis_fontsize)
    cbar.set_label(r"$E_{\mathrm{C=O}}$ per molecule [$10^2$ cm$^{-1}$]", fontsize=12)
    fig.subplots_adjust(left=0.13, right=0.82, bottom=0.12, top=0.91)

    def update(frame_id):
        image.set_data(amp[frame_id].reshape(shape))
        time_text.set_text(f"t = {t_ps[frame_id]:.3f} ps")
        return image, time_text, center_marker

    ani = FuncAnimation(fig, update, frames=frame_ids, interval=1000 / fps, blit=False)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    if output_file.suffix.lower() == ".gif":
        writer = PillowWriter(fps=fps)
    else:
        writer = FFMpegWriter(fps=fps, bitrate=2400)
    ani.save(output_file, writer=writer, dpi=dpi)
    plt.close(fig)
    print(f"Saved animation to {output_file}")

def plot_2d_coenergy():

    fig, axes = clp.initialize(2, 3, width=3*4.3, height=4.3*2*0.86, LaTeX=True, fontsize=12, return_fig_args=True)
    data_list = ["./data/mxl_144*144/coenergy_36_0.002.npy", "./data/mxl_144*144/coenergy_36_0.005.npy"]
    time_list = [50, 150, 300]
    time_labels = [r"$t=%.1f$ ps" % (t/50) for t in time_list]
    label_list = {0 : ["(a)","(b)","(c)"], 1: ["(d)","(e)","(f)"]}
    amp_list = [r"$E_0=2\times10^{-3}$ a.u.", r"$E_0=5\times10^{-3}$ a.u."]

    for x0 in range(2):
        sp = np.load(data_list[x0]) * (219474.63 / 36) * 0.01
        vmin = np.percentile(sp[time_list[0]], 80)
        vmax = np.percentile(sp[time_list[0]], 99)
        for y0 in range(3):
            extent = [0 , 144, 0, 144]
            pos = axes[x0, y0].imshow(sp[time_list[y0]],
                     aspect='equal', extent=extent, origin="lower",
                    cmap=cm.hot, interpolation='nearest',
                    norm=LogNorm(vmin=vmin, vmax=vmax))
            axes[x0, y0].set_box_aspect(1)
            clp.plotone([], [], axes[x0, y0], colors=["c-","c-"], ylabel=r"$L_y$ \ position \ [$\mu\rm{m}$]" if y0==0 else None, xlabel=r"$L_x$ \ position \ [$\mu\rm{m}$]" if x0==1 else None, showlegend=False)
            axes[x0, y0].set_yticks([0,72, 144])
            axes[x0, y0].set_yticklabels([r"$0$", r"$200$", r"$400$"])
            axes[x0, y0].set_xticks([0,72, 144])
            axes[x0, y0].set_xticklabels([r"$0$", r"$200$", r"$400$"])
            axes[x0, y0].text(0.98, 0.18, time_labels[y0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[x0, y0].text(0.08, 0.98, label_list[x0][y0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")
            axes[x0, y0].text(0.98, 0.08, amp_list[x0], transform=axes[x0, y0].transAxes, fontsize=12, fontweight='bold', va='top', ha='right', color="w")

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
            axes[x0, y0].add_patch(center_marker)
            if y0 == 2:
                cbar_ax = axes[x0, y0].inset_axes([1.04, 0.0, 0.055, 1.0])
                cbar = fig.colorbar(pos, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=12)
                cbar.set_label(r"$E_{\mathrm{C=O}}$ per molecule [$10^2$ cm$^{-1}$]", fontsize=12)
                if x0 == 0:
                    cbar.set_ticks([3, 4])
                    cbar.set_ticklabels([r"$3$", r"$4$"])
                else:
                    cbar.set_ticks([3,4,5,6,7])
                    cbar.set_ticklabels([r"$3$", r"$4$", r"$5$", r"$6$", r"$7$"])

    for j in range(2):
        axes[j, 0].text(0.6, 0.95, r"$k_{\parallel} \ =$"+rf"$ \ 450 \ $"+r"$\rm{cm}^{-1}$", transform=axes[j, 0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")
        axes[j, 0].text(0.6, 0.85, r"$W_{\rm ph}=$"+rf"$ \ 0.61 \ $", transform=axes[j, 0].transAxes, fontsize=12, fontweight='bold', va='top', ha='left', color="w")

    plt.rcParams["axes.axisbelow"] = False
    fig.subplots_adjust(left=0.06, right=0.94, bottom=0.10, top=0.94, wspace=0.18, hspace=0.20)
    fig_width, fig_height = fig.get_size_inches()
    image_group_shift = 0.055
    image_col_gap = 0.027
    for row in range(2):
        base_pos = axes[row, 1].get_position()
        square_width = base_pos.height * fig_height / fig_width
        for col in range(3):
            axes[row, col].set_position([
                base_pos.x0 + image_group_shift + (col - 1) * (square_width + image_col_gap),
                base_pos.y0,
                square_width,
                base_pos.height,
            ])
    clp.adjust(savefile=f"./figures/figS1_2d_coenergy.png")

def plot_vg():
    fig, axes = clp.initialize(1, 2, width=4.3*2, height=4.3*0.618*1.1, LaTeX=True, fontsize=12, return_fig_args=True)
    with h5py.File(f"./data/mxl_144*144/multimode_cavmd_144_eq.h5", "r") as f:
        data = {key: f[key][:] for key in f.keys()}
    sp_full = np.zeros((5000,144))
    for i in range(144):
        x, sp_full[:,i] = ir_spectrum(data['qc'][:,i,1], 2)
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
    clp.adjust(savefile=f"./figures/figS2_vg.png")

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
        mt_list, mmsd_list = get_adaptive_emsd(uplist[j], coeff[j], amp_list=["0.001", "0.003", "0.005", "0.007"])
        clp.plotone(mt_list, mmsd_list, axes[2,j], ylabel=r"$\mathrm{MSD}$ \ [$10^{4} \ \mu\rm{m}^2$]" if j == 0 else None, xlabel="time [ps]",
                    showlegend=True if j==0 else False, legendloc=(0.25,0.3), legendFontSize=8, 
                    xlim=[0,mt_list[0][-1]], ylim=[0,0.05*(j+1)],
                    labels=amp_list, colorMap=plt.cm.hot, colorMap_endpoint=0.6, alpha=0.4)
        quad_t = mt_list[0][(mt_list[0] >= 1.0) & (mt_list[0] <= 4.5)]
        linear_t = mt_list[0][(mt_list[0] >= 2.0) & (mt_list[0] <= 5.0)]
        quad_fit_params = [[0.003969807815718676, -0.002321773431376476, 0.0006899936493188445],
                           [0.006300141003230975, -0.008397181870203245, 0.004568356292778106],
                           [0.01990514613643173, -0.023893159803765197, 0.015850866748061753]]
        linear_fit_params = [[0.016093165436323338, -0.020920098506520733],
                             [0.039898507445015934, -0.08864441131004144],
                             [0.13881807081968064, -0.3094249434532433]]
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
        with h5py.File(f"./data/mxl_144_scanamp/multimode_cavmd_144_24_{amp[y0]}_neq.h5", "r") as feq:
            eqref = {key: feq[key][:] for key in feq.keys()}
        eqref = np.sum(np.reshape(eqref["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2    
        vmax = np.max(np.max(eqref))
        vmin = vmax * 0.01 if y0 == 0 else vmax * 0.005
        for x0 in range(3):
            with h5py.File(f"./data/mxl_144_scanamp/multimode_cavmd_144_{uplist[x0]}_{amp[y0]}_neq.h5", "r") as f:
                data = {key: f[key][:] for key in f.keys()}
            sp_avg = np.sum(np.reshape(data["effective_efield"][:, :, 1], (10000, 144, -1)), axis=2)**2
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
    clp.adjust(savefile=f"./figures/figS3_1d_efield_else.png")

def plot_1d_efield_ndependence():

    axes = clp.initialize(2, 4, width=4.3*4, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12)
    label_list = {0 : ["(a)","(b)","(c)","(d)"], 1: ["(e)","(f)","(g)","(h)"]}
    amp_list = [r"$7\times10^{-3}$"]
    n_list = [144, 288, 576, 1152]

    for j in range(4):
        mt_list, mmsd_list = get_adaptive_emsd(36, 0.11638023, amp_list=["0.007"], folder="./data/mxl_n_dependence/", ngrid=n_list[j])
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

    refdata = np.load("./data/mxl_n_dependence/multimode_cavmd_1152_36_0.007_neq.npy")
    refsp = np.mean(np.reshape(refdata, (10000, 144, -1)), axis=2)**2
    vmax = np.max(np.max(refsp))
    vmin = vmax * 0.01
    for x0 in range(4):
        data = np.load(f"./data/mxl_n_dependence/multimode_cavmd_{n_list[x0]}_36_0.007_neq.npy")
        sp_avg = np.mean(np.reshape(data, (10000, 144, -1)), axis=2)**2
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
    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS4_1d_efield_numberdependence.png")

def plot_addition_mmsd():

    axes = clp.initialize(2, 4, width=4*4.3, height=4.3*0.618*1.1*2, LaTeX=True, fontsize=12, sharex=True, sharey=True)
    up_list = np.array([24, 27, 30, 48])
    k_parallel_list = [300, 337.5, 375, 600]
    amp = ["0.001", "0.007"]
    amp_list = [r"$E_0=1\times10^{-3}$ a.u.", r"$E_0=7\times10^{-3}$ a.u."]
    label_list = {0 : ["(a)","(b)","(c)","(d)"], 1 : ["(e)","(f)","(g)","(h)"]}

    for y0 in range(2):
        for x0 in range(4):
            sp_avg = np.load(f"./data/mxl_144_co/coenergy_144_{up_list[x0]}_{amp[y0]}.npy")
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
                cbar.set_ticks([2.5, 3, 3.5, 4])
                cbar.set_ticklabels([r"$2.5$", r"$3.0$", r"$3.5$", r"$4.0$"])
            else:
                cbar.set_ticks([1, 2, 3, 4])
                cbar.set_ticklabels([r"$1.0$", r"$2.0$", r"$3.0$", r"$4.0$"])
            cbar.ax.tick_params(labelsize=12)
            if x0 == 3: cbar.set_label(r"C=O energy [$10^{%d} \ $cm$^{-1}$]"%(y0+2), fontsize=12)
            axes[y0,x0].hlines(144*0.05, 0, 20, colors="c", linestyles="--", lw=1.5)
            axes[y0,x0].hlines(144*0.35, 0, 20, colors="c", linestyles="--", lw=1.5)

    plt.rcParams["axes.axisbelow"] = False
    plt.tight_layout()
    clp.adjust(savefile=f"./figures/figS5_addition_mmsd.png")

if __name__ == "__main__":
    plot_2d_efield_final()
    plot_1d_efield()
    plot_mmsd()
    plot_neq_spectrum_with_model()
    #plot_efieldy_animation()
    #plot_coenergy_animation()
    plot_2d_coenergy()
    plot_vg()
    plot_1d_efield_else()
    plot_1d_efield_ndependence()
    plot_addition_mmsd()
