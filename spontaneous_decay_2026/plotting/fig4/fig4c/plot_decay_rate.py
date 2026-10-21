import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
plt.rcParams.update({
    # ----------- 基础字体 -----------
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],       # Nature 通常使用 Arial（无衬线字体）
    'font.size': 8,                     # 正文字号通常为 8 pt（Nature 图要求小而清晰）

    # ----------- 坐标轴 -----------
    'axes.linewidth': 0.8,              # 坐标轴线条宽度
    'axes.labelsize': 10,                # 轴标签字号
    'axes.titlesize': 10,                # 标题字号
    'axes.labelpad': 2,                 # 标签与轴之间的距离
    'axes.titlepad': 4,                 # 标题与图之间的距离
    'axes.labelcolor': 'black',
    'axes.edgecolor': 'black',

    # ----------- 刻度 -----------
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

    # ----------- 线条与图形 -----------
    'lines.linewidth': 0.8,
    'lines.markersize': 3,
    'lines.markeredgewidth': 0.3,
    'patch.linewidth': 0.6,

    # ----------- 图例 -----------
    'legend.frameon': False,
    'legend.fontsize': 10,
    'legend.handlelength': 1.2,
    'legend.handletextpad': 0.3,
    'legend.borderaxespad': 0.2,
    'legend.labelspacing': 0.2,

    # ----------- 网格与背景 -----------
    'axes.grid': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',

    # ----------- 图像尺寸与保存 -----------
    'figure.figsize': [3.15, 2.36],     # 单列图大小（约 8 cm × 6 cm）
    'savefig.dpi': 600,                 # 高分辨率导出
    'savefig.bbox': 'tight',            # 紧凑布局
    'savefig.pad_inches': 0.02,
    'savefig.transparent': True,

    # ----------- 数学字体 -----------
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
})

rescaling_factor = 1e1
def plot_from_file():
    data = np.loadtxt('rate')
    cavity_frequency = data[:,0]
    decay_rate = data[:,1]
    #decay_rate = decay_rate / rescaling_factor**2 * 1e3
    #plt.figure(figsize=(3.15,2.36))
    plt.figure(figsize=(2.86,2.36))
    """
    ax1.plot(
        amplitude,
        decay_rate_vacuum,
        '-o',
        color='black',
        markerfacecolor='none',
        markeredgewidth=1.5,
        linewidth=1.5,
        markersize=5,
        label='in vacuum',
        zorder=5
    )"""
    plt.plot(cavity_frequency, decay_rate, '-o', color='slateblue', markersize=5, markerfacecolor='none', markeredgewidth=1.5, linewidth=1.5, label="cavity", zorder=10, clip_on=False)
    plt.axhline(y=0.088513, linestyle='--', linewidth=1.5, color='grey', label='vacuum')
    plt.axvline(x=2.591332, linestyle='--', linewidth=1.5, color='coral', label=r"$\bar{\omega}_{\rm m}$")
    plt.xlabel(r"$\omega_{\mathrm{c}}$ [eV]")
    plt.ylabel(r"decay rate [ns$^{-1}$]")
    #plt.title("C2H4 TDDFT Spectrum Comparison")
    plt.xlim(2.0, 3.2)
    plt.xticks(np.arange(2.0, 3.2, 0.4))
    plt.ylim(0.0, 0.225)
    plt.yticks(np.arange(0.0, 0.225, 0.05))
    #plt.legend(bbox_to_anchor=(0.60, 1.00))
    plt.legend(loc=1, borderaxespad=0.0)
    plt.tight_layout()
    plt.savefig('decay_rate.png', dpi=400)
    plt.show()

if __name__ == "__main__":
    plot_from_file()
