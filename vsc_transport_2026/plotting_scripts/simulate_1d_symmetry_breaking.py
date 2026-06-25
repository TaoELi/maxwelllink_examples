"""
One-dimensional polariton spectrum with spatial translational-symmetry breaking.

This script is extracted from PolaritonSpectrumCalcUpdate.py and keeps only the
1D Gaussian case.  The model follows the coupled exciton-photon Hamiltonian used
for in-plane translational symmetry breaking in 064308_1_5.0209212.pdf:

    H = sum_x omega_m |x><x| + sum_k omega_k |k><k|
        + sum_{x,k} [g_k(x) |x><k| + c.c.]

with

    omega_k = sqrt(omega_perp^2 + omega_parallel(k)^2),
    g_k(x) = g0 sqrt(rho(x)) exp(i k x),
    rho(x) = density profile centered at center.

The non-uniform density rho(x) breaks 1D translational symmetry, so different
in-plane photon momenta can mix through the non-uniform molecular distribution.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

try:
    import cupy as cp

    CUPY_AVAILABLE = True
except ImportError:
    cp = np
    CUPY_AVAILABLE = False


@dataclass
class Cavity1DParams:
    ngrid: int = 1080
    nmode: int = 200
    omega_perp: float = 2320.0
    omega_m: float = 2320.0
    domega: float | None = None
    g0: float | None = None
    sigma: float = 0.15
    center: float = 0.50
    density_profile: str = "gaussian"
    density_depth: float = 0.90

    def resolved(self) -> "Cavity1DParams":
        if self.domega is None:
            self.domega = 2000.0 / self.nmode
        if self.g0 is None:
            self.g0 = 2.0 * np.sqrt(1080.0 / self.ngrid)
        return self


class GaussianSymmetryBreaking1D:
    def __init__(self, params: Cavity1DParams):
        self.params = params.resolved()
        self.x_grid = np.arange(self.params.ngrid, dtype=float) / self.params.ngrid
        self.k_grid = self._build_k_grid(self.params.nmode)
        self.omega_parallel = np.abs(self.k_grid / (2.0 * np.pi) * self.params.domega)
        self.omega_k = np.sqrt(self.params.omega_perp**2 + self.omega_parallel**2)
        self.density = self._build_density()
        self.g_distribution = np.sqrt(self.density)
        self.coupling = self._build_coupling()

        self.translation_error = self._check_plane_wave_orthogonality()
        self.hamiltonian = self._build_hamiltonian()
        self.eigenvalues, self.eigenvectors = self._diagonalize()
        self.photon_weight = self._photon_weight()

    @staticmethod
    def _build_k_grid(nmode: int) -> np.ndarray:
        if nmode % 2 != 0:
            raise ValueError("nmode must be even to match the +/- k grid used in the source script.")
        mode_numbers = np.arange(-nmode // 2, nmode // 2, dtype=float)
        return 2.0 * np.pi * mode_numbers

    def _build_density(self) -> np.ndarray:
        sigma = self.params.sigma
        center = self.params.center
        if sigma <= 0.0:
            raise ValueError("sigma must be positive.")

        if self.params.density_profile == "gaussian":
            gaussian = np.exp(-0.5 * (self.x_grid - center) ** 2 / sigma**2)
            density = gaussian / (sigma * np.sqrt(2.0 * np.pi))
            density /= np.mean(density)
        elif self.params.density_profile == "inverted_gaussian":
            gaussian = np.exp(-0.5 * (self.x_grid - center) ** 2 / sigma**2)
            depth = self.params.density_depth
            if not 0.0 <= depth < 1.0:
                raise ValueError("density_depth must satisfy 0 <= density_depth < 1 for inverted_gaussian.")
            gaussian_dip = (gaussian - np.min(gaussian)) / (np.max(gaussian) - np.min(gaussian))
            density = 1.0 - depth * gaussian_dip
        elif self.params.density_profile == "inverted_square":
            half_width = sigma
            depth = self.params.density_depth
            if center - half_width < 0.0 or center + half_width > 1.0:
                raise ValueError("sigma is the inverted_square half-width and center +/- sigma must stay within [0, 1].")
            if not 0.0 <= depth < 1.0:
                raise ValueError("density_depth must satisfy 0 <= density_depth < 1 for inverted_square.")
            square = np.abs(self.x_grid - center) <= half_width
            density = np.ones_like(self.x_grid)
            density[square] = 1.0 - depth
        else:
            raise ValueError("density_profile must be 'gaussian', 'inverted_gaussian', or 'inverted_square'.")
        return density

    def _density_label(self) -> str:
        return self.params.density_profile.replace("_", " ")

    def _build_coupling(self) -> np.ndarray:
        phase = np.exp(1j * np.outer(self.k_grid, self.x_grid))
        return self.params.g0 * phase * self.g_distribution[np.newaxis, :]

    def _check_plane_wave_orthogonality(self) -> float:
        phase = np.exp(1j * np.outer(self.k_grid, self.x_grid))
        overlap = phase @ np.conjugate(phase.T) / self.params.ngrid
        return float(np.mean(np.abs(overlap - np.eye(self.params.nmode))))

    def _build_hamiltonian(self) -> np.ndarray:
        start = time.perf_counter()
        ngrid = self.params.ngrid
        nmode = self.params.nmode
        size = ngrid + nmode
        hamiltonian = np.zeros((size, size), dtype=np.complex128)
        hamiltonian[:ngrid, :ngrid] = np.eye(ngrid) * self.params.omega_m
        hamiltonian[ngrid:, ngrid:] = np.diag(self.omega_k)
        hamiltonian[:ngrid, ngrid:] = self.coupling.T
        hamiltonian[ngrid:, :ngrid] = np.conjugate(self.coupling)
        print(f"Hamiltonian constructed in {time.perf_counter() - start:.3f} s")
        return hamiltonian

    def _diagonalize(self) -> tuple[np.ndarray, np.ndarray]:
        start = time.perf_counter()
        hamiltonian = cp.asarray(self.hamiltonian, dtype=cp.complex64)
        eigenvalues, eigenvectors = cp.linalg.eigh(hamiltonian)
        if CUPY_AVAILABLE:
            eigenvalues = cp.asnumpy(eigenvalues)
            eigenvectors = cp.asnumpy(eigenvectors)
        print(f"Diagonalization finished in {time.perf_counter() - start:.3f} s")
        return eigenvalues, eigenvectors

    def _photon_weight(self) -> np.ndarray:
        weight = np.abs(self.eigenvectors) ** 2
        return np.sum(weight[self.params.ngrid :, :], axis=0)

    def make_spectrum(
        self,
        resolution: float = 1.0,
        linewidth: float = 5.0,
        omega_min: float | None = None,
        omega_max: float | None = None,
    ) -> tuple[np.ndarray, tuple[float, float, float, float], float]:
        omega_parallel_min = float(np.min(self.omega_parallel))
        omega_parallel_max = float(np.max(self.omega_parallel))
        if omega_min is None:
            omega_min = float(np.min(self.eigenvalues) - 60.0)
        if omega_max is None:
            omega_max = float(np.max(self.eigenvalues) + 30.0)

        nx = int(np.rint((omega_parallel_max - omega_parallel_min) / self.params.domega)) + 1
        ny = int(np.ceil((omega_max - omega_min) / resolution)) + 1
        spectrum = np.zeros((nx, ny))

        photon_weight_by_mode = np.abs(self.eigenvectors[self.params.ngrid :, :]) ** 2
        x_bins = np.clip(
            np.rint((self.omega_parallel - omega_parallel_min) / self.params.domega).astype(int),
            0,
            nx - 1,
        )

        if linewidth > 0.0:
            omega_axis = omega_min + np.arange(ny) * resolution
            lineshape = np.exp(-0.5 * ((omega_axis[:, np.newaxis] - self.eigenvalues[np.newaxis, :]) / linewidth) ** 2)
            lineshape /= np.sum(lineshape, axis=0, keepdims=True)
            for mode_index, x_bin in enumerate(x_bins):
                spectrum[x_bin] += photon_weight_by_mode[mode_index] @ lineshape.T
        else:
            y_bins = np.clip(
                np.rint((self.eigenvalues - omega_min) / resolution).astype(int),
                0,
                ny - 1,
            )
            for mode_index, x_bin in enumerate(x_bins):
                np.add.at(spectrum[x_bin], y_bins, photon_weight_by_mode[mode_index])

        rabi = self._estimate_rabi(spectrum, resolution)
        spectrum += 1e-8
        spectrum = spectrum.T
        extent = (omega_parallel_min, omega_parallel_max, omega_min, omega_max)
        return spectrum, extent, rabi

    @staticmethod
    def _estimate_rabi(spectrum: np.ndarray, resolution: float) -> float:
        if spectrum.shape[0] < 2:
            return 0.0
        resonance_row = spectrum[1]
        first = int(np.argmax(resonance_row))
        row_without_first = np.array(resonance_row)
        row_without_first[first] = 0.0
        second = int(np.argmax(row_without_first))
        return abs(first - second) * resolution

    def plot(self, spectrum: np.ndarray, extent: tuple[float, float, float, float], output: Path) -> None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4.3, 4.3 * 0.618 * 1.25))
        fig.subplots_adjust(left=0.16, right=0.86, bottom=0.18, top=0.95)
        spectrum_xmax = min(extent[1], 1800.0)
        x_ticks = np.linspace(extent[0], spectrum_xmax, 4)

        vmax = float(np.max(spectrum))
        vmin = vmax * 1e-4
        image = ax.imshow(
            spectrum,
            aspect="auto",
            extent=extent,
            origin="lower",
            interpolation="nearest",
            cmap=cm.hot,
            norm=LogNorm(vmin=vmin, vmax=vmax),
        )

        omega_p = np.linspace(extent[0], extent[1], spectrum.shape[1])
        omega_m = np.ones_like(omega_p) * self.params.omega_m
        omega_c = np.sqrt(omega_p**2 + self.params.omega_perp**2)
        ax.plot(omega_p, omega_m, "w--", lw=1.2)
        ax.plot(omega_p, omega_c, "g--", lw=1.2)

        ax.set_xlim(extent[0], spectrum_xmax)
        ax.set_xlabel(r"$\omega_{\parallel}$ [cm$^{-1}$]")
        ax.set_ylabel(r"$\omega_{\rm pol}$ [cm$^{-1}$]")
        ax.set_xticks(x_ticks)
        ax.set_yticks([2200.0, 2400.0, 2600.0, 2800.0, 3000.0])
        ax.tick_params(color="c", labelsize="medium", width=2, direction="in", bottom=True, top=True, left=True, right=True)

        ax.text(0.96, 0.72, "cavity photon", transform=ax.transAxes, color="g", fontsize=12, ha="right")
        ax.text(
            0.96,
            0.36,
            "molecular mode",
            transform=ax.transAxes,
            color="w",
            fontsize=12,
            ha="right",
        )
        ax.text(
            0.96,
            0.94,
            rf"{self._density_label()}, $\sigma={self.params.sigma:g}$",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="right",
            color="c",
        )
        plt.rcParams["axes.axisbelow"] = False

        cbar_ax = fig.add_axes([0.88, 0.18, 0.025, 0.77])
        fig.colorbar(image, cax=cbar_ax)
        fig.savefig(output, dpi=220)
        plt.close(fig)

    def plot_density_distribution(self, output: Path) -> None:
        output.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4.3, 4.3 * 0.618 * 1.25))
        fig.subplots_adjust(left=0.16, right=0.96, bottom=0.18, top=0.92)

        ax.plot(self.x_grid, self.density, color="tab:blue", lw=2.0, label=r"$\rho(x)$")
        ax.plot(
            self.x_grid,
            self.g_distribution,
            color="tab:red",
            lw=1.5,
            ls="--",
            label=r"$\sqrt{\rho(x)}$",
        )
        ax.axhline(1.0, color="0.35", lw=1.0, ls=":", label=r"$\rho=1$")
        ax.axvline(self.params.center, color="0.5", lw=1.0, ls=":")

        ax.set_xlabel(r"$x/L_x$")
        ax.set_ylabel("density distribution")
        ax.set_xlim(float(np.min(self.x_grid)), float(np.max(self.x_grid)))
        ax.set_ylim(bottom=0.0)
        ax.tick_params(labelsize="medium", width=1.5, direction="in", bottom=True, top=True, left=True, right=True)
        ax.legend(frameon=False, loc="best", fontsize=10)
        ax.text(
            0.04,
            0.96,
            rf"{self._density_label()}, $\sigma={self.params.sigma:g}$",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )
        fig.savefig(output, dpi=220)
        plt.close(fig)

    def save_plot_data(self, spectrum: np.ndarray, extent: tuple[float, float, float, float], output: Path) -> None:
        output.parent.mkdir(parents=True, exist_ok=True)
        omega_axis = np.linspace(extent[2], extent[3], spectrum.shape[0])
        omega_parallel_axis = np.linspace(extent[0], extent[1], spectrum.shape[1])
        omega_parallel_12p5_ir = spectrum[:, 0]
        np.savez(
            output,
            spectrum=spectrum,
            extent=np.array(extent),
            omega_axis=omega_axis,
            omega_parallel_axis=omega_parallel_axis,
            omega_parallel_12p5_ir=omega_parallel_12p5_ir,
            first_up_spectrum=omega_parallel_12p5_ir,
            omega_m=self.params.omega_m,
            omega_perp=self.params.omega_perp,
            sigma=self.params.sigma,
            g0=self.params.g0,
            density_profile=self.params.density_profile,
            density_depth=self.params.density_depth,
            x_grid=self.x_grid,
            density=self.density,
            g_distribution=self.g_distribution,
        )


def main() -> None:
    # Modify simulation parameters here.
    base_params = dict(
        ngrid=288,
        nmode=288,
        omega_perp=2320.0,
        omega_m=2320.0,
        domega=12.5,
        g0=5.625,
        sigma=0.07,
        center=0.20,
        density_depth=0.70,
    )
    resolution = 1.0
    linewidth = 2.0

    density_profile = "inverted_gaussian"
    params = Cavity1DParams(**base_params, density_profile=density_profile)
    output_tag = params.density_profile
    figure = Path(f"figures/1d_{output_tag}_symmetry_breaking.png")
    density_figure = Path(f"figures/1d_{output_tag}_density_distribution.png")
    plot_data = Path(f"data/mxl_144_qc/model_1d_{output_tag}_spectrum.npz")

    print(f"\nRunning {params.density_profile} density profile")
    calc = GaussianSymmetryBreaking1D(params)
    spectrum, extent, rabi = calc.make_spectrum(resolution=resolution, linewidth=linewidth)
    normalized_spectrum = spectrum / np.max(spectrum)
    calc.plot(normalized_spectrum, extent, figure)
    calc.plot_density_distribution(density_figure)
    calc.save_plot_data(normalized_spectrum, extent, plot_data)

    print(f"Plane-wave orthogonality mean error: {calc.translation_error:.3e}")
    print(f"Estimated Rabi splitting: {rabi:.3f} cm^-1")
    print(f"Saved figure: {figure}")
    print(f"Saved density figure: {density_figure}")
    print(f"Saved plot data: {plot_data}")


if __name__ == "__main__":
    main()
