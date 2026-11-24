"""
Figuras de verificacao: Fortran x Python (series diarias e erro relativo ao
longo do tempo) para vazao e QSS nos postos.

Uso: python verification/plots.py <dir_fortran> <dir_python> <saida.png>
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"


def load(path, skip_hdr=1):
    with open(path, encoding="latin-1") as f:
        hdr = f.readline().split()
    a = np.loadtxt(path, skiprows=skip_hdr)
    return hdr[1:], a[:, 1:]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def main(fd, pd, out):
    codes, Qf = load(f"{fd}/VAZAO.HIG")
    _, Qp = load(f"{pd}/VAZAO.HIG")
    _, Sf = load(f"{fd}/MINI_QSS.txt")
    _, Sp = load(f"{pd}/MINI_QSS.txt")
    t = 1974 + np.arange(Qf.shape[0]) / 365.25
    k = len(codes) - 1                      # posto mais a jusante
    fig, axs = plt.subplots(2, 2, figsize=(12, 7), dpi=130)
    for row, (F, P, name, unit) in enumerate(((Qf, Qp, "Vazao", "m3/s"), (Sf, Sp, "QSS", "t/dia"))):
        ax = axs[row, 0]
        ax.plot(t, F[:, k], color=BLUE, lw=1.6, label="Fortran (original)")
        ax.plot(t, P[:, k], color=ORANGE, lw=0.8, label="Python")
        ax.set_title(f"{name} na minibacia {codes[k]} ({unit})", color=INK, fontsize=10, loc="left")
        ax.legend(frameon=False, fontsize=8)
        style(ax)
        ax = axs[row, 1]
        floor = 1.0 if name == "Vazao" else 10.0
        rel = np.abs(P - F) / np.maximum(np.abs(F), floor)
        ax.semilogy(t, np.maximum(rel.max(axis=1), 1e-9), color=BLUE, lw=0.8, label="maximo entre postos")
        ax.semilogy(t, np.maximum(np.median(rel, axis=1), 1e-9), color=ORANGE, lw=0.8, label="mediana entre postos")
        ax.set_title(f"{name}: erro relativo |Py-F|/max(|F|,{floor:g}) por dia", color=INK, fontsize=10, loc="left")
        ax.legend(frameon=False, fontsize=8)
        style(ax)
    fig.tight_layout()
    fig.savefig(out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
