"""
Metricas de QSS simulada x observada para uma ou mais rodadas.

Uso:
    python verification/metricas_qss.py <QSS_observado.txt> <rodada1_dir> [<rodada2_dir> ...] [--ini 365]

Le MINI_QSS.txt de cada rodada (t/dia, colunas = minibacias dos postos) e o arquivo
observado (formato de entrada do MGB-SED: DIA MES ANO + uma coluna por posto, -1 = falha).
Considera os dias a partir de --ini (o MGB-SED usa 365 em SED_Fob_Calib), como no Fortran.

Metricas por posto (apenas dias com observacao >= 0):
  n, NSE, KGE (Gupta et al. 2009), PBIAS (%), NSE(log10), razao media sim/obs.
"""
import argparse
import os

import numpy as np


def load_sim(d):
    p = os.path.join(d, "MINI_QSS.txt")
    with open(p, encoding="latin-1") as f:
        codes = f.readline().split()[1:]
    a = np.loadtxt(p, skiprows=1)
    return codes, a[:, 1:]


def load_obs(p):
    a = np.loadtxt(p, skiprows=1)
    return a[:, 3:]


def metrics(o, s):
    m = o >= 0
    o = o[m].astype(float); s = s[m].astype(float)
    n = o.size
    if n < 3:
        return dict(n=n)
    nse = 1 - ((s - o) ** 2).sum() / ((o - o.mean()) ** 2).sum()
    r = np.corrcoef(o, s)[0, 1]
    kge = 1 - np.sqrt((r - 1) ** 2 + (s.std() / o.std() - 1) ** 2 + (s.mean() / o.mean() - 1) ** 2)
    pbias = 100 * (s.sum() - o.sum()) / o.sum()
    eps = 1.0
    lo, ls = np.log10(o + eps), np.log10(s + eps)
    nsel = 1 - ((ls - lo) ** 2).sum() / ((lo - lo.mean()) ** 2).sum()
    return dict(n=n, NSE=nse, KGE=kge, PBIAS=pbias, NSElog=nsel, r=r, ratio=s.mean() / o.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("obs")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--ini", type=int, default=365)
    a = ap.parse_args()
    obs = load_obs(a.obs)
    sims = [load_sim(d) for d in a.runs]
    codes = sims[0][0]
    i0 = a.ini - 1
    print("posto  " + "  ".join(f"{os.path.basename(d.rstrip('/')):>44s}" for d in a.runs))
    print("       " + "  ".join(f"{'n':>5s}{'NSE':>8s}{'KGE':>8s}{'PBIAS%':>9s}{'NSElog':>8s}{'sim/obs':>8s}" for _ in a.runs))
    for k, c in enumerate(codes):
        line = f"{c:>6s} "
        for _, S in sims:
            mt = metrics(obs[i0:, k], S[i0:, k])
            if "NSE" not in mt:
                line += f"{mt['n']:>5d}" + " " * 41 + "  "
            else:
                line += (f"{mt['n']:>5d}{mt['NSE']:8.3f}{mt['KGE']:8.3f}{mt['PBIAS']:9.1f}"
                         f"{mt['NSElog']:8.3f}{mt['ratio']:8.2f}  ")
        print(line)


if __name__ == "__main__":
    main()
