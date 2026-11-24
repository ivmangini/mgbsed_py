"""
Comparacao das saidas da versao Python com as saidas originais do Fortran.

Uso:
    python verification/compare.py <dir_fortran> <dir_python> [--nrows N] [--json out.json]

Para cada arquivo, le as colunas numericas (ignora colunas de data/indice) em
blocos (arquivos de ~470 MB) e calcula:
  * erro absoluto maximo e medio
  * erro relativo (|py-f|/max(|f|, piso)) - mediana, p99, maximo
  * correlacao, NSE do Python em relacao ao Fortran, vies de soma (%)
  * fracao de valores identicos na precisao gravada
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

# (arquivo, linhas de cabecalho, numero de colunas iniciais a ignorar, pasta alternativa)
FILES = [
    ("VAZAO.HIG", 1, 1),
    ("NOSOLO.txt", 0, 1),
    ("QPROP.HIG", 0, 1),
    ("BACIAS.HIG", 0, 0),
    ("VAZAO_MINI.HIG", 1, 1),
    ("QSUP.txt", 1, 4),
    ("QCEL.txt", 1, 4),
    ("MUSLE_TOTAL.txt", 1, 4),
    ("MUSLE_areia.txt", 1, 4),
    ("MUSLE_silte.txt", 1, 4),
    ("MUSLE_argila.txt", 1, 4),
    ("CONC_areia.txt", 1, 3),
    ("CONC_silte.txt", 1, 3),
    ("CONC_argila.txt", 1, 3),
    ("DEPdin_RIO.txt", 0, 3),
    ("MINI_areia.txt", 1, 1),
    ("MINI_silte.txt", 1, 1),
    ("MINI_argila.txt", 1, 1),
    ("MINI_QSS.txt", 1, 1),
    ("SEDRIO_SAIDA.txt", 1, 1),
    ("SEDRIO_DEPOSITO.txt", 1, 1),
    ("SEDRIO_EROSAO.txt", 1, 1),
    ("BAL_PLAN.txt", 1, 1),
    ("AJUSTESED.fob", 0, 0),
    ("AJUSTE.fob", 0, 0),
]


def find(dirs, name):
    for d in dirs:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    return None


def _parse(line, skip):
    t = line.split()
    out = []
    for x in t[skip:]:
        try:
            out.append(float(x))
        except ValueError:
            out.append(np.nan)          # '*****' ou texto
    return out


def iter_blocks(path, header, skip, nrows=None, block=500):
    with open(path, "r", encoding="latin-1") as f:
        for _ in range(header):
            f.readline()
        buf = []
        n = 0
        for line in f:
            if not line.strip():
                continue
            buf.append(_parse(line, skip))
            n += 1
            if len(buf) >= block:
                yield buf
                buf = []
            if nrows is not None and n >= nrows:
                break
        if buf:
            yield buf


def compare_file(pf, pp, header, skip, nrows=None, floor=None):
    st = dict(n=0, nnan_f=0, nnan_p=0, maxabs=0.0, sumabs=0.0, sf=0.0, sp=0.0,
              sff=0.0, spp=0.0, sfp=0.0, sse=0.0, nexact=0, rel=[], rows=0, maxabs_loc=None,
              ncols_f=None, ncols_p=None)
    gf = iter_blocks(pf, header, skip, nrows)
    gp = iter_blocks(pp, header, skip, nrows)
    row0 = 0
    for bf, bp in zip(gf, gp):
        # linhas "quebradas" (formatos com quebra de registro) -> achata
        try:
            F = np.array(bf, dtype=float); P = np.array(bp, dtype=float)
        except ValueError:
            F = np.concatenate([np.array(r, float) for r in bf])[None, :]
            P = np.concatenate([np.array(r, float) for r in bp])[None, :]
        if F.shape != P.shape:
            return dict(error=f"shape diferente {F.shape} vs {P.shape} (bloco na linha {row0})")
        st["ncols_f"] = F.shape[1]
        st["rows"] += F.shape[0]
        nf = np.isnan(F); npn = np.isnan(P)
        st["nnan_f"] += int(nf.sum()); st["nnan_p"] += int(npn.sum())
        m = ~(nf | npn)
        f = F[m]; p = P[m]
        d = np.abs(p - f)
        st["n"] += f.size
        if d.size:
            k = int(np.argmax(np.where(m, np.abs(P - F), -1)))
            if d.max() > st["maxabs"]:
                st["maxabs"] = float(d.max())
                r, c = divmod(k, F.shape[1])
                st["maxabs_loc"] = (row0 + r + 1, c + 1, float(F.flat[k]), float(P.flat[k]))
        st["sumabs"] += float(d.sum())
        st["sf"] += float(f.sum()); st["sp"] += float(p.sum())
        st["sff"] += float((f * f).sum()); st["spp"] += float((p * p).sum()); st["sfp"] += float((f * p).sum())
        st["sse"] += float(((p - f) ** 2).sum())
        st["nexact"] += int((d == 0).sum())
        st["nbig"] = st.get("nbig", 0) + int((d > 0.01 * np.maximum(np.abs(f), 1.0)).sum())
        fl = floor if floor is not None else 1e-3
        rel = d / np.maximum(np.abs(f), fl)
        if rel.size:
            # amostra para quantis
            if rel.size > 20000:
                rel = rel[np.random.default_rng(0).choice(rel.size, 20000, replace=False)]
            st["rel"].append(rel)
        row0 += F.shape[0]
    n = st["n"]
    if n == 0:
        return dict(error="sem dados")
    mf = st["sf"] / n; mp = st["sp"] / n
    vf = st["sff"] / n - mf * mf; vp = st["spp"] / n - mp * mp
    cov = st["sfp"] / n - mf * mp
    corr = cov / np.sqrt(vf * vp) if vf > 0 and vp > 0 else float("nan")
    nse = 1 - st["sse"] / (vf * n) if vf > 0 else float("nan")
    rel = np.concatenate(st["rel"]) if st["rel"] else np.array([0.0])
    return dict(rows=st["rows"], cols=st["ncols_f"], n=n, nan_fortran=st["nnan_f"], nan_python=st["nnan_p"],
                max_abs=st["maxabs"], max_abs_at=st["maxabs_loc"], mean_abs=st["sumabs"] / n,
                rel_median=float(np.median(rel)), rel_p99=float(np.quantile(rel, 0.99)),
                rel_max_sample=float(rel.max()),
                frac_identical=st["nexact"] / n,
                frac_err_gt_1pct=st.get("nbig", 0) / n, corr=float(corr), nse=float(nse),
                sum_bias_pct=100 * (st["sp"] - st["sf"]) / st["sf"] if st["sf"] != 0 else float("nan"),
                mean_fortran=mf, mean_python=mp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fortran_dirs", nargs="+")
    ap.add_argument("--py", required=True)
    ap.add_argument("--nrows", type=int, default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    res = {}
    for name, header, skip in FILES:
        if a.only and name not in a.only.split(","):
            continue
        pf = find(a.fortran_dirs, name)
        pp = os.path.join(a.py, name)
        if pf is None or not os.path.exists(pp):
            res[name] = dict(error="arquivo ausente (fortran=%s, python=%s)" % (pf is not None, os.path.exists(pp)))
            print(f"{name:22s} {res[name]['error']}")
            continue
        r = compare_file(pf, pp, header, skip, a.nrows)
        res[name] = r
        if "error" in r:
            print(f"{name:22s} ERRO: {r['error']}")
        else:
            print(f"{name:22s} n={r['n']:>10d} ident={r['frac_identical']:.4f} maxabs={r['max_abs']:.4g} "
                  f"relmed={r['rel_median']:.2e} relp99={r['rel_p99']:.2e} corr={r['corr']:.8f} "
                  f"gt1%={r['frac_err_gt_1pct']:.2e} bias={r['sum_bias_pct']:+.4f}% nanF/P={r['nan_fortran']}/{r['nan_python']} at={r['max_abs_at']}",
                  flush=True)
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1, default=float)


if __name__ == "__main__":
    main()
