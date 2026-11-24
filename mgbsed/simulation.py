"""
Laco de tempo e saidas da simulacao.

Traducao de: MODELO.f90 e SIMULA.f90 (modo ICALIB=0).
"""
from __future__ import annotations

import os
import time

import numpy as np

from . import hydrology as hy
from . import routing as rt
from . import sediment as sed
from .fortran_io import caldat, fmt_A, fmt_I, fmt_F, fmt_ES, format_F_matrix
from .objective import fobj, sed_fob_calib


# ---------------------------------------------------------------------------
# Escritor de arquivos grandes (uma linha por passo de tempo)
# ---------------------------------------------------------------------------
class RowWriter:
    """Grava linhas '(nI10, mF15.3)' com bufferizacao."""

    def __init__(self, path, header, nint, w=15, d=3, chunk=200, enabled=True):
        self.enabled = enabled
        self.w, self.d = w, d
        self.nint = nint
        self.chunk = chunk
        self.ints = []
        self.vals = []
        if enabled:
            self.f = open(path, "w", newline="\r\n", encoding="latin-1")
            if header is not None:
                self.f.write(header + "\n")

    def add(self, ints, vals):
        if not self.enabled:
            return
        self.ints.append(ints)
        self.vals.append(np.asarray(vals, dtype=np.float64))
        if len(self.vals) >= self.chunk:
            self.flush()

    def flush(self):
        if not self.enabled or not self.vals:
            return
        rows = format_F_matrix(np.vstack(self.vals), self.w, self.d)
        out = []
        for ints, r in zip(self.ints, rows):
            out.append("".join(fmt_I(v, 10) for v in ints) + r)
        self.f.write("\n".join(out) + "\n")
        self.ints, self.vals = [], []

    def close(self):
        if self.enabled:
            self.flush()
            self.f.close()


def _hdr(labels, codes, wl=10, wc=15):
    return "".join(fmt_A(s, wl) for s in labels) + "".join(fmt_I(c, wc) for c in codes)


def _open_text(M, name):
    return open(os.path.join(M.OUTPUT_DIRECTORY, name), "w", newline="\r\n", encoding="latin-1")


# ---------------------------------------------------------------------------
# MODELO - laco do tempo
# ---------------------------------------------------------------------------
def modelo(M, writers=None, progress=True):
    R = M.R
    NT, NC = M.NT, M.NC
    M.IT = 0
    INmini = np.zeros(3); OUTmini = np.zeros(3)
    BALAreia = BALSilte = BALArgila = 0.0
    M.log("SIMULANDO")
    t0 = time.time()
    W = writers or {}
    rain = M.rainfile
    filsol = W.get("NOSOLO")
    filsol2 = W.get("fort25001")
    while M.IT < NT:
        M.IT += 1
        IT = M.IT
        it = IT - 1
        if progress and IT % 100 == 0:
            el = time.time() - t0
            print(f"{100 * IT // NT:4d}%  IT={IT}  {el / 60:.1f} min  (estimado total {el / IT * NT / 60:.1f} min)",
                  flush=True)
        M.JDIA = M.IDINI + int((IT + M.HORAINI - 1) / (86400. / M.DTP))
        M.IMES, M.IDIA, M.IANO = caldat(M.JDIA)
        M.DIAH[it] = M.IDIA; M.MESH[it] = M.IMES; M.ANOH[it] = M.IANO
        M.HORAH[it] = IT % 24 + M.HORAINI
        M.ITCHUVA = IT
        M.P[:] = rain.read(M.ITCHUVA)                   # LECHUVA
        M.PM[:] = M.PM + M.P
        M.TONTEM[:] = M.TA
        hy.leclima(M)
        hy.celula(M)
        if M.SED_INDEX == 1:
            rows = sed.sed_bacia(M)
            if M.ICALIB == 0:
                ints = (IT, M.IDIA, M.IMES, M.IANO)
                W["MUSLE_TOTAL"].add(ints, rows["PSC1"])
                W["MUSLE_areia"].add(ints, rows["PSC2"])
                W["MUSLE_silte"].add(ints, rows["PSC3"])
                W["MUSLE_argila"].add(ints, rows["PSC4"])
                W["QSUP"].add(ints, rows["QSUP"])
                W["QCEL"].add(ints, rows["QCEL"])
        if M.ICALIB == 0 and filsol is not None:
            for f, jc, jb in ((filsol, 0, 0), (filsol2, 1, 1)):
                f.write(fmt_I(IT, 6) + "".join(fmt_F(v, 10, 4) for v in (
                    M.P[jc], M.W[jc, jb], M.SI[jc, jb], M.ET[jc, jb], M.CAF[jc, jb],
                    M.QBAS[jc], M.QINT[jc], M.QSUP[jc])) + "\n")
        if M.hdFLAG0 == 0:
            rt.rede(M)
        M.QCEL2[:] = M.QBAS + M.QINT + M.QSUP
        if M.QITUDO is not None:
            M.QITUDO[:, it] = M.QCEL2
        if M.QBASTUDO is not None:
            M.QBASTUDO[:, it] = M.QBAS
        if M.hdFLAG0 > 0:
            from .inertial import flood_inercial
            flood_inercial(M)
        if M.SED_INDEX == 1:
            srows = sed.sed_rede(M)
            if M.ICALIB == 0:
                ints = (M.IDIA, M.IMES, M.IANO)
                W["CONC_areia"].add(ints, srows["CA"])
                W["CONC_silte"].add(ints, srows["CS"])
                W["CONC_argila"].add(ints, srows["CC"])
                W["DEPdin_RIO"].add(ints, srows["DEP"])
        M.QR[:, it] = M.QJ2[M.IQOBS - 1]
        if M.ICALIB == 0:
            M.QRB[:, it] = M.QJ2[M.IEXUT - 1]
            M.QRG[:, it] = M.QJ2[M.IHIDG - 1]
            if M.hdFLAG0 > 0:
                M.HRG[:, it] = M.Hfl[M.IHIDG - 1]
            _store_sed(M, it, also_conc=True)
            M.Qmini[:, it] = M.QJ2[M.CONTIC - 1] if M.SED_INDEX == 1 else 0
            MWP = 0
            M.QB[it] = M.QJ2[MWP] * M.PJB2[MWP]
            M.QBI[it] = M.QB[it] + M.QJ2[MWP] * M.PJI2[MWP]
            M.QBIS[it] = M.QBI[it] + M.QJ2[MWP] * M.PJS2[MWP]
        else:
            if M.SED_INDEX == 1:
                _store_sed(M, it, also_conc=False)

    if M.SED_INDEX == 1:
        _write_sedrio(M)
        if M.ICALIB == 0:
            INmini[:] = M.iENTRA.sum(axis=0)
            BALAreia = float((M.CSJ2[:, 0].astype(np.float64) * M.VolTREC2).sum())
            BALSilte = float((M.CSJ2[:, 1].astype(np.float64) * M.VolTREC2).sum())
            BALArgila = float((M.CSJ2[:, 2].astype(np.float64) * M.VolTREC2).sum())
            OUTmini[0] += M.iSAI[:, 0].sum() + BALAreia
            OUTmini[1] += M.iSAI[:, 1].sum() + (M.inFL[:, 1].sum() + M.outFL[:, 1].sum()) + BALSilte
            OUTmini[2] += M.iSAI[:, 2].sum() + (M.inFL[:, 2].sum() + M.outFL[:, 2].sum()) + BALArgila
            M.mass_balance = dict(IN=INmini.copy(), OUT=OUTmini.copy(),
                                  BAL=(BALAreia, BALSilte, BALArgila))
            M.log("VERFICACAO DE CONSERVACAO DE SEDIMENTOS")
            M.log("ENTRADA (Ton.)  = ", INmini)
            M.log("SAIDA   (Ton.)  = ", OUTmini)
            with np.errstate(divide="ignore", invalid="ignore"):
                M.log("ERRO    (%)        = ", 100. * (INmini - OUTmini) / INmini)
            M.log(" BALAreia (tons), BALSilte (tons), BALArgila (tons)")
            M.log(BALAreia, BALSilte, BALArgila)
    M.elapsed_modelo = time.time() - t0


def _store_sed(M, it, also_conc):
    """Grandezas de sedimento nos postos (CSSCAL, SEDstore, CAREIA...)."""
    if M.SED_INDEX != 1:
        return
    R = M.R
    e6 = R(10.) ** R(6.)
    q = M.QJ2
    C = M.CSJ2
    for kss in range(M.AUXSS):
        idx = M.ISSOBS[kss] - 1
        M.CSSCAL[kss][:, it] = _sedvar(M, C[idx], q[idx], e6)
        if also_conc:
            idx = M.ISS[kss] - 1
            M.CAREIA[:, it] = C[idx, 0] * e6
            M.CSILTE[:, it] = C[idx, 1] * e6
            M.CARGIL[:, it] = C[idx, 2] * e6
            M.SEDstore[kss][:, it] = _sedvar(M, C[idx], q[idx], e6)


def _sedvar(M, C, q, e6):
    R = M.R
    if M.FLAG_SEDVAR == 1:
        return (C[:, 1] + C[:, 2]) * e6
    if M.FLAG_SEDVAR == 2:
        return (C[:, 1] + C[:, 2]) * e6 * q * R(0.0864)
    if M.FLAG_SEDVAR == 3:
        return (C[:, 0] + C[:, 1] + C[:, 2]) * e6 * q * R(0.0864)
    return np.zeros(C.shape[0], R)


def _write_sedrio(M):
    n = M.nSEDmini
    with _open_text(M, "SEDRIO_SAIDA.txt") as fs, _open_text(M, "SEDRIO_DEPOSITO.txt") as fd, \
            _open_text(M, "SEDRIO_EROSAO.txt") as fe, _open_text(M, "BAL_PLAN.txt") as fb:
        fs.write(fmt_A("IC", 10) + "".join(fmt_A(s, 16) for s in ("SAI_A", "SAI_S", "SAI_C")) + "\n")
        fd.write(fmt_A("IC", 10) + "".join(fmt_A(s, 16) for s in ("DEP_A", "DEP_S", "DEP_C")) + "\n")
        fe.write(fmt_A("IC", 10) + "".join(fmt_A(s, 16) for s in ("ERO_A", "ERO_S", "ERO_C")) + "\n")
        fb.write(fmt_A("IC", 10) + "".join(fmt_A(s, 16) for s in (
            "IN_FLPsand", "IN_FLPsilt", "IN_FLPclay", "OUT_FLPsand", "OUT_FLPsilt",
            "OUT_FLPclay", "QFL_in", "QFL_out", "DFL")) + "\n")
        for i in range(n):
            code = M.CONTIC[i]
            fs.write(fmt_I(code, 10) + "".join(fmt_ES(v, 16, 7) for v in M.iSAI[i]) + "\n")
            fd.write(fmt_I(code, 10) + "".join(fmt_ES(v, 16, 7) for v in M.DEPT[i]) + "\n")
            fe.write(fmt_I(code, 10) + "".join(fmt_ES(v, 16, 7) for v in M.EROT[i]) + "\n")
            vals = list(M.inFL[i]) + list(M.outFL[i]) + [M.BalQFL[i, 0], M.BalQFL[i, 1],
                                                          np.float32(M.DFL[i, 1] + M.DFL[i, 2])]
            fb.write(fmt_I(code, 10) + "".join(fmt_ES(v, 16, 7) for v in vals) + "\n")


# ---------------------------------------------------------------------------
# SIMULA - simulacao padrao (ICALIB = 0)
# ---------------------------------------------------------------------------
def simula(M, progress=True):
    R = M.R
    M.log("Calling initial conditions...")
    hy.condinic(M)
    if M.SED_INDEX == 1:
        sed.sed_inicial(M)
    big = M.write_big_outputs
    out = M.OUTPUT_DIRECTORY
    W = {}
    if M.ICALIB == 0:
        if M.SED_INDEX == 1:
            codes = M.CONTIC
            h4 = _hdr(("IT", "DIA", "MES", "ANO"), codes)
            for nm in ("MUSLE_TOTAL", "MUSLE_areia", "MUSLE_silte", "MUSLE_argila", "QSUP", "QCEL"):
                W[nm] = RowWriter(os.path.join(out, nm + ".txt"), h4, 4, enabled=big)
            W["DEPdin_RIO"] = RowWriter(os.path.join(out, "DEPdin_RIO.txt"), None, 3, enabled=big)
            h3 = _hdr(("DIA", "MES", "ANO"), codes)
            for nm in ("CONC_areia", "CONC_silte", "CONC_argila"):
                W[nm] = RowWriter(os.path.join(out, nm + ".txt"), h3, 3, enabled=big)
    W["NOSOLO"] = _open_text(M, "NOSOLO.txt")
    W["fort25001"] = _open_text(M, "fort.25001")      # unidade 25001 nao aberta no Fortran
    if M.SED_INDEX == 1:
        M.CSSCAL = [M.zeros(int(M.NOBSS[k]), M.NT) for k in range(M.AUXSS)]
        nss = int(M.NUMSS[M.AUXCALIBSS - 1])
        M.SEDstore = [M.zeros(int(M.NUMSS[k]), M.NT) for k in range(M.AUXSS)]
        M.CAREIA = M.zeros(nss, M.NT); M.CSILTE = M.zeros(nss, M.NT); M.CARGIL = M.zeros(nss, M.NT)

    modelo(M, W, progress=progress)

    for k, w in W.items():
        w.close()
    if M.ICALIB == 0 and M.SED_INDEX == 1:
        nss = int(M.NUMSS[M.AUXCALIBSS - 1])
        names = ["MINI_areia.txt", "MINI_silte.txt", "MINI_argila.txt"]
        if M.FLAG_SEDVAR == 1:
            names.append("MINI_CSS.txt")
        elif M.FLAG_SEDVAR == 2:
            names.append("MINI_QSS.txt")
        elif M.FLAG_SEDVAR == 3:
            names.append("MINI_QStot.txt")
        data = [M.CAREIA, M.CSILTE, M.CARGIL, M.SEDstore[0]]
        hdr = fmt_A("IT", 10) + "".join(fmt_I(c, 15) for c in M.ISS[0][:nss])
        for nm, arr in zip(names, data):
            with _open_text(M, nm) as f:
                f.write(hdr + "\n")
                rows = format_F_matrix(arr[:nss].T, 15, 3)
                f.write("\n".join(fmt_I(i + 1, 10) + r for i, r in enumerate(rows)) + "\n")

    fobj(M)
    faju = _open_text(M, "AJUSTE.fob")
    fajus = _open_text(M, "AJUSTESED.fob")
    nh = M.NUMHIDG

    def wrap(vals, per, w, d):
        vals = list(vals)
        lines = []
        for i in range(0, max(len(vals), 1), per):
            lines.append("".join(fmt_F(v, w, d) for v in vals[i:i + per]))
        return "\n".join(lines) + "\n"

    if M.SED_INDEX == 1:
        sed_fob_calib(M)
        if M.ICALIB == 0:
            for L in range(M.NOBS):
                faju.write(wrap([M.R2[L], M.R2L[L], M.ERRV[L]], nh, 10, 3))
            for kss in range(M.AUXSS):
                for K in range(int(M.NOBSS[kss])):
                    fajus.write(wrap([M.SSKGE[K], M.R2SS[K], M.ERRVS[K]], int(M.NOBSS[kss]), 10, 3))
    faju.write(wrap(M.R2[:M.NOBS], nh, 10, 3))
    faju.write(wrap(M.R2L[:M.NOBS], nh, 10, 3))
    faju.write(wrap(M.ERRV[:M.NOBS], nh, 10, 3))
    faju.close(); fajus.close()
    # Adicional (nao existe no Fortran): metricas hidrologicas antes de serem
    # sobrescritas por SED_Fob_Calib. Mesmo formato de AJUSTE.fob.
    with _open_text(M, "AJUSTE_HIDRO_PY.fob") as f:
        f.write("! Gerado pela versao Python: R2 (NSE), R2L (NSE log), ERRV (%) por posto (FOBJ),\n")
        f.write("! valores que no Fortran sao sobrescritos por SED_Fob_Calib antes da gravacao.\n")
        for L in range(M.NOBS):
            f.write(fmt_I(M.IQOBS[L], 10) + "".join(fmt_F(v, 10, 3) for v in (
                M.R2_hidro[L], M.R2L_hidro[L], M.ERRV_hidro[L])) + "\n")

    with _open_text(M, "BACIAS.HIG") as fb, _open_text(M, "QPROP.HIG") as fp, \
            _open_text(M, "VAZAO.HIG") as fh:
        fh.write(fmt_A("IT", 10) + "".join(fmt_I(c, 15) for c in M.IHIDG) + "\n")
        rows = format_F_matrix(M.QRG.T, 15, 3)
        fh.write("\n".join(fmt_I(i + 1, 10) + r for i, r in enumerate(rows)) + "\n")
        ib = 0
        bac = np.vstack([M.WBM[ib], M.PM2[ib], M.SIM[ib], M.EM[ib], M.DSUM[ib], M.DINM[ib], M.DBAM[ib]]).T
        fb.write("\n".join(format_F_matrix(bac, 10, 5)) + "\n")
        prp = np.vstack([M.QB, M.QBI, M.QBIS]).T
        rows = format_F_matrix(prp, 15, 3)
        fp.write("\n".join(fmt_I(i + 1, 10) + r for i, r in enumerate(rows)) + "\n")

    M.PMB[:] = 0.0
    M.PM[:] = (M.PM / R(M.NT)) * (R(365.) * R(24.) * R(3600.) / R(M.DTP))
    kpm_int = np.zeros(M.NB, dtype=np.int64)
    for ic in range(M.NC):
        ibb = M.IBAC[ic]
        M.PMB[ibb] = M.PMB[ibb] + M.PM[ic] * M.ACEL[ic]
        kpm_int[ibb] = int(kpm_int[ibb] + M.ACEL[ic])      # KPM e INTEGER: truncamento a cada soma
    for ibb in range(M.NB):
        M.PMB[ibb] = M.PMB[ibb] / R(kpm_int[ibb])
        M.log("Sub-basin ", ibb + 1, "   Annual Rainfall", f"{M.PMB[ibb]:.3f}")

    if M.SED_INDEX == 1:
        with _open_text(M, "VAZAO_MINI.HIG") as f:
            f.write(fmt_A("IT", 10) + "".join(fmt_I(c, 15) for c in M.CONTIC) + "\n")
            if big:
                for i0 in range(0, M.NT, 500):
                    blk = M.Qmini[:, i0:i0 + 500].T
                    rows = format_F_matrix(blk, 15, 3)
                    f.write("\n".join(fmt_I(i0 + i + 1, 10) + r for i, r in enumerate(rows)) + "\n")
