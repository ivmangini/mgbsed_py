"""
Funcoes objetivo.

Traducao de: FOBJ.f90 e SED_Fob_Calib.f90.

As somas sao feitas sequencialmente em precisao REAL (np.cumsum), como nos
lacos DO do Fortran, para reproduzir o arredondamento em precisao simples.
"""
from __future__ import annotations

import warnings

import numpy as np


def _seqsum(x, R):
    x = np.asarray(x, dtype=R)
    if x.size == 0:
        return R(0.0)
    return np.cumsum(x, dtype=R)[-1]


def fobj(M):
    R = M.R
    iniaju = 360
    M.R2[:] = -100.0
    M.R2L[:] = -100.0
    M.ERRV[:] = 1000.0
    NT = M.NT
    with np.errstate(divide="ignore", invalid="ignore"):
        for kb in range(M.NOBS):
            if iniaju < NT:
                qo = M.QOBS[kb, iniaju - 1:NT]
                qc = M.QR[kb, iniaju - 1:NT]
                m = qo >= 0.0
                qo = qo[m]; qc = qc[m]
                ITVAL = int(m.sum())
                SOMAOBS = _seqsum(qo, R)
                SOMACAL = _seqsum(qc, R)
                lo = np.log(qo + R(0.0001)); lc = np.log(qc + R(0.0001))
                SOMALOGS = _seqsum(lo, R)
                XMOBS = SOMAOBS / R(ITVAL) if ITVAL else R(np.nan)
                XLOGS = SOMALOGS / R(ITVAL) if ITVAL else R(np.nan)
                SQDES = _seqsum((qo - qc) ** R(2.0), R)
                SOMQ = _seqsum((qo - XMOBS) ** R(2.0), R)
                SQLOG = _seqsum((lo - lc) ** R(2.0), R)
                SOLOG = _seqsum((lo - XLOGS) ** R(2.0), R)
                if ITVAL > 0:
                    M.R2[kb] = R(1.) - SQDES / SOMQ
                    M.R2L[kb] = R(1.) - SQLOG / SOLOG
                    M.ERRV[kb] = (SOMACAL / SOMAOBS - R(1.0)) * R(100.)
    M.IMES = 0
    if M.ICALIB == 0:
        # Fortran: calibFLAG=1.0 e NF nao inicializado (=0) -> POND=0 e VFO=0
        M.CalibFLAG = np.ones((M.NOBS, M.NFP), dtype=R)
        nf = getattr(M, "NF", 0)
    else:
        nf = M.NF
    POND = np.zeros((M.NOBS, 3), dtype=R)
    for kb in range(M.NOBS):
        for j in range(nf):
            if M.CalibFLAG[kb, j] > 0.0:
                POND[kb, j] = M.CalibFLAG[kb, j] / M.CalibFLAG[:M.NOBS, j].sum()
    M.VFO[:] = 0.0
    for kb in range(M.NOBS):
        if POND[kb, 0] > 0.0:
            M.VFO[0] += (R(1.) - M.R2[kb]) * POND[kb, 0]
        if POND[kb, 1] > 0.0:
            M.VFO[1] += (R(1.) - M.R2L[kb]) * POND[kb, 1]
        if POND[kb, 2] > 0.0:
            M.VFO[2] += np.abs(M.ERRV[kb]) * POND[kb, 2]
    # Copia (adicional ao Fortran): SED_Fob_Calib sobrescreve R2/R2L/ERRV.
    M.R2_hidro = M.R2.copy(); M.R2L_hidro = M.R2L.copy(); M.ERRV_hidro = M.ERRV.copy()
    M.log(f"{M.VFO[0]} {M.VFO[1]} {M.VFO[2]}")


def sed_fob_calib(M):
    R = M.R
    NT = M.NT
    nobs_max = M.NMAXHID
    M.SSCORRELtemp = np.zeros(nobs_max, R); M.SSCORRELesp = np.zeros(nobs_max, R)
    M.SSCORRELtudo = np.zeros(nobs_max, R)
    # (sic) estas atribuicoes apagam as metricas hidrologicas calculadas em FOBJ
    M.R2[:] = -100.0
    M.R2L[:] = -100.0
    M.ERRV[:] = 1000.0
    M.R2SS = np.full(nobs_max, -100.0, R)
    M.RMSE = np.zeros(nobs_max, R)
    M.ERRVS = np.zeros(nobs_max, R)
    M.SSKGE = np.full(nobs_max, -100.0, R)
    M.KGEa = np.full(nobs_max, -100.0, R)
    M.CoefB = np.full(nobs_max, -100.0, R)
    M.DCPERM = np.full(nobs_max, 1000.0, R)
    M.SSCORRELAC = np.zeros(nobs_max, R)
    kss = M.AUXCALIBSS - 1
    iniaju = 365
    nobs = int(M.NOBSS[kss])
    SSDES = np.zeros(nobs, R)
    Resp = np.zeros((nobs, 2), R)
    obs_all = []
    cal_all = []
    M.VFO[:] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        for kb in range(nobs):
            so = M.SSOBS[kss][kb, iniaju - 1:NT]
            sc = M.CSSCAL[kss][kb, iniaju - 1:NT]
            m = so >= 0.0
            so = so[m]; sc = sc[m]
            SITVAL = int(m.sum())
            SSOMAOBS = _seqsum(so, R)
            SSOMACAL = _seqsum(sc, R)
            obs_all.append(so); cal_all.append(sc)
            # curva de permanencia (acumulada entre postos, como no Fortran)
            QSOROBS = np.sort(np.concatenate(obs_all))
            QSORCALC = np.sort(np.concatenate(cal_all))
            auxcont = QSOROBS.size
            ipos = (90 * auxcont) // 100
            if ipos < 1 and kb == 0:
                # No Fortran, QSOROBS(0) esta fora do vetor: aborta em Debug
                # (/check:bounds) e le memoria indefinida em Release.
                warnings.warn("SED_Fob_Calib: menos de 2 observacoes de sedimento a partir do "
                              f"dia {iniaju}; curva de permanencia indefinida no Fortran "
                              "(aqui Q10 = Q50 = 0).", RuntimeWarning, stacklevel=2)
            Q10 = (QSOROBS[ipos - 1] if ipos >= 1 else R(0), QSORCALC[ipos - 1] if ipos >= 1 else R(0))
            ipos = (50 * auxcont) // 100
            Q50 = (QSOROBS[ipos - 1] if ipos >= 1 else R(0), QSORCALC[ipos - 1] if ipos >= 1 else R(0))
            M.DCPERM[kb] = R(100) * ((Q50[1] - Q10[1]) / (Q50[0] - Q10[0]) - R(1))
            SSXOBS = SSOMAOBS / R(SITVAL) if SITVAL else R(np.nan)
            CSSXCAL = SSOMACAL / R(SITVAL) if SITVAL else R(np.nan)
            Resp[kb, 0] = SSXOBS; Resp[kb, 1] = CSSXCAL
            SSDPOBS = _seqsum((so - SSXOBS) ** R(2.0), R)
            SSDPCAL = _seqsum((sc - CSSXCAL) ** R(2.0), R)
            SSCOVAR = _seqsum((so - SSXOBS) * (sc - CSSXCAL), R)
            SSDES[kb] = _seqsum((so - sc) ** R(2.0), R)
            SOMS = _seqsum((so - SSXOBS) ** R(2.0), R)
            if SITVAL > 0:
                M.R2SS[kb] = R(1.) - SSDES[kb] / SOMS
                M.ERRVS[kb] = (SSOMACAL / SSOMAOBS - R(1.0)) * R(100.)
                M.RMSE[kb] = (SSDES[kb] / R(SITVAL)) ** R(0.5)
                M.CoefB[kb] = (CSSXCAL - SSXOBS) / ((SSDPOBS / R(SITVAL)) ** R(0.5))
            M.SSCORRELAC[kb] = SSCOVAR / ((SSDPCAL * SSDPOBS) ** R(0.5))
            M.SSCORRELtemp[kss] = M.SSCORRELtemp[kss] + M.SSCORRELAC[kb] / R(nobs)
            r = float(M.SSCORRELAC[kb])
            M.SSKGE[kb] = 1 - ((1 - r) ** 2 + (1 - float(CSSXCAL / SSXOBS)) ** 2 +
                               (1 - float((SSDPCAL / SSDPOBS) ** R(0.5))) ** 2) ** 0.5
            M.KGEa[kb] = float((SSDPCAL / SSDPOBS) ** R(0.5))
    if M.ICALIB == 3:
        for kb in range(nobs):
            M.VFO[0] += (1 - M.R2SS[kb])
            M.VFO[1] += abs(M.ERRVS[kb])
            M.VFO[2] += (1 - M.SSKGE[kb])
    M.log("VFO(1),VFO(2),VFO(3)")
    M.log(f"{M.VFO[0]} {M.VFO[1]} {M.VFO[2]}")
