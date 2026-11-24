"""
Propagacao inercial local (rios + planicies).

Traducao de: flood_TOPO.f90, flood_TAB.f90, flood_ini.f90, flood_inercial.f90,
flood_timestep.f90, flood_discharge.f90, flood_continuity.f90, flood_write.f90,
hunt.f90, FINT.f90, e do ramo inercial de SED_REDE.F90.

ATENCAO: este caminho (IND_FLOODROUTING=1) nao e usado na configuracao do rio
Ivai fornecida (INERTIAL=0) e, portanto, NAO foi verificado contra saidas do
Fortran. A traducao segue o codigo linha a linha; os lacos por minibacia de
flood_discharge e flood_continuity nao tem dependencia entre minibacias dentro
do mesmo sub-passo e foram vetorizados.
"""
from __future__ import annotations

import os

import numpy as np


def flood_topo(M):
    """MINIMONT(J,1) = n de minibacias de montante; MINIMONT(J,2:) = codigos (1-based)."""
    M.MINIMONT[:] = 0
    for ic in range(M.NC):
        j = M.CELJUS[ic]
        if j < 0:
            continue
        j0 = j - 1
        M.MINIMONT[j0, 0] += 1
        k = M.MINIMONT[j0, 0]
        if 1 + k >= M.MINIMONT.shape[1]:
            raise RuntimeError("MINIMONT: mais de 9 minibacias de montante (limite do Fortran)")
        M.MINIMONT[j0, k] = ic + 1


def flood_tab(M):
    for ic in range(M.NC):
        M.ZTAB[0, ic] = M.ZFUNDOFL[ic] - M.HRIO[ic]
        M.VTAB[0, ic] = 0.0
        M.ATAB[0, ic] = 0.0
        M.ZTAB[1, ic] = M.ZFUNDOFL[ic]
        M.VTAB[1, ic] = float(M.BRIO[ic]) * M.HRIO[ic] * float(M.SRIO[ic]) * 1000.0
        M.ATAB[1, ic] = 0.0
        M.ZTAB[2, ic] = M.ZFL[0, ic]
        M.VTAB[2, ic] = M.VTAB[1, ic] + 0.5 * (M.AFL[0, ic] + 0.0) * (M.ZTAB[2, ic] - M.ZTAB[1, ic]) * 1000000.0
        M.ATAB[2, ic] = M.AFL[0, ic]
        for k in range(M.NPFL[ic] - 1):
            M.ZTAB[k + 3, ic] = M.ZFL[k + 1, ic]
            VOLINCR = np.float32((M.ZFL[k + 1, ic] - M.ZFL[k, ic]) * 1000000.0 *
                                 ((M.AFL[k, ic] + M.AFL[k + 1, ic]) * 0.5))      # REAL (sic)
            M.VTAB[k + 3, ic] = M.VTAB[k + 2, ic] + float(VOLINCR)
            M.ATAB[k + 3, ic] = M.AFL[k + 1, ic]


def fint(X, Y, N, ABC):
    """FINT.f90 - interpolacao linear com extrapolacao no ultimo intervalo."""
    for i in range(1, N):
        if ABC < X[i]:
            j = i - 1
            return Y[j] + (Y[i] - Y[j]) * (ABC - X[j]) / (X[i] - X[j])
        elif ABC == X[i]:
            return Y[i]
        elif i == N - 1:
            j = i - 1
            return Y[j] + (Y[i] - Y[j]) * (ABC - X[j]) / (X[i] - X[j])
    return Y[0]


def flood_ini(M):
    M.dtflood0 = float(M.DTP)
    for ic in range(M.NC):
        ALFA = (M.nMan[ic] * (float(M.BRIO[ic]) ** (2.0 / 3.0))) ** 0.6 / ((float(M.DECL[ic]) ** 0.5) ** 0.6)
        BETA = 0.6
        M.Q2fl[ic] = M.QREF[ic]
        M.Hfl[ic] = ALFA * (M.Q2fl[ic] ** BETA) / float(M.BRIO[ic]) + 0.00001
        M.Yfl[ic] = M.ZTAB[0, ic] + M.Hfl[ic]
        M.Vol1[ic] = M.Hfl[ic] * float(M.BRIO[ic]) * float(M.SRIO[ic] * np.float32(1000.))
        M.Vol2[ic] = M.Vol1[ic]
        n = M.NPFL[ic] + 2
        M.Area2[ic] = fint(M.VTAB[:, ic], M.ATAB[:, ic], n, M.Vol2[ic])
        M.Area2[ic] = min(M.Area2[ic], float(M.ACEL[ic]))
        M.jtab[ic] = 1
    # (sic) Q2face=QREF(IC) atribui o vetor inteiro com o ultimo QREF
    M.Q2face = np.full(M.NC, float(M.QREF[M.NC - 1]))


def flood_timestep(M):
    hmax = np.maximum(M.Hfl, 0.001)
    dt = M.alpha * M.SRIO.astype(np.float64) * 1000. / (M.g * hmax) ** 0.5
    dtteste = min(float(M.DTP), float(dt.min()))
    return min(float(M.DTP), dtteste)


def flood_discharge(M, dtflood):
    g = M.g
    z1 = M.ZTAB[0]
    y1 = M.Hfl + z1
    jus = M.jus
    has = jus >= 0
    jj = np.where(has, jus, 0)
    z2_out = np.minimum(z1, z1 - M.DECL.astype(np.float64) * M.SRIO.astype(np.float64) * 1000)
    z2 = np.where(has, M.ZTAB[0, jj], z2_out)
    y2 = np.where(has, M.Hfl[jj] + z2, M.Hfl + z2)
    hflow = np.maximum(np.maximum(y2, y1) - np.maximum(z2, z1), 0.0)
    sr = (M.SRIO * np.float32(1000.)).astype(np.float64)
    dxflow = np.where(has, (sr + sr[jj]) / 2., sr)
    bflow = M.BRIO.astype(np.float64)
    xMan = M.nMan
    q0 = M.Q2fl / bflow
    Sflow = -(y1 - y2) / dxflow
    with np.errstate(divide="ignore", invalid="ignore"):
        q = (q0 - (g * dtflood * hflow * Sflow)) / (1. + g * dtflood * hflow * xMan * xMan * np.abs(q0) /
                                                   (hflow ** (10.0 / 3.0)))
        q = q * bflow
        lim = bflow * (hflow ** 1.5) * 9.81 ** 0.5
    q = np.where(q > 0, np.minimum(q, lim), np.maximum(q, -lim))
    q = np.where(hflow > 0.0, q, 0.0)
    M.Q2fl[:] = q
    M.QJ2[:] = q
    if M.ICALIB == 3 and getattr(M, "Q_SC", None) is not None:
        M.Q_SC[:, M.IT - 1] = M.QJ2


def _hunt_vec(M, x):
    """hunt.f90 vetorizado. Retorna (Area, Cota) interpoladas em VTAB.

    Resultado de hunt = maior j (1-based) com VTAB(j) <= x, com os ajustes de
    igualdade nas extremidades. Busca "quente" a partir de jtab (como o hunt
    original): onde o intervalo anterior ainda contem x, nao ha nova busca.
    """
    NC = M.NC
    if getattr(M, "_hunt_cache", None) is None:
        n = M.NPFL + 2
        V = M.VTAB
        big = np.arange(V.shape[0])[:, None] >= n[None, :]
        Vp = np.where(big, np.inf, V).T.copy()             # (NC, 102)
        Vp_ext = np.concatenate([np.full((NC, 1), -np.inf), Vp, np.full((NC, 1), np.inf)], axis=1)
        M._hunt_cache = (n, Vp, Vp_ext, np.arange(NC))
    n, Vp, Vp_ext, rows = M._hunt_cache
    V = M.VTAB
    j = M.jtab.copy()
    # contagem de VTAB <= x : valido se Vp_ext[j] <= x < Vp_ext[j+1]
    ok = (Vp_ext[rows, j] <= x) & (x < Vp_ext[rows, np.minimum(j + 1, Vp_ext.shape[1] - 1)])
    if not ok.all():
        b = ~ok
        j[b] = (Vp[b] <= x[b, None]).sum(axis=1)
    cnt = j.copy()
    jlo = np.where(x == V[n - 1, rows], n - 1, cnt)
    jlo = np.where(x == V[0], 1, jlo)
    out = []
    j0 = np.clip(jlo - 1, 0, None)
    j1 = np.clip(jlo, 0, V.shape[0] - 1)
    j5 = np.clip(jlo - 1 - 5, 0, None)
    v0 = V[j0, rows]; v1 = V[j1, rows]
    with np.errstate(divide="ignore", invalid="ignore"):
        inc = jlo + (1 / (v1 - v0)) * (x - v0)
        for Y in (M.ATAB, M.ZTAB):
            y0 = Y[j0, rows]
            y_in = y0 + (Y[j1, rows] - y0) * (inc - jlo)
            y_ex = y0 + (y0 - Y[j5, rows]) / 5 * jlo
            out.append(np.where(jlo < n, y_in, y_ex))
    M.jtab[:] = cnt          # guarda a contagem (ponto de partida da proxima busca)
    return out[0], out[1]


def flood_continuity(M, dtflood, tflood):
    NC = M.NC
    # soma das vazoes de montante (MINIMONT), na ordem das colunas
    up = np.zeros(NC)
    for k in range(1, M.MINIMONT.shape[1]):
        m = M.MINIMONT[:, 0] >= k
        if not m.any():
            break
        up[m] += M.Q2fl[M.MINIMONT[m, k] - 1]
    E0 = M.E0agua.astype(np.float64); P = M.P.astype(np.float64)
    DTP = float(M.DTP)
    Vol2 = M.Vol1 + dtflood * (up + M.QCEL2.astype(np.float64) - M.Q2fl) \
        - (E0 * dtflood * M.Area2 * 1000000.0 / (DTP * 1000.)) \
        + (P * dtflood * M.Area2 * 1000000.0 / (DTP * 1000.))
    Vol2 = np.maximum(Vol2, 0.0)
    M.Vol2[:] = Vol2
    area, yj = _hunt_vec(M, Vol2)
    M.Area2[:] = np.minimum(area, M.ACEL.astype(np.float64))
    y2 = np.maximum(yj, M.ZTAB[0]) + 0.001
    M.Hfl[:] = y2 - M.ZTAB[0]
    M.Yfl[:] = y2
    M.Vol1[:] = M.Vol2
    if tflood / M.dtflood0 == 1:
        M.AFLTUDO[0, M.IT - 1] += np.float32(M.Area2.sum())


def flood_inercial(M):
    tflood = 0.0
    M.AFLTUDO[0, M.IT - 1] = 0.0
    while tflood < M.dtflood0:
        dtflood = flood_timestep(M)
        dtflood = min(dtflood, M.dtflood0 - tflood)
        tflood = tflood + dtflood
        flood_discharge(M, dtflood)
        flood_continuity(M, dtflood, tflood)
    M.QJ2[:] = M.Q2fl
    if M.QTUDO is not None:
        M.QTUDO[:, M.IT - 1] = M.Q2fl
    if M.YTUDO is not None:
        M.YTUDO[M.IT - 1, :] = M.Yfl


def flood_write(M):
    """flood_write.f90: grava Yfl de todas as minibacias em Output/flood/Yfl_#####.txt."""
    d = os.path.join(M.OUTPUT_DIRECTORY, "flood")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"Yfl_{M.IT:05d}.txt"), "w") as f:
        for v in M.Yfl:
            f.write(f"{v:16.6f}\n")


def sed_hydraulics_inertial(M, act, HmJ, Sf, Vm, Uat, QtFL, HFLSD, Qlim):
    """Ramo 'ELSE' (hdFLAG0>0) de SED_REDE: geometria hidraulica do trecho (vetorizado).

    z2/y2/bflow sao variaveis de modulo: na minibacia exutorio (CELJUS=-1) o
    Fortran reaproveita os valores da minibacia anterior processada no laco
    (sic). Isso e reproduzido com "forward fill" na ordem de IC.
    """
    R = M.R
    NC = M.NC
    DTP = R(M.DTP)
    idx = np.nonzero(act)[0]
    if idx.size == 0:
        return
    j = M.jus[idx]
    has = M.CELJUS[idx] != -1
    jj = np.where(has, j, 0)
    z1 = M.ZTAB[0, idx]
    y1 = M.Hfl[idx] + z1
    bf = (M.BRIO[idx].astype(np.float64) + M.BRIO[jj].astype(np.float64)) / 2.
    z2 = M.ZTAB[0, jj]
    y2 = M.Hfl[jj] + z2
    # forward fill para exutorios (valor da minibacia anterior no laco; 0 se nenhuma)
    pos = np.where(has, np.arange(idx.size), -1)
    last = np.maximum.accumulate(pos)
    src = np.maximum(last, 0)
    okp = last >= 0
    bf = np.where(okp, bf[src], 0.0); z2 = np.where(okp, z2[src], 0.0); y2 = np.where(okp, y2[src], 0.0)
    hflow = np.maximum(np.maximum(y2, y1) - np.maximum(z2, z1), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        q0 = M.Q2fl[idx] / bf
        xMan = M.nMan[idx]
        h = hflow.astype(R)
        HmJ[idx] = h
        Sf[idx] = np.maximum((xMan * xMan * np.abs(q0) * q0) / (h.astype(np.float64) ** (10.0 / 3.0)), 0.00001)
        M.VolTREC2[idx] = h * M.BRIO[idx] * M.SRIO[idx] * R(1000.)
        M.VTOTAL[idx, 0] = M.VolTREC2[idx]
        Vm[idx] = np.where(M.QJ2[idx] > Qlim, M.QJ2[idx] / (M.BRIO[idx] * h), R(0.0))
        RH = (h * M.BRIO[idx]) / (R(2.0) * h + M.BRIO[idx])
        Uat[idx] = fm_sqrt(R(9.81) * RH * Sf[idx])
        AtFL = np.maximum(0., (M.Area2[idx] * 1000000.) - (M.BRIO[idx] * M.SRIO[idx] * R(1000.)).astype(np.float64))
        M.VFL2[idx] = np.maximum(0., M.Vol2[idx] - M.VTAB[2, idx])
        QtFL[idx] = (M.VFL2[idx] - M.VFL1[idx]) / DTP
        M.VTOTAL[idx, 1] = M.VFL2[idx]
        HFLSD[idx] = M.VFL2[idx] / AtFL.astype(R)
    if M.ICALIB == 3:
        t = M.IT - 1
        M.HmJx_SC[idx, t] = HmJ[idx]
        M.SfTRECx_SC[idx, t] = Sf[idx]
        M.VOLTREC2_SC[idx, t] = M.VolTREC2[idx]
        M.VFL2_SC[idx, t] = M.VFL2[idx]
        M.HFLSD_SC[idx, t] = HFLSD[idx]


def fm_sqrt(x):
    from .fmath import sqrt
    return sqrt(x)
