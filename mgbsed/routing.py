"""
Propagacao na rede de drenagem por Muskingum-Cunge.

Traducao de: REDE.f90 e MUSK.f90.

Estrategia
----------
No Fortran, REDE percorre as minibacias em ordem crescente de indice (todas as
minibacias de montante tem indice menor que o de jusante). Cada minibacia com
rio executa MUSK: uma recorrencia de Muskingum-Cunge em NTRECH sub-trechos x
NTC sub-intervalos, com MAX(.,0) em cada ponto da grade.

Aqui:
1. As minibacias sao agrupadas por "nivel topologico" (todas as de montante de
   uma minibacia estao em niveis menores). Minibacias do mesmo nivel sao
   independentes e sao processadas juntas (vetorizacao).
2. Sem o MAX(.,0), a recorrencia e linear nas entradas (condicao inicial
   QRIOINI e contorno QCONTORM). Para cada minibacia pre-calcula-se a matriz G
   que leva as entradas a grade completa QC. Se a grade linear nao tem valores
   negativos, o MAX(.,0) nunca atua e o resultado e identico ao da recorrencia
   original (a menos de arredondamento). Se houver qualquer negativo, a
   minibacia e recalculada com a recorrencia original (exata), ponto a ponto.
"""
from __future__ import annotations

import numpy as np

from . import fmath as fm


def topo_levels(M):
    """Nivel topologico de cada minibacia (0 = cabeceira)."""
    lev = np.zeros(M.NC, int)
    for i in range(M.NC):
        j = M.jus[i]
        if j >= 0:
            lev[j] = max(lev[j], lev[i] + 1)
    M.level = lev
    M.levels = [np.nonzero(lev == L)[0] for L in range(lev.max() + 1)]


def _musk_coefs(M):
    """Coeficientes C1,C2,C3 de REDE (constantes no tempo) em precisao REAL."""
    R = M.R
    NC = M.NC
    C1 = M.zeros(NC); C2 = M.zeros(NC); C3 = M.zeros(NC)
    for ic in range(NC):
        if M.NSUBT[ic] <= 0:
            continue
        DTCAL = M.DT[ic]
        DX = M.SRIO[ic] * R(1000.) / R(M.NSUBT[ic])
        AK = DX / M.CEL[ic]
        XXX = R(0.5) * (R(1.0) - (M.QREF[ic] / (M.BRIO[ic] * M.DECL[ic] * M.CEL[ic] * DX)))
        if XXX < 0.2:
            XXX = R(0.2)
            AK = DTCAL / (R(3.125) * fm.pw(XXX, R(1.25)))
        elif XXX < 0.4:
            pass
        elif XXX < 0.5:
            pass
        else:
            XXX = R(0.5)
            AK = DTCAL
        DEN = R(2.) * AK * (R(1.) - XXX) + DTCAL
        C1[ic] = (R(2.) * AK * XXX + DTCAL) / DEN
        C2[ic] = (DTCAL - R(2.) * AK * XXX) / DEN
        C3[ic] = (R(2.) * AK * (R(1.) - XXX) - DTCAL) / DEN
    return C1, C2, C3


def prepare_routing(M):
    """Pre-calculos (equivalentes a constantes de REDE/MUSK)."""
    topo_levels(M)
    if M.DT[M.NC - 1] < 0.00001:
        raise RuntimeError("Something wrong in time-step on rede.f90")
    M.NTMUSK = int(M.DTP / float(M.DT[M.NC - 1]))
    M.MC_C1, M.MC_C2, M.MC_C3 = _musk_coefs(M)
    riv = M.NSUBT > 0
    M.NTC = np.zeros(M.NC, int)
    M.NTC[riv] = (M.DTP / M.DT[riv].astype(np.float64)).astype(int)
    if riv.any() and (M.NTC[riv] != M.NTMUSK).any():
        raise NotImplementedError("NTC diferente entre minibacias (nao ocorre com DT=DTP/24)")
    T = M.NTMUSK
    nmax = int(M.NSUBT.max()) if riv.any() else 0
    M.MC_nmax = nmax
    # Matriz linear G (nc, nmax*T, nmax+1+T): grade QC[1..n][1..T] <- [QRIOINI(0..nmax), b(1..T)]
    nin = nmax + 1 + T
    G = np.zeros((M.NC, nmax * T, nin))
    for ic in np.nonzero(riv)[0]:
        n = int(M.NSUBT[ic])
        c1, c2, c3 = (float(M.MC_C1[ic]), float(M.MC_C2[ic]), float(M.MC_C3[ic]))
        # simulacao simbolica: cada valor da grade e um vetor de coeficientes
        QC = np.zeros((n + 1, T + 1, nin))
        for r in range(n + 1):
            QC[r, 0, r] = 1.0                       # QC(r,1) = QRIOINI(r)
        for t in range(1, T + 1):
            QC[0, t, nmax + t] = 1.0                # QC(1,t) = b(t)
        for t in range(T):
            for r in range(n):
                QC[r + 1, t + 1] = c1 * QC[r, t] + c2 * QC[r, t + 1] + c3 * QC[r + 1, t]
        for r in range(1, n + 1):
            G[ic, (r - 1) * T:r * T, :] = QC[r, 1:, :]
    M.MC_G = G


def _musk_exact(M, ic, qcm, qrioini):
    """MUSK.f90 original (recorrencia com MAX(.,0)), em precisao REAL."""
    R = M.R
    n = int(M.NSUBT[ic]); T = int(M.NTMUSK)
    c1, c2, c3 = M.MC_C1[ic], M.MC_C2[ic], M.MC_C3[ic]
    QC = np.zeros((n + 1, T + 1), dtype=R)
    QPERDAS = R(0.0)
    for itc in range(T + 1):
        QC[0, itc] = max(qcm[itc] - M.EVQ[ic] - QPERDAS, R(0.0))
    for itr in range(n + 1):
        QC[itr, 0] = qrioini[itr]
    contorj = np.zeros(T + 1, dtype=R)
    for itc in range(T):
        for itr in range(n):
            v = c1 * QC[itr, itc] + c2 * QC[itr, itc + 1] + c3 * QC[itr + 1, itc]
            QC[itr + 1, itc + 1] = max(v, R(0.0))
        contorj[itc + 1] = QC[n, itc + 1]
    return QC, contorj


def rede(M):
    R = M.R
    NC = M.NC
    DTP = R(M.DTP)
    T = M.NTMUSK
    it = M.IT - 1
    nmax = M.MC_nmax
    jus = M.jus

    M.QM1[:NC] = M.QM2[:NC]
    M.QJ1[:] = M.QJ2
    M.QM2[:] = 0.0
    M.QJ2[:] = 0.0
    M.QCEL1[:] = M.QCEL2
    M.QCONTORM[:] = 0.0
    M.QCONTORJ[:] = 0.0

    act = (M.sbtFLAG != 1)
    M.QCEL2[act] = (M.QBAS + M.QINT + M.QSUP)[act]
    compute = act & (M.HDFLAG != 1)

    lt = np.arange(T + 1, dtype=R)                 # (LT-1)
    subst = {int(c): i for i, c in enumerate(M.ISUBST)} if M.NUMSUBST > 0 else {}

    for cells in M.levels:
        cells = cells[compute[cells]]
        if cells.size == 0:
            continue
        riv = cells[M.NSUBT[cells] > 0]
        nor = cells[M.NSUBT[cells] <= 0]

        # ---------------- minibacias com rio ---------------------------
        if riv.size:
            qbas, qint, qsup = M.QBAS[riv], M.QINT[riv], M.QSUP[riv]
            qm_up = M.QM2[riv]
            M.PMB2[riv] = M.PMB2[riv] * qm_up
            M.PMI2[riv] = M.PMI2[riv] * qm_up
            M.PMS2[riv] = M.PMS2[riv] * qm_up
            M.QM2[riv] = ((qm_up + qsup) + qint) + qbas
            qc1 = M.QCEL1[riv]; qc2 = M.QCEL2[riv]
            M.QCONTORM[riv, :T + 1] = M.QCONTORM[riv, :T + 1] + (
                qc1[:, None] + lt[None, :] * (qc2 - qc1)[:, None] / R(T))
            qm2 = M.QM2[riv]
            M.PMB2[riv] = (M.PMB2[riv] + qbas) / qm2
            M.PMI2[riv] = (M.PMI2[riv] + qint) / qm2
            M.PMS2[riv] = (M.PMS2[riv] + qsup) / qm2
            M.VRB[riv] = M.VRB[riv] + M.PMB2[riv] * qm2 * DTP
            M.VRI[riv] = M.VRI[riv] + M.PMI2[riv] * qm2 * DTP
            M.VRS[riv] = M.VRS[riv] + M.PMS2[riv] * qm2 * DTP

            # ---- MUSK ----
            qcm = M.QCONTORM[riv, :T + 1]
            b = np.maximum(qcm[:, 1:] - M.EVQ[riv][:, None], R(0.0))   # QC(1,2..T+1)
            x = np.zeros((riv.size, nmax + 1 + T))
            x[:, :nmax + 1] = M.QRIOINI[riv, :nmax + 1]
            x[:, nmax + 1:] = b
            grid = np.einsum("kij,kj->ki", M.MC_G[riv], x).reshape(riv.size, nmax, T)
            nn = M.NSUBT[riv]
            rows = np.arange(riv.size)
            # zera linhas alem de NTRECH (padding) para o teste de negatividade
            rmask = (np.arange(nmax)[None, :] < nn[:, None])
            neg = ((grid < 0.0) & rmask[:, :, None]).any(axis=(1, 2))
            contorj = grid[rows, nn - 1, :].astype(R)                # QC(NTRECH+1, 2..T+1)
            new_rio = grid[:, :, T - 1].astype(R)                       # QC(2..NTRECH+1, T+1)
            M.QCONTORJ[riv, 0] = M.QJ1[riv]
            M.QCONTORJ[riv, 1:T + 1] = contorj
            qrio_old = M.QRIOINI[riv, :nmax + 1].copy()
            M.QRIOINI[riv, 0] = qcm[:, T]
            M.QRIOINI[riv, 1:nmax + 1] = np.where(rmask, new_rio, M.QRIOINI[riv, 1:nmax + 1])
            M.QJ2[riv] = contorj[:, T - 1]
            if neg.any():
                for k in np.nonzero(neg)[0]:
                    ic = riv[k]
                    n = int(nn[k])
                    QC, cj = _musk_exact(M, ic, qcm[k], qrio_old[k])
                    M.QCONTORJ[ic, 1:T + 1] = cj[1:]
                    M.QRIOINI[ic, 0] = qcm[k, T]
                    M.QRIOINI[ic, 1:n + 1] = QC[1:n + 1, T]
                    M.QJ2[ic] = QC[n, T]
                M.n_musk_exact = getattr(M, "n_musk_exact", 0) + int(neg.sum())

            qj2 = M.QJ2[riv]
            VT = M.VRB[riv] + M.VRI[riv] + M.VRS[riv]
            M.PJB2[riv] = M.VRB[riv] / VT
            M.PJI2[riv] = M.VRI[riv] / VT
            M.PJS2[riv] = M.VRS[riv] / VT
            M.VRB[riv] = np.maximum(R(0.0), M.VRB[riv] - M.PJB2[riv] * qj2 * DTP)
            M.VRI[riv] = np.maximum(R(0.0), M.VRI[riv] - M.PJI2[riv] * qj2 * DTP)
            M.VRS[riv] = np.maximum(R(0.0), M.VRS[riv] - M.PJS2[riv] * qj2 * DTP)

            # contribuicao para jusante
            j = jus[riv]
            has = j >= 0
            if has.any():
                src = riv[has]; jj = j[has]
                q = M.QJ2[src]
                cb = M.PJB2[src] * q; ci = M.PJI2[src] * q; cs = M.PJS2[src] * q
                bnd = M.QCONTORJ[src, :T + 1].copy()
                if subst:
                    for kk, ic in enumerate(src):
                        if M.IT > 1 and (ic + 1) in subst:
                            isb = subst[ic + 1]
                            if M.QLIDO[isb, it] >= 0.0:
                                # vazao lida substitui a calculada (QM2(JUS)=QLIDO)
                                bnd[kk] = (M.QLIDO[isb, it - 1] + lt * (M.QLIDO[isb, it] - M.QLIDO[isb, it - 1]) / R(T))
                _accumulate_downstream(M, jj, q, cb, ci, cs, bnd, src, subst, it)

        # ---------------- minibacias sem rio ---------------------------
        if nor.size:
            qbas, qint, qsup = M.QBAS[nor], M.QINT[nor], M.QSUP[nor]
            q = ((qbas + qint) + qsup)
            M.QJ2[nor] = q
            with np.errstate(invalid="ignore", divide="ignore"):
                M.PJB2[nor] = qbas / q
                M.PJI2[nor] = qint / q
                M.PJS2[nor] = qsup / q
            j = jus[nor]
            has = j >= 0
            if has.any():
                src = nor[has]; jj = j[has]
                qq = ((M.QSUP[src] + M.QINT[src]) + M.QBAS[src])
                bnd = (M.QJ1[src][:, None] + lt[None, :] * (M.QJ2[src] - M.QJ1[src])[:, None] / R(T))
                _accumulate_downstream(M, jj, qq, M.QBAS[src], M.QINT[src], M.QSUP[src], bnd,
                                       None, {}, it)

        if M.QTUDO is not None:
            M.QTUDO[cells, it] = M.QJ2[cells]


def _accumulate_downstream(M, jj, q, cb, ci, cs, bnd, src, subst, it):
    """Equivalente as linhas de REDE que somam a contribuicao de IC em JUS:
        PMx2(JUS) = (PMx2(JUS)*QM2(JUS) + c) / (QM2(JUS) + q)
        QM2(JUS) += q ;  QCONTORM(JUS,:) += contorno
    Varias contribuicoes para o mesmo JUS no mesmo nivel sao somadas antes."""
    T = M.NTMUSK
    R = M.R
    uj, inv = np.unique(jj, return_inverse=True)
    Sq = np.zeros(uj.size, dtype=np.float64); np.add.at(Sq, inv, q)
    Sb = np.zeros(uj.size); np.add.at(Sb, inv, cb)
    Si = np.zeros(uj.size); np.add.at(Si, inv, ci)
    Ss = np.zeros(uj.size); np.add.at(Ss, inv, cs)
    Sbnd = np.zeros((uj.size, T + 1)); np.add.at(Sbnd, inv, bnd)
    qm_old = M.QM2[uj]
    qm_new = (qm_old + Sq).astype(R)
    if subst and src is not None:
        for kk, ic in enumerate(src):
            if M.IT > 1 and (ic + 1) in subst:
                isb = subst[ic + 1]
                if M.QLIDO[isb, it] >= 0.0:
                    pos = np.searchsorted(uj, jj[kk])
                    qm_new[pos] = M.QLIDO[isb, it]
    with np.errstate(invalid="ignore", divide="ignore"):
        M.PMB2[uj] = ((M.PMB2[uj] * qm_old + Sb) / qm_new).astype(R)
        M.PMI2[uj] = ((M.PMI2[uj] * qm_old + Si) / qm_new).astype(R)
        M.PMS2[uj] = ((M.PMS2[uj] * qm_old + Ss) / qm_new).astype(R)
    M.QM2[uj] = qm_new
    M.QCONTORM[uj, :T + 1] = (M.QCONTORM[uj, :T + 1] + Sbnd).astype(R)
