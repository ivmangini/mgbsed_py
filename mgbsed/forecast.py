"""
Modo de previsao (ICALIB=2) com assimilacao de dados.

Traducao de: PREVISAO.f90, ATUALIZA.f90, leprev.f90.

ATENCAO (nao verificado: nao ha saidas de referencia deste modo).
Adaptacoes necessarias (o original e especifico de uma instalacao Windows):
* Caminhos fixos 'c:\\MGB-hora\\output\\...' e 'c:\\MGB-hora\\input\\Previsao.meta'
  foram substituidos pelas pastas de saida/entrada do modelo.
* ARQHID(K) nunca e lido em LEFIX (linha comentada no original); os arquivos
  de previsao recebem o nome 'Previsao_<codigo da minibacia>.txt'.
* A opcao IOPCHU=2 (LEPREV) le um arquivo binario de acesso direto (unidade
  999) que o original nunca abre; aqui o arquivo esperado e 'PREVISAO.PBI'.
"""
from __future__ import annotations

import os

import numpy as np

from . import hydrology as hy
from . import routing as rt
from .fortran_io import FortranFile, caldat


def _alloc_prev(M):
    for nm in ("WPREV", "SIPREV"):
        setattr(M, nm, np.zeros_like(M.W))
    for nm in ("VBASPREV", "VINTPREV", "VSUPPREV", "TAPREV", "QJ2PREV", "QJ1PREV",
               "QCEL1PREV", "QCEL2PREV"):
        setattr(M, nm, np.zeros(M.NC, M.R))
    for nm in ("QM2PREV", "PMB2PREV", "PMI2PREV", "PMS2PREV", "PJB2PREV", "PJI2PREV", "PJS2PREV"):
        setattr(M, nm, np.zeros(M.NC + 1, M.R))
    M.QRIOINIPREV = np.zeros_like(M.QRIOINI)
    M.QCONTORMPREV = np.zeros_like(M.QCONTORM)
    M.QCONTORJPREV = np.zeros_like(M.QCONTORJ)
    M.QPREV = np.full((M.NUMHIDG, 100), -1.0, M.R)


def atualiza(M):
    """ATUALIZA.f90 - correcao dos reservatorios e vazoes pelos erros recentes."""
    R = M.R
    NB = M.NB
    it = M.IT - 1
    POSTORIO = np.zeros(NB, int); POSTOBAS = np.zeros(NB, int)
    POSTOINT = np.zeros(NB, int); POSTOSUP = np.zeros(NB, int)
    POSTORIO[0] = 1; POSTOBAS[0] = 1; POSTOINT[0] = 1; POSTOSUP[0] = 1
    EXPOENTEBAC = np.full(NB, 1.0, R)
    EXPOBAS = np.full(NB, 0.5, R); EXPOINT = np.full(NB, 0.5, R); EXPOSUP = np.full(NB, 0.8, R)
    PLIMBAS = np.full(NB, 0.5, R); PLIMINT = np.full(NB, 0.2, R); PLIMSUP = np.full(NB, 0.2, R)
    NDIASFBAS = np.full(NB, 5, int); NDIASFINT = np.full(NB, 2, int); NDIASFSUP = np.full(NB, 1, int)
    NTMUSK = M.NTMUSK
    n = max(NB, M.NOBS)
    FCORBAS = np.ones(n, R); FCORSUP = np.ones(n, R); FCORINT = np.ones(n, R)

    def fator(i, nd):
        so = sc = R(0.0)
        for iww in range(1, nd + 1):
            t = it + 1 - iww
            if M.QOBS[i, t] > 0.0:
                so += M.QOBS[i, t]; sc += M.QR[i, t]
        return so / sc if sc > 0.0 else R(1.0)

    for i in range(M.NOBS):
        K = M.IBAC[M.IQOBS[i] - 1]
        FCORSUP[i] = fator(i, NDIASFSUP[i])
        FCORINT[i] = fator(i, NDIASFINT[i])
        FCORBAS[i] = fator(i, NDIASFBAS[i])
        FCORBAS[K] = max(min(FCORBAS[K], 1.2), 0.8)
        FCORSUP[K] = max(min(FCORSUP[K], 5.0), 0.2)
        FCORINT[K] = max(min(FCORINT[K], 1.2), 0.8)
    for i in range(M.NOBS):
        FCORBAS[i] = max(min(FCORBAS[i], 1.2), 0.8)
        FCORSUP[i] = max(min(FCORSUP[i], 5.0), 0.2)
        FCORINT[i] = max(min(FCORINT[i], 1.2), 0.8)
    for ic in range(M.NC):
        NTRECH = M.NSUBT[ic]
        i = M.IBAC[ic]
        if POSTOBAS[i] > 0 and M.PJB2[ic] > PLIMBAS[i]:
            f = FCORBAS[POSTOBAS[i] - 1] ** EXPOBAS[i]
            M.VBAS[ic] = M.VBAS[ic] * f * M.PJB2[ic] + M.VBAS[ic] * (1. - M.PJB2[ic])
            M.VBASPREV[ic] = M.VBASPREV[ic] * f * M.PJB2[ic] + M.VBASPREV[ic] * (1. - M.PJB2[ic])
        if POSTOINT[i] > 0 and M.PJI2[ic] > PLIMINT[i]:
            f = FCORINT[POSTOINT[i] - 1] ** EXPOINT[i]
            M.VINT[ic] = M.VINT[ic] * f * M.PJI2[ic] + M.VINT[ic] * (1. - M.PJI2[ic])
            M.VINTPREV[ic] = M.VINTPREV[ic] * f * M.PJI2[ic] + M.VINTPREV[ic] * (1. - M.PJI2[ic])
        if POSTOSUP[i] > 0 and M.PJS2[ic] > PLIMSUP[i]:
            f = FCORSUP[POSTOSUP[i] - 1] ** EXPOSUP[i]
            M.VSUP[ic] = M.VSUP[ic] * f * M.PJS2[ic] + M.VSUP[ic] * (1. - M.PJS2[ic])
            M.VSUPPREV[ic] = M.VSUPPREV[ic] * f * M.PJS2[ic] + M.VSUPPREV[ic] * (1. - M.PJS2[ic])
        for a in (M.VSUP, M.VINT, M.VBAS, M.VSUPPREV, M.VINTPREV, M.VBASPREV):
            a[ic] = max(a[ic], 0.0)
        if POSTORIO[i] > 0 and M.OD[ic] > 1:
            K = M.IQOBS[POSTORIO[i] - 1] - 1
            fc = FCORSUP[POSTORIO[i] - 1]
            PESOBAC = (M.ACUR[ic] / M.ACUR[K]) ** EXPOENTEBAC[i]
            M.QJ2PREV[ic] = (M.QJ2[ic] * fc) * PESOBAC + M.QJ2[ic] * (1. - PESOBAC)
            M.QJ2[ic] = M.QJ2PREV[ic]
            for itr in range(NTRECH + 1):
                M.QRIOINI[ic, itr] = M.QRIOINI[ic, itr] * fc * PESOBAC + M.QRIOINI[ic, itr] * (1. - PESOBAC)
                M.QRIOINIPREV[ic, itr] = M.QRIOINI[ic, itr]
            for lt in range(NTMUSK + 1):
                M.QCONTORM[ic, lt] = (M.QCONTORM[ic, lt] * fc) * PESOBAC + M.QCONTORM[ic, lt] * (1. - PESOBAC)


def previsao(M):
    from .fortran_io import fmt_I, fmt_F
    R = M.R
    hy.condinic(M)
    _alloc_prev(M)
    M.ATIVABACIA[:] = 1
    ITINI = 1
    ISUBSTAUX = M.ISUBST.copy()
    M.QTUDO = np.full((M.NC, M.NT), -1.0, R)
    M.QITUDO = np.full((M.NC, M.NT), -1.0, R)
    M.QRG[:] = -1.0
    out = M.OUTPUT_DIRECTORY
    fprev = []
    for k in range(M.NUMHIDG):
        f = open(os.path.join(out, f"Previsao_{M.IHIDG[k]}.txt"), "w")
        f.write("   DIA   MES   ANO  HORA        Q_(m3/s)\n")
        fprev.append(f)
    meta = FortranFile(os.path.join(M.INPUT_DIRECTORY, "Previsao.meta"))
    meta.skip(7)
    IANOPREV, IMESPREV, IDIAPREV, IHORAPREV, IHORIZ = [int(x) for x in meta.read_list(5)]
    IOPCHU = 1
    M.log("Start forecast at: Year, Month, Day, Hour, Horizon")
    M.log(IANOPREV, IMESPREV, IDIAPREV, IHORAPREV, IHORIZ)
    M.log("Simulating and performing data assimilation...")
    prevfile = None
    M.IT = ITINI - 1
    ITEMP = M.NT
    while M.IT < M.NT - IHORIZ + 1:
        M.IT += 1
        IT = M.IT; it = IT - 1
        if IT == ITINI:
            M.SI[:] = 0.0
            M.W[:, :] = M.WM[M.IBAC, :] * R(0.40)
            ib = M.IBAC
            M.VBAS[:] = M.QESP[ib] * M.ACEL * M.CB[ib] * R(3600.)
            M.VINT[:] = 0.0; M.VSUP[:] = 0.0
            M.TA[:] = 20.0
            M.QPREV[:] = -1.0
        else:
            M.W[:] = M.WPREV; M.VBAS[:] = M.VBASPREV; M.VINT[:] = M.VINTPREV; M.VSUP[:] = M.VSUPPREV
            M.TA[:] = M.TAPREV; M.QM2[:] = M.QM2PREV; M.QJ2[:] = M.QJ2PREV; M.SI[:] = M.SIPREV
            M.QJ1[:] = M.QJ1PREV; M.QRIOINI[:] = M.QRIOINIPREV; M.QCONTORM[:] = M.QCONTORMPREV
            M.QCONTORJ[:] = M.QCONTORJPREV; M.QCEL1[:] = M.QCEL1PREV; M.QCEL2[:] = M.QCEL2PREV
            M.PMB2[:] = M.PMB2PREV; M.PMI2[:] = M.PMI2PREV; M.PMS2[:] = M.PMS2PREV
            M.PJB2[:] = M.PJB2PREV; M.PJI2[:] = M.PJI2PREV; M.PJS2[:] = M.PJS2PREV
        # (sic) PREVISAO usa (IT+HORAINI) e nao (IT+HORAINI-1)
        M.JDIA = M.IDINI + int((IT + M.HORAINI) / (86400. / M.DTP))
        M.IMES, M.IDIA, M.IANO = caldat(M.JDIA)
        M.DIAH[it] = M.IDIA; M.MESH[it] = M.IMES; M.ANOH[it] = M.IANO
        M.HORAH[it] = IT % 24 + M.HORAINI
        M.P[:] = M.rainfile.read(IT)
        M.TONTEM[:] = M.TA
        hy.leclima(M)
        hy.celula(M)
        M.ISUBST = ISUBSTAUX
        rt.rede(M)
        M.QR[:, it] = M.QJ2[M.IQOBS - 1]
        M.QRB[:, it] = M.QJ2[M.IEXUT - 1]
        M.QRG[:, it] = M.QJ2[M.IHIDG - 1]
        MWP = M.NC - 1
        M.QB[it] = M.QJ2[MWP] * M.PJB2[MWP]
        M.QBI[it] = M.QB[it] + M.QJ2[MWP] * M.PJI2[MWP]
        M.QBIS[it] = M.QBI[it] + M.QJ2[MWP] * M.PJS2[MWP]
        M.WPREV[:] = M.W; M.VBASPREV[:] = M.VBAS; M.VINTPREV[:] = M.VINT; M.VSUPPREV[:] = M.VSUP
        M.TAPREV[:] = M.TA; M.QM2PREV[:] = M.QM2; M.QJ2PREV[:] = M.QJ2; M.SIPREV[:] = M.SI
        M.QJ1PREV[:] = M.QJ1; M.QRIOINIPREV[:] = M.QRIOINI; M.QCONTORMPREV[:] = M.QCONTORM
        M.QCONTORJPREV[:] = M.QCONTORJ; M.QCEL1PREV[:] = M.QCEL1; M.QCEL2PREV[:] = M.QCEL2
        M.PMB2PREV[:] = M.PMB2; M.PMI2PREV[:] = M.PMI2; M.PMS2PREV[:] = M.PMS2
        M.PJB2PREV[:] = M.PJB2; M.PJI2PREV[:] = M.PJI2; M.PJS2PREV[:] = M.PJS2
        if IT > 48:
            atualiza(M)
        M.JDIA = M.IDINI + int((IT + M.HORAINI) / (86400. / M.DTP))
        M.IMES, M.IDIA, M.IANO = caldat(M.JDIA)
        IHORA = IT % 24 + M.HORAINI
        if (M.IDIA == IDIAPREV and M.IMES == IMESPREV and M.IANO == IANOPREV and IHORA == IHORAPREV):
            M.log("Realizando a previsao...")
            ITEMP = IT
            for ITP in range(IT + 1, IT + IHORIZ + 1):
                IANTE = ITP - ITEMP
                M.JDIA = M.IDINI + int((ITP + M.HORAINI) / (86400. / M.DTP))
                mm, dd, yy = caldat(M.JDIA)
                M.IMES, M.IDIA, M.IANO = mm, dd, yy
                DIA_ANTE, MES_ANTE, ANO_ANTE, HORA_ANTE = dd, mm, yy, ITP % 24 + M.HORAINI
                if IOPCHU == 1:
                    M.P[:] = M.rainfile.read(ITP)
                elif IOPCHU == 2:
                    if prevfile is None:
                        from .io_input import RainFile
                        prevfile = RainFile(os.path.join(M.INPUT_DIRECTORY, "PREVISAO.PBI"), M.NC)
                    M.P[:] = prevfile.read(ITP - ITEMP)
                elif IOPCHU == 3:
                    M.P[:] = 0.0
                else:
                    raise RuntimeError("OPCAO DE CHUVA DESCONHECIDA")
                M.IT = ITP
                M.TONTEM[:] = M.TA
                hy.leclima(M)
                hy.celula(M)
                rt.rede(M)
                M.QR[:, ITP - 1] = M.QJ2[M.IQOBS - 1]
                M.QPREV[:, IANTE - 1] = M.QJ2[M.IHIDG - 1]
                for k in range(M.NUMHIDG):
                    fprev[k].write("".join(fmt_I(v, 6) for v in (DIA_ANTE, MES_ANTE, ANO_ANTE, HORA_ANTE)) +
                                   fmt_F(M.QPREV[k, IANTE - 1], 16, 6) + "\n")
            for k in range(M.NUMHIDG):
                fprev[k].write("".join(fmt_I(v, 6) for v in (DIA_ANTE, MES_ANTE, ANO_ANTE, HORA_ANTE)) +
                               fmt_F(-1.0, 16, 6) + "\n")
            M.IT = ITEMP
            break
    for f in fprev:
        f.close()
    with open(os.path.join(out, "VAZAO.HIG"), "w") as fh, open(os.path.join(out, "QPROP.HIG"), "w") as fp:
        for it in range(M.NT):
            fh.write(fmt_I(it + 1, 10) + "".join(fmt_F(v, 15, 3) for v in M.QRG[:, it]) + "\n")
            fp.write(fmt_I(it + 1, 10) + "".join(fmt_F(v, 15, 3) for v in (M.QB[it], M.QBI[it], M.QBIS[it])) + "\n")
    # binarios de acesso direto (RECL=NC palavras)
    with open(os.path.join(out, "QTUDO_ATUAL.MGB"), "wb") as f1, \
            open(os.path.join(out, "QITUDO_ATUAL.MGB"), "wb") as f3, \
            open(os.path.join(out, "QTUDO_PREV.MGB"), "wb") as f2:
        for it in range(M.NT):
            if it + 1 <= ITEMP:
                f1.write(M.QTUDO[:, it].astype("<f4").tobytes())
                f3.write(M.QITUDO[:, it].astype("<f4").tobytes())
            else:
                f2.write(M.QTUDO[:, it].astype("<f4").tobytes())
    for k in range(M.NUMHIDG):
        with open(os.path.join(out, f"Hidrograma_{M.IHIDG[k]}.txt"), "w") as f:
            f.write("   DIA   MES   ANO  HORA        Q_(m3/s)\n")
            for it in range(M.NT):
                f.write("".join(fmt_I(v, 6) for v in (M.DIAH[it], M.MESH[it], M.ANOH[it], M.HORAH[it])) +
                        fmt_F(M.QRG[k, it], 16, 6) + "\n")
