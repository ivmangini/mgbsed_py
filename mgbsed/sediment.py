"""
Modulo de sedimentos (MGB-SED).

Traducao de: SED_INICIAL.F90, SED_VELOCITY.f90 (FALLVELOCITY), SED_BACIA.F90,
SED_REDE.F90, SED_PROPAG.F90, SED_CT.f90, SED_NEWTON_H.f90.

Estrategia
----------
* SED_BACIA (MUSLE por URH) e vetorizada em (NC,NU); as somas por URH seguem
  a ordem do Fortran.
* SED_REDE/SED_PROPAG: a unica dependencia entre minibacias e a carga de
  montante CARGM(JUS). As minibacias sao processadas por nivel topologico
  (todas do mesmo nivel juntas), o que preserva exatamente a dependencia.

Precisao: as variaveis REAL(8) do Fortran (PSC, QSAUX, CARGM, DEPTREC, EROSTREC,
DEPT, EROT, iENTRA, iSAI, DEPaux, EROSaux) sao float64; as demais usam M.R.

Particularidades do Fortran reproduzidas (ver README, secao "Comportamentos
herdados do codigo original"):
* SED_CT (Meyer-Peter & Muller, flag_TC=1): ``Kr=26/(D90**(1/6))`` e
  ``(G/9.82)**(1/3)`` usam divisao inteira (1/6 = 0, 1/3 = 0), logo Kr=26 e o
  denominador vale 1. Quando a base de ``**1.5`` e negativa o resultado e NaN,
  e as comparacoes com NaN sao falsas (ramo de erosao, que remobiliza todo o
  deposito CSdep da areia).
* SED_BACIA: QCELaux usa QCEL2 do passo anterior (REDE ainda nao rodou).
"""
from __future__ import annotations

import numpy as np

from .fortran_io import caldat
from . import fmath as fm


# ---------------------------------------------------------------------------
# SED_INICIAL + FALLVELOCITY
# ---------------------------------------------------------------------------
def fallvelocity(M):
    R = M.R
    M.WSP = M.zeros(3)
    for i in range(3):
        FatM = R(53.5) * fm.exp(R(-0.65) * M.Fform[i])
        FatN = R(5.65) * fm.exp(R(-2.5) * M.Fform[i])
        FatP = R(0.7) + R(0.9) * M.Fform[i]
        DCar = M.DMP[i] * fm.pw((R(2.65) - R(1.)) * R(9.81) / fm.pw(M.VISC, R(2)), R(1.) / R(3.))
        w = fm.pw(fm.sqrt(R(0.25) + fm.pw(R(4.) * FatN * fm.pw(DCar, R(3)) / (R(3.) * fm.pw(FatM, R(2))),
                                         R(1.) / FatP)) - R(0.5), FatP)
        M.WSP[i] = ((FatM * M.VISC) / (FatN * M.DMP[i])) * w


def alloca_vars_sed(M):
    """ALLOCA_VARS_SED(0)."""
    NB, NC = M.NB, M.NC
    M.ALFsed = M.zeros(NB); M.BETsed = M.zeros(NB); M.AuxTKS = M.zeros(NB)
    M.CSSbase = M.zeros(NB, NB)
    M.CARGAaux = M.zeros(NC)
    M.VEH = M.zeros(NC); M.BEH = M.zeros(NC)


def sed_inicial(M):
    R = M.R
    NC = M.NC
    M.CTS = M.zeros(NC, 3)
    M.Fform = np.array([0.7, 0.7, 0.7], dtype=R)
    M.DMP = np.array([0.0008, 0.000016, 0.000001], dtype=R)
    M.CSM1 = M.zeros(NC, 3); M.CSM2 = M.zeros(NC, 3)
    M.CSJ1 = M.zeros(NC, 3); M.CSJ2 = M.zeros(NC, 3)
    M.VolTREC1 = M.zeros(NC); M.VolTREC2 = M.zeros(NC)
    M.QSJ2 = M.zeros(NC, 3)
    M.FracS = M.zeros(NC, 3)
    M.DEPTREC = np.zeros((NC, 3)); M.DEPT = np.zeros((NC, 3))
    M.EROSTREC = np.zeros((NC, 3)); M.EROT = np.zeros((NC, 3))
    M.CARGM = np.zeros((NC, 3))
    M.CSdep = np.zeros((NC, 3), dtype=np.float32)            # REAL(4)
    M.iENTRA = np.zeros((NC, 3)); M.iSAI = np.zeros((NC, 3))
    M.inFL = M.zeros(NC, 3); M.outFL = M.zeros(NC, 3)
    M.BalQFL = M.zeros(NC, 2); M.VTOTAL = M.zeros(NC, 2)
    M.VFL1 = M.zeros(NC); M.VFL2 = M.zeros(NC)
    M.CFL1 = M.zeros(NC, 3); M.CFL2 = M.zeros(NC, 3)
    M.DFL = M.zeros(NC, 3)
    M.VISC = fm.pw(R(10.), R(-6))
    fallvelocity(M)
    M.FracS[:] = R(1.) / R(3.)
    ib = M.IBAC
    sel = (ib <= M.SUBfim - 1) & (ib >= M.SUBini - 1)
    M.sed_active = sel
    M.CONTIC = np.nonzero(sel)[0] + 1                       # codigos 1-based
    M.nSEDmini = M.CONTIC.size
    M.iSEDaux = M.nSEDmini
    # posicao "iSEDaux" de cada minibacia ativa
    M.sedpos = -np.ones(NC, int)
    M.sedpos[sel] = np.arange(M.nSEDmini)
    SfTREC = np.maximum(M.DECL, R(0.00001))
    if M.hdFLAG0 == 0:
        HmJ = fm.pw((M.RUGMAN * M.QJ2) / (M.BRIO * fm.pw(SfTREC, R(0.5))), R(3.0) / R(5.0))
    else:
        HmJ = fm.pw((M.RUGMAN * np.abs(M.QJ2)) / (M.BRIO * fm.pw(SfTREC, R(0.5))), R(3.0) / R(5.0))
    M.VolTREC2[sel] = (HmJ * M.BRIO * M.SRIO * R(1000.))[sel]
    M.CSJ2[sel] = 0.0
    M.FracS[sel] = 0.0
    M.Qmini = M.zeros(M.nSEDmini, M.NT)


# ---------------------------------------------------------------------------
# SED_BACIA - producao de sedimentos nas minibacias (MUSLE)
# ---------------------------------------------------------------------------
# Valores padrao da MUSLE (Williams): os mesmos do bloco "!VALORES PADRAO" de
# SED_BACIA.F90 (linhas 130-133, comentado no Fortran).
ALF_PADRAO = 11.8
BET_PADRAO = 0.56
AUXTKS_PADRAO = 1.0


def set_musle_coefs_simulation(M):
    """Coeficientes alfa, beta e AuxTKS por sub-bacia usados em SED_BACIA (ICALIB=0).

    No Fortran original esses valores sao escritos no proprio codigo-fonte
    (SED_BACIA.F90). Aqui:
      * se a pasta de entrada tiver PAR_MUSLE_SED.txt, os valores sao lidos dele
        (io_input.le_par_musle_sed);
      * senao, usam-se os valores padrao da MUSLE (11,8; 0,56; 1,0) em todas as
        sub-bacias, os mesmos do bloco "VALORES PADRAO" comentado no Fortran.
    """
    NB = M.NB
    if getattr(M, "PAR_MUSLE_SED", None) is not None:
        alf, bet, atks = M.PAR_MUSLE_SED
        M.ALFsed = alf.astype(M.R).copy()
        M.BETsed = bet.astype(M.R).copy()
        M.AuxTKS = atks.astype(M.R).copy()
    else:
        M.ALFsed = M.full(NB, ALF_PADRAO)
        M.BETsed = M.full(NB, BET_PADRAO)
        M.AuxTKS = M.full(NB, AUXTKS_PADRAO)


def sed_bacia(M, calib_dsup=None, calib_tks=None):
    """SED_BACIA (ICALIB=0) e SED_CALIB_BACIA (calib_dsup/calib_tks dados).

    Retorna dict com as linhas a gravar (PSCaux, QSUPaux, QCELaux).
    """
    R = M.R
    NC, NU = M.NC, M.NU
    DTP = R(M.DTP)
    if M.ICALIB == 0:
        set_musle_coefs_simulation(M)
    M.PSC[:] = 0.0
    mm, dd, yy = caldat(M.IDINI + M.IT - 1)
    M.IMES, M.IDIA, M.IANO = mm, dd, yy
    ib = M.IBAC
    cell_ok = (M.sbtFLAG != 1) & M.sed_active
    valid = (M.NPXU > 0) & (M.PUSO >= R(0.0001)) & cell_ok[:, None]

    ALF = M.ALFsed[ib][:, None]; BET = M.BETsed[ib][:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        if calib_dsup is None:
            DSUP = M.QSAUX[:, 2:].astype(R)
        else:
            DSUP = calib_dsup.astype(R)
        if M.ICALIB == 3 and calib_dsup is None and getattr(M, "DSUP_SC", None) is not None:
            M.DSUP_SC[:, M.IT - 1, :][valid] = DSUP[valid]          # antes da correcao DSUP<0
        npx = np.where(M.NPXU > 0, M.NPXU, 1)
        APIX = (M.ACEL[:, None] * (M.PUSO / R(100.))) / npx.astype(R)
        APIX = APIX * R(100.)                                   # ha
        FCTE = ALF * fm.pw(APIX, BET) * M.Kusle[ib] * M.Cusle[ib] * M.Pusle[ib] * M.Rgros[ib] * M.LSAcu
        APIX = APIX / R(100.)                                   # km2
        DSUP = np.where(DSUP < 0.0, R(0.0), DSUP)
        P = M.P[:, None]
        if calib_dsup is None:
            QPIC = np.where(P > 0.0, ((DSUP / P) * (P / R(24.)) * APIX) / R(360.0), R(0.0))
        else:   # SED_CALIB_BACIA usa /3.6 (sic)
            QPIC = np.where(P > 0.0, ((DSUP / P) * (P / R(24.)) * APIX) / R(3.6), R(0.0))
        add = FCTE * fm.pw(DSUP, BET) * fm.pw(QPIC, BET)
    M.SLU[valid] = (M.SLU + add)[valid]

    c = cell_ok
    sS = M.SLU.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        FRAC = np.where((sS != 0.0)[:, None], M.SLU / sS[:, None], R(0.0)).astype(R)
    if M.ICALIB == 3 and calib_tks is None and getattr(M, "TKS_SC", None) is not None:
        M.TKS_SC[c] = M.QSAUX[c, 0]
    if calib_tks is None:
        TKSed = (M.AuxTKS[ib] * M.QSAUX[:, 0]).astype(R)
    else:
        TKSed = (M.AuxTKS[ib] * calib_tks).astype(R)
    with np.errstate(divide="ignore", invalid="ignore"):
        QSU = M.SLU / TKSed[:, None]
    SLXU = M.SLU - QSU * DTP
    negx = SLXU < 0.0
    QSU = np.where(negx, M.SLU / DTP, QSU)
    SLUn = np.where(negx, R(0.0), SLXU)
    M.QSU[c] = QSU[c]
    M.SLU[c] = SLUn[c]
    M.PSU[c] = (M.QSU * DTP)[c]
    psc1 = (M.QSU.sum(axis=1) * DTP)
    M.PSC[c, 0] = psc1[c]
    # fracoes granulometricas (ate NU-1, como no Fortran) - PSC e REAL(8)
    PSC2 = np.zeros(NC); PSC3 = np.zeros(NC); PSC4 = np.zeros(NC)
    pa2 = np.zeros(NC, dtype=np.float32); pa3 = np.zeros(NC, dtype=np.float32); pa4 = np.zeros(NC, dtype=np.float32)
    psc1_8 = M.PSC[:, 0]
    psc1_4 = psc1_8.astype(np.float32)                         # PSCaux e REAL(4)
    for iu in range(NU - 1):
        fr = FRAC[:, iu]
        ma = M.Mareia[ib, iu]; ms = M.Msilte[ib, iu]; mg = M.Margila[ib, iu]
        PSC2 = PSC2 + psc1_8 * fr * ma / 100
        PSC3 = PSC3 + psc1_8 * fr * ms / 100
        PSC4 = PSC4 + psc1_8 * fr * mg / 100
        pa2 = (pa2 + psc1_4 * fr * ma / np.float32(100)).astype(np.float32)
        pa3 = (pa3 + psc1_4 * fr * ms / np.float32(100)).astype(np.float32)
        pa4 = (pa4 + psc1_4 * fr * mg / np.float32(100)).astype(np.float32)
    M.PSC[c, 1] = PSC2[c]; M.PSC[c, 2] = PSC3[c]; M.PSC[c, 3] = PSC4[c]
    return {
        "PSC1": psc1_4[c], "PSC2": pa2[c], "PSC3": pa3[c], "PSC4": pa4[c],
        "QSUP": M.QSUP[c].astype(np.float32), "QCEL": M.QCEL2[c].astype(np.float32),
    }


# ---------------------------------------------------------------------------
# SED_CT - capacidade de transporte (vetorizada)
# ---------------------------------------------------------------------------
def sed_ct(M, ic, Uat, VmTREC, SfTREC, HmJ):
    """Retorna CTS (k,3) para as minibacias ``ic`` (array)."""
    R = M.R
    QJ2 = M.QJ2[ic]; BRIO = M.BRIO[ic]
    k = ic.size
    CTS = M.zeros(k, 3)
    flag = M.flag_TC
    DMP = M.DMP
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        if flag == 3:     # Rottner
            VEH = R(3.2808) * QJ2 / (BRIO * HmJ)
            HEH = R(3.2808) * HmJ
            grEH = R(32.2175); SG = R(2.65); GEH = R(168.70)
            BEH = BRIO * R(3.2808)
            for i in range(3):
                # (2/3) e divisao inteira = 0  -> x**0 = 1
                t23 = (DMP[i] * R(3.2808) / HEH) ** 0
                c = BEH * GEH * (((SG - 1) * grEH * (HEH ** 3)) ** R(0.5)) * \
                    (((VEH / (((SG - 1) * grEH * HEH) ** R(0.5))) * (R(0.667) * t23 + R(0.14)) - R(0.778) * t23) ** 3)
                c = c * R(0.0454) / (10 ** 3)
                CTS[:, i] = c / QJ2
            M.VEH[ic] = VEH; M.BEH[ic] = BEH
        if flag == 4:     # Schoklitsch
            for i in range(3):
                QoS = (R(0.00532) * DMP[i] * R(39.37)) / (SfTREC ** (4 // 3))
                c = (R(86.7) * (SfTREC ** R(1.5)) * (R(35.31) * QJ2 - BRIO * R(3.2808) * QoS)) / ((R(0.008) * R(39.37)) ** R(0.5))
                c = c * R(0.0454) / (10 ** 3)
                CTS[:, i] = c / QJ2
        if flag == 5:     # Ackers & White (areia)
            SG = R(2.65); alfaAW = R(10); D35 = R(0.0003)
            RhAW = (BRIO * HmJ) / (BRIO + 2 * HmJ)
            UXAW = (R(9.81) * RhAW * SfTREC) ** R(0.5)
            DXAW = D35 * (((R(9.82) * (SG - 1)) / (M.VISC ** 2)) ** (1 // 3))
            nAW = R(1) - R(0.56) * np.log10(DXAW)
            AAW = R(0.23) / (DXAW ** R(0.5)) + R(0.14)
            mAW = R(9.66) / DXAW + R(1.34)
            CA = R(10) ** (R(2.86) * np.log10(DXAW) - (np.log10(DXAW) ** 2) - R(3.53))
            Fgr = (UXAW / (R(9.82) * D35 * (SG - 1))) * ((VmTREC / (R(5.657) * np.log10(alfaAW * HmJ / D35))) ** (R(1) - nAW))
            Ggr = CA * ((Fgr / AAW - 1) ** mAW)
            c = (10 ** 6) * Ggr * SG * D35 / (HmJ * ((UXAW / VmTREC) ** nAW))
            c = c / (10 ** 6)
            CTS[:, 0] = c; CTS[:, 1] = c; CTS[:, 2] = c
        if flag == 6:     # Engelund & Hansen
            VEH = R(3.2808) * QJ2 / (BRIO * HmJ)
            HEH = R(3.2808) * HmJ
            grEH = R(32.2175); GEH = R(168.70)
            BEH = BRIO * R(3.2808)
            for i in range(3):
                c = (BEH * R(0.05) * GEH * (VEH ** 2) * (HEH ** R(0.5)) * (SfTREC ** R(1.5))) / \
                    ((DMP[i] * R(3.2808)) * grEH * (R(1.65) ** 2))
                c = c * R(0.0454) / 1000
                CTS[:, i] = c / QJ2
            M.VEH[ic] = VEH; M.BEH[ic] = BEH
        if flag == 1:     # Meyer-Peter & Muller
            Kr = R(26) / (R(1) ** 0)          # D90**(1/6) com 1/6 inteiro = 0  ->  Kr = 26
            Ks = R(1) / M.RUGMAN[ic]
            GSS = R(1.65)
            G = R(1)
            for i in range(3):
                QLT = QJ2 * R(0.6)
                base = (((G * (QLT / QJ2) * fm.pw(Ks / Kr, R(1.5))) * HmJ * SfTREC - R(0.047) * GSS * DMP[i]) /
                        R(1.0))                     # (G/9.82)**(1/3) com 1/3 inteiro = 0 -> 1
                c = BRIO * fm.pw(base, R(1.5))
                CTS[:, i] = c / QJ2
        if flag == 2:     # Yang
            for i in range(3):
                REX = np.maximum(Uat * DMP[i] / M.VISC, R(1.2))
                UcWs = np.where(REX < R(70.), R(0.66) + R(2.5) / (fm.log10(REX) - R(0.06)), R(2.05))
                if DMP[i] <= 0.002:
                    FMY = R(5.435) - R(0.286) * fm.log10(M.WSP[i] * DMP[i] / M.VISC) - R(0.457) * fm.log10(Uat / M.WSP[i])
                    FNY = R(1.799) - R(0.409) * fm.log10(M.WSP[i] * DMP[i] / M.VISC) - R(0.314) * fm.log10(Uat / M.WSP[i])
                else:
                    FMY = R(6.681) - R(0.633) * fm.log10(M.WSP[i] * DMP[i] / M.VISC) - R(4.816) * fm.log10(Uat / M.WSP[i])
                    FNY = R(2.874) - R(0.305) * fm.log10(M.WSP[i] * DMP[i] / M.VISC) - R(0.282) * fm.log10(Uat / M.WSP[i])
                aux = VmTREC * SfTREC / M.WSP[i] - UcWs * SfTREC
                LOG = FMY + FNY * fm.log10(np.where(aux > 0, aux, R(1)))
                c = fm.pw(R(10.), LOG) / R(1.0e6) * M.FracS[ic, i]
                CTS[:, i] = np.where(aux > 0.0, c, R(0.0))
    return CTS.astype(R)


# ---------------------------------------------------------------------------
# SED_REDE + SED_PROPAG
# ---------------------------------------------------------------------------
def sed_rede(M):
    """SED_REDE.F90 (com SED_PROPAG vetorizado por nivel). Retorna as linhas
    a gravar (CONC_*, DEPdin)."""
    R = M.R
    NC = M.NC
    DTP = R(M.DTP)
    M.CSM1[:] = M.CSM2
    M.CSJ1[:] = M.CSJ2
    M.VolTREC1[:] = M.VolTREC2
    M.CFL1[:] = M.CFL2
    M.VFL1[:] = M.VFL2
    M.VFL2[:] = 0.0
    M.CARGM[:] = 0.0
    M.CSJ2[:] = 0.0
    M.CFL2[:] = 0.0
    M.VTOTAL[:] = 0.0
    Qlim = R(0.01)

    if M.hdFLAG0 > 0:
        M.QM1[:NC] = M.QM2[:NC]
        q = M.Q2fl
        qm = np.zeros(NC)
        for ic in range(NC):
            ne = M.MINIMONT[ic, 0]
            if ne > 0:
                qm[ic] = q[M.MINIMONT[ic, 1:1 + ne] - 1].sum()
        M.QM2[:NC] = qm

    act = M.sed_active
    # grandezas hidraulicas por minibacia (dependem apenas de QJ2 - ja conhecido)
    HmJ = M.zeros(NC); Sf = M.zeros(NC); Vm = M.zeros(NC); Uat = M.zeros(NC)
    QtFL = M.zeros(NC); HFLSD = M.zeros(NC)
    with np.errstate(divide="ignore", invalid="ignore"):
        if M.hdFLAG0 == 0:
            Sf[:] = np.maximum(M.DECL, R(0.00001))
            HmJ[:] = fm.pw((M.RUGMAN * M.QJ2) / (M.BRIO * fm.pw(Sf, R(0.5))), R(3.0) / R(5.0))
            M.VolTREC2[act] = (HmJ * M.BRIO * M.SRIO * R(1000.))[act]
            M.VTOTAL[act, 0] = M.VolTREC2[act]
            Vm[:] = np.where(M.QJ2 > Qlim, M.QJ2 / (M.BRIO * HmJ), R(0.0))
            RHm = HmJ
            Uat[:] = fm.sqrt(R(9.81) * RHm * Sf)
            M.VFL1[act] = 0.0
            M.VFL2[act] = 0.0
        else:
            from .inertial import sed_hydraulics_inertial
            sed_hydraulics_inertial(M, act, HmJ, Sf, Vm, Uat, QtFL, HFLSD, Qlim)

    for cells in M.levels:
        cells = cells[act[cells]]
        if cells.size == 0:
            continue
        _sed_propag(M, cells, HmJ[cells], Vm[cells], Sf[cells], Uat[cells], Qlim,
                    QtFL[cells], HFLSD[cells])

    if M.hdFLAG0 > 0:
        M.QJ1[:] = M.QJ2

    pos = M.CONTIC - 1
    return {
        "CA": (M.CSJ2[pos, 0] * R(10.) ** R(6.)).astype(np.float32),
        "CS": (M.CSJ2[pos, 1] * R(10.) ** R(6.)).astype(np.float32),
        "CC": (M.CSJ2[pos, 2] * R(10.) ** R(6.)).astype(np.float32),
        # (sic) Fortran grava CSdep(IC,1) para IC=1..iSEDaux (indice da minibacia)
        "DEP": M.CSdep[:M.nSEDmini, 0].astype(np.float32),
    }


def _sed_propag(M, ic, HmJ, VmJ, SfTREC, UatJ, Qlim, QtFL, HFLSD):
    R = M.R
    DTP = R(M.DTP)
    DTP8 = float(M.DTP)
    PONDt = R(1.0)
    porosi = 0.0
    densiSED = 2.0
    QM2 = M.QM2[ic]; QJ2 = M.QJ2[ic]; QM1 = M.QM1[ic]; QJ1 = M.QJ1[ic]
    Vol1 = M.VolTREC1[ic]; Vol2 = M.VolTREC2[ic]
    PSC = M.PSC[ic]
    riv = M.NSUBT[ic] > 0
    qpos = QJ2 > Qlim
    qmpos = QM2 > Qlim
    CSM2 = M.zeros(ic.size, 3)
    CSJ2 = M.zeros(ic.size, 3)
    CSJ1 = M.CSJ1[ic]; CSM1 = M.CSM1[ic]
    CFL1 = M.CFL1[ic]; VFL1 = M.VFL1[ic]; VFL2X = M.VFL2[ic]
    CFL2 = M.CFL2[ic].copy()
    DFL = M.DFL[ic].copy()
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # =================== concentracao de entrada (balanco) ==============
        for i in range(3):
            CSM2[:, i] = np.where(qmpos, (M.CARGM[ic, i] / (QM2 * DTP)).astype(R), R(0.0))
        # --- com rio ---
        # areia
        EQ1 = PONDt * QM2 * CSM2[:, 0]
        EQ2 = (R(1.0) - PONDt) * (QJ1 * CSJ1[:, 0] - QM1 * CSM1[:, 0])
        EQ3 = R(0.0)
        EQ4 = Vol1 * (CSJ1[:, 0]) / DTP
        EQ5 = Vol2 / DTP + PONDt * QJ2
        c = ((((PSC[:, 1] / DTP8 + EQ1) - EQ2) - EQ3) + EQ4) / EQ5
        c = c.astype(R)
        c = np.where(c < 0.0, R(0.0), c)
        cr_sand = np.where(qpos, c, R(0.0))
        cr = [cr_sand, None, None]
        for i in (1, 2):
            EQ1 = PONDt * QM2 * (CSM2[:, i])
            EQ2 = (R(1.0) - PONDt) * (QJ1 * (CSJ1[:, i]) - QM1 * (CSM1[:, i]))
            EQ3 = (Vol1 * (CSJ1[:, i]) / DTP)
            EQ4 = (Vol2 / DTP + PONDt * QJ2)
            # ramo QtFL >= 0
            EQ5a = R(0.5) * QtFL * (R(0.5) * (CSM2[:, i] + CSJ1[:, i]))
            EQ6a = (R(0.5) * QtFL)
            ca = np.abs(((((EQ1 - EQ2) + EQ3) + (PSC[:, i + 1] / DTP8)) - EQ5a) / (EQ4 + EQ6a)).astype(R)
            ca = np.where(ca < 0.0, R(0.0), ca)
            # ramo QtFL < 0
            EQ5b = R(0.5) * np.minimum(QtFL * CFL1[:, i], CFL1[:, i] * VFL1 / DTP)
            cb = (((((EQ1 - EQ2) + EQ3) + (PSC[:, i + 1] / DTP8)) - EQ5b) / (EQ4 + R(0.0))).astype(R)
            cb = np.where(cb < 0.0, R(0.0), cb)
            brA = qpos & (QtFL >= 0.0)
            brB = qpos & (QtFL < 0.0)
            ci = np.where(brA, ca, np.where(brB, cb, R(0.0))).astype(R)
            cr[i] = ci
            # planicie (so atua no modelo inercial)
            GFLa = np.maximum(CFL1[:, i] * VFL1 + R(0.5) * QtFL * ((ci) + (CSM2[:, i])) * DTP, R(0.0))
            GFLb = np.maximum(CFL1[:, i] * VFL1 + QtFL * CFL1[:, i] * DTP, R(0.0))
            GFL = np.where(brA, GFLa, GFLb)
            LIMd = np.minimum(M.WSP[i] * DTP / HFLSD, R(1.0))
            alt = np.where(brA, GFL, CFL1[:, i] * VFL1)
            DFLaux = np.where(HFLSD > 0.0, GFL * LIMd,
                              np.where(CFL1[:, i] * VFL1 > 0.0, alt, R(0.0)))
            cfl2 = np.where(VFL2X > 0.0, np.maximum((GFL - DFLaux) / VFL2X, R(0.0)), R(0.0))
            sel = (brA | brB) & riv
            CFL2[:, i] = np.where(sel & (brA | (VFL2X > 0.0)), cfl2, CFL2[:, i])
            DFL[:, i] = np.where(sel, DFL[:, i] + DFLaux, DFL[:, i])
            # acumuladores de troca com a planicie
            p = M.sedpos[ic]
            m_in = riv & (QtFL > 0.0)
            m_out = riv & (QtFL < 0.0)
            M.inFL[p[m_in], i] += (R(0.5) * QtFL * (ci + CSM2[:, i]) * DTP)[m_in]
            M.outFL[p[m_out], i] += (QtFL * CFL1[:, i] * DTP)[m_out]
            if i == 1:
                pass
        # BalQFL acumula uma vez por particula i=2,3 (como no laco do Fortran)
        p = M.sedpos[ic]
        m_in = riv & (QtFL > 0.0)
        m_out = riv & (QtFL < 0.0)
        M.BalQFL[p[m_in], 0] += (2 * QtFL)[m_in]
        M.BalQFL[p[m_out], 1] += (2 * QtFL)[m_out]

        # --- sem rio ---
        cn = [np.where(qpos, (PSC[:, i + 1] / (QJ2 * DTP)).astype(R), R(0.0)) for i in range(3)]
        for i in range(3):
            CSJ2[:, i] = np.where(riv, cr[i], cn[i])
        M.CFL2[ic] = CFL2
        M.DFL[ic] = DFL

        # =================== capacidade de transporte =======================
        FracS = M.zeros(ic.size, 3)
        sumc = CSJ2.sum(axis=1)
        fr_ok = qpos & (CSJ2[:, 0] > 0.0)
        for i in range(3):
            FracS[:, i] = np.where(fr_ok, CSJ2[:, i] / sumc, R(0.0))
        M.FracS[ic] = FracS
        CTS = M.zeros(ic.size, 3)
        if qpos.any():
            sub = np.nonzero(qpos)[0]
            CTS[sub] = sed_ct(M, ic[sub], UatJ[sub], VmJ[sub], SfTREC[sub], HmJ[sub])
        M.CTS[ic] = CTS

        # =================== deposicao / erosao ============================
        VQ = np.where(riv, Vol2 + QJ2 * DTP, QJ2 * DTP)        # volume de referencia
        jus = M.jus[ic]
        p = M.sedpos[ic]
        for i in range(3):
            DEPaux = np.zeros(ic.size)
            EROSaux = np.zeros(ic.size)
            csj = CSJ2[:, i].copy()
            cts = CTS[:, i]
            csdep = M.CSdep[ic, i].copy()
            dep_m = qpos & (csj >= cts)
            ero_m = qpos & ~(csj >= cts)                     # inclui CTS = NaN
            dmp = M.DMP[i]
            # deposicao
            if dmp <= 0.000062:
                pass                                          # finos passam
            else:
                if dmp <= 0.000125:
                    DeT = HmJ
                elif dmp <= 0.00025:
                    DeT = HmJ / R(2.5)
                else:
                    DeT = HmJ / R(11.24)
                LIMd = np.minimum(M.WSP[i] * DTP / DeT, R(1.0))
                dep = ((csj - cts) * VQ * LIMd).astype(np.float64)
                dep2 = ((csj - cts) * VQ * (R(1.0) - LIMd) + cts * VQ).astype(R)
                if dmp <= 0.00025:
                    newc = cts
                else:
                    newc = (dep2 / VQ).astype(R)
                DEPaux = np.where(dep_m, dep, DEPaux)
                csj = np.where(dep_m, newc, csj)
                M.DEPTREC[ic, i] = np.where(dep_m, M.DEPTREC[ic, i] + dep, M.DEPTREC[ic, i])
                csdep = np.where(dep_m, (csdep + dep).astype(np.float32), csdep)
                if i == 0:
                    h2 = M.HRIO2[ic]
                    m = dep_m & riv
                    M.HRIO1[ic] = np.where(m, h2, M.HRIO1[ic])
                    M.HRIO2[ic] = np.where(m, np.maximum(0.0, h2 - dep / ((1 - porosi) * densiSED * M.SRIO[ic] * M.BRIO[ic])), h2)
            # erosao
            if dmp > 0.000062:
                LIMe = np.minimum((R(1.0) + fm.exp(R(-1.0))) - fm.exp(-M.SRIO[ic] * R(1000.) / (R(30.) * HmJ)), R(1.0))
                ero = ((cts - csj) * VQ * LIMe).astype(np.float64)
                ok = ero <= csdep                             # NaN -> False
                c_ok = (csj + ero / VQ).astype(R)
                c_all = (csj + csdep / VQ).astype(R)
                er_a = ero_m & ok
                er_b = ero_m & ~ok
                csj = np.where(er_a, c_ok, np.where(er_b, c_all, csj))
                EROSaux = np.where(er_a, ero, np.where(er_b, csdep.astype(np.float64), EROSaux))
                M.EROSTREC[ic, i] = np.where(ero_m, M.EROSTREC[ic, i] + EROSaux, M.EROSTREC[ic, i])
                csdep = np.where(er_a, (csdep - ero).astype(np.float32),
                                 np.where(er_b, np.float32(0.0), csdep))
                if i == 0:
                    h2 = M.HRIO2[ic]
                    m = ero_m & riv
                    M.HRIO1[ic] = np.where(m, h2, M.HRIO1[ic])
                    M.HRIO2[ic] = np.where(m, h2 + EROSaux / ((1 - porosi) * densiSED * M.SRIO[ic] * M.BRIO[ic]), h2)
            # reinicia geometria se deposito zerado (ramo de erosao, com rio)
            m = ero_m & riv & (csdep == 0.0)
            M.HRIO1[ic] = np.where(m, M.HRIO[ic], M.HRIO1[ic])
            M.HRIO2[ic] = np.where(m, M.HRIO[ic], M.HRIO2[ic])
            # sem vazao: tudo deposita
            nq = ~qpos
            dep_nq = np.where(riv, PSC[:, i + 1] + (CSM2[:, i] * np.abs(QM2) * DTP).astype(np.float64)
                              + (CSJ1[:, i] * Vol1).astype(np.float64),
                              PSC[:, i + 1])
            DEPaux = np.where(nq, dep_nq, DEPaux)
            csdep = np.where(nq, (csdep + dep_nq).astype(np.float32), csdep)
            M.DEPTREC[ic, i] = np.where(nq, M.DEPTREC[ic, i] + dep_nq, M.DEPTREC[ic, i])
            csj = np.where(nq, R(0.0), csj)

            M.CSdep[ic, i] = csdep
            CSJ2[:, i] = csj
            QSJ2 = np.where(qpos, csj * QJ2, R(0.0)).astype(R)
            M.QSJ2[ic, i] = QSJ2
            h = jus >= 0
            if h.any():
                np.add.at(M.CARGM[:, i], jus[h], (QSJ2 * DTP)[h].astype(np.float64))
            ent = np.where(riv, M.CARGM[ic, i] + PSC[:, i + 1] + EROSaux, PSC[:, i + 1] + EROSaux)
            M.iENTRA[p, i] += ent
            M.iSAI[p, i] += (csj * QJ2 * DTP).astype(np.float64) + DEPaux
        M.CSM2[ic] = CSM2
        M.CSJ2[ic] = CSJ2
        M.DEPT[p] = M.DEPTREC[ic]
        M.EROT[p] = M.EROSTREC[ic]


# ---------------------------------------------------------------------------
# NEWTON_H (nao chamada no fluxo de simulacao do MGB-SED v1.0)
# ---------------------------------------------------------------------------
def newton_h(M, ic, QmTRECx, SfTRECx):
    R = M.R
    B = M.BRIO[ic]; n = M.RUGMAN[ic]
    HK = R(0.01) * B
    ERROH = R(9999.)
    itera = 0
    TOL = R(0.000001)
    Hm = HK
    while ERROH > TOL:
        itera += 1
        QK = ((B * HK) ** (R(5.) / R(3.))) * (SfTRECx ** (R(1.) / R(2.))) / (((B + R(2.) * HK) ** (R(2.) / R(3.))) * n)
        Hm = HK - (R(1.) - QmTRECx / QK) * (R(3.) * HK * (B + R(2.) * HK)) / (R(5.) * B + R(6.) * HK)
        ERROH = np.sqrt((Hm - HK) ** R(2.))
        HK = Hm
        if itera % 30 == 0:
            TOL = TOL * 10
    return Hm
