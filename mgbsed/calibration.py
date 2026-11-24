"""
Calibracao automatica.

Traducao de: CALIBRA.f90 (vazia no original), CalibParam.f90, SEMENTE.f90,
ran1.f90, ALLOCA_CALIB.f90, SED_CALIBRA.f90 (CALIBRA_SED, MOCOM-UA),
SED_CalibParam.f90, SED_CALIB_MODELO.f90, SED_CALIB_BACIA.f90,
SED_CALIB_REDE.f90.

ATENCAO (nao verificado: nao ha saidas de referencia de calibracao).
Problemas do codigo original reproduzidos fielmente e registrados no log:
* CALIBRA (hidrologia) e uma sub-rotina vazia na v1.0.
* CALIBRA_SED: na geracao da populacao inicial a funcao objetivo nunca e
  avaliada (FO permanece zerada); com FO identicas, nenhum ponto e dominado,
  IRMAX=1 e o laco evolutivo termina na primeira iteracao.
* SED_CALIB_REDE depende de Q_SC, HmJx_SC, SfTRECx_SC, VOLTREC2_SC, VFL2_SC e
  HFLSD_SC, que so sao preenchidos no ramo inercial (IND_FLOODROUTING=1), e de
  MINIMONT (so criado no modo inercial).
"""
from __future__ import annotations

import time

import numpy as np

from . import hydrology as hy
from . import io_input as io
from . import sediment as sed
from .fortran_io import caldat
from .objective import sed_fob_calib


# ---------------------------------------------------------------------------
# Gerador aleatorio (Numerical Recipes ran1) e semente
# ---------------------------------------------------------------------------
class Ran1:
    IA = 16807; IM = 2147483647; AM = 1. / 2147483647; IQ = 127773; IR = 2836
    NTAB = 32; NDIV = 1 + (2147483647 - 1) // 32; EPS = 1.2e-7; RNMX = 1. - 1.2e-7

    def __init__(self):
        self.iv = [0] * self.NTAB
        self.iy = 0

    def __call__(self, idum_box):
        """idum_box: lista [idum] (passagem por referencia)."""
        idum = idum_box[0]
        if idum <= 0 or self.iy == 0:
            idum = max(-idum, 1)
            for j in range(self.NTAB + 8, 0, -1):
                k = int(idum / self.IQ)
                idum = self.IA * (idum - k * self.IQ) - self.IR * k
                if idum < 0:
                    idum += self.IM
                if j <= self.NTAB:
                    self.iv[j - 1] = idum
            self.iy = self.iv[0]
        k = int(idum / self.IQ)
        idum = self.IA * (idum - k * self.IQ) - self.IR * k
        if idum < 0:
            idum += self.IM
        j = 1 + self.iy // self.NDIV
        self.iy = self.iv[j - 1]
        self.iv[j - 1] = idum
        idum_box[0] = idum
        return float(np.float32(min(np.float32(self.AM) * np.float32(self.iy), np.float32(self.RNMX))))


def semente():
    """SEMENTE: ISEED = TIME() (segundos desde 1970)."""
    return int(time.time())


# ---------------------------------------------------------------------------
# Hidrologia
# ---------------------------------------------------------------------------
def calibra(M):
    """CALIBRA.f90 - vazia na versao MGB-SED v1.0."""
    M.log("AVISO: a sub-rotina CALIBRA (calibracao hidrologica) e vazia no MGB-SED v1.0.")


def calibparam(M):
    """CalibParam.f90 - aplica multiplicadores PARX aos parametros 'OLD'."""
    M.WM[:] = M.WMOLD; M.B[:] = M.BOLD; M.KINS[:] = M.KIOLD; M.KBAS[:] = M.KBOLD
    M.PLAM[:] = M.PLAMOLD; M.CAP[:] = M.CAPOLD; M.WC[:] = M.WCOLD
    M.CS[:] = M.CSOLD; M.CI[:] = M.CIOLD; M.CB[:] = M.CBOLD
    cur = {1: M.WM, 2: M.B, 3: M.KBAS, 4: M.KINS, 5: M.PLAM, 6: M.CAP, 7: M.WC}
    old = {1: M.WMOLD, 2: M.BOLD, 3: M.KBOLD, 4: M.KIOLD, 5: M.PLAMOLD, 6: M.CAPOLD, 7: M.WCOLD}
    L = 0
    for i1 in range(M.NU):
        for i2 in range(1, 8):
            pc = M.p_Calib[i1, i2 - 1]
            if pc == 0:
                continue
            if pc == 1:
                L += 1
                cur[i2][:, i1] = old[i2][:, i1] * M.PARX[L - 1]
            elif pc < 0:
                i3 = -pc - 1
                cur[i2][:, i1] = old[i2][:, i1] * cur[i2][:, i3] / old[i2][:, i3]
    if M.p_Calib[M.NU, 0] == 1:
        L += 1; M.CS[:] = M.CSOLD * M.PARX[L - 1]
    if M.p_Calib[M.NU + 1, 0] == 1:
        L += 1; M.CI[:] = M.CIOLD * M.PARX[L - 1]
    if M.p_Calib[M.NU + 2, 0] == 1:
        L += 1; M.CB[:] = M.CBOLD * M.PARX[L - 1]
    for i1 in range(M.NCONGEL):
        ib = M.IBCONGEL[i1] - 1
        for k in cur:
            cur[k][ib, :] = old[k][ib, :]
        M.CS[ib] = M.CSOLD[ib]; M.CI[ib] = M.CIOLD[ib]; M.CB[ib] = M.CBOLD[ib]


# ---------------------------------------------------------------------------
# Sedimentos
# ---------------------------------------------------------------------------
def alloca_calib(M):
    R = M.R
    NS, NF, NPAR, NC, NT, NU = M.NS, M.NF, M.NPAR, M.NC, M.NT, M.NU
    M.IRUIM = np.zeros(NS, int)
    M.FO = np.zeros((NS, NF), R)
    M.IPARET = np.zeros(NS, int)
    M.XPAR = np.zeros(NPAR, R)
    M.SOMAPAR = np.zeros(NPAR, R); M.REFLEX = np.zeros(NPAR, R); M.CONTRA = np.zeros(NPAR, R)
    M.PPAR = np.zeros((NPAR, NS), R)
    M.PAR = np.zeros((NPAR, NS), R)
    M.PROB = np.zeros(NS, R)
    M.VMIN = np.zeros(NF, R)
    M.PARX = np.zeros(NPAR, R)
    M.DSUP_SC = np.zeros((NC, NT, NU), R)
    M.Q_SC = np.zeros((NC, NT), R)
    M.VOLTREC2_SC = np.zeros((NC, NT), R)
    M.VFL2_SC = np.zeros((NC, NT), R)
    M.SfTRECx_SC = np.zeros((NC, NT), R)
    M.HmJx_SC = np.zeros((NC, NT), R)
    M.HFLSD_SC = np.zeros((NC, NT), R)
    M.TKS_SC = np.zeros(NC, R)


def calibparam_sed(M):
    M.ALFsed[0] = M.PARX[0]; M.BETsed[0] = M.PARX[1]; M.AuxTKS[0] = M.PARX[2]
    for k in range(1, M.NB):
        M.ALFsed[k] = M.PARX[3 * k]; M.BETsed[k] = M.PARX[3 * k + 1]; M.AuxTKS[k] = M.PARX[3 * k + 2]


def _prepare_run(M):
    calibparam_sed(M)
    hy.condinic(M)
    io.sed_param(M)
    sed.sed_inicial(M)


def sed_calib_modelo(M):
    """SED_CALIB_MODELO.f90."""
    R = M.R
    M.IT = 0
    while M.IT < M.NT:
        M.IT += 1
        IT = M.IT; it = IT - 1
        M.JDIA = M.IDINI + int((IT + M.HORAINI - 1) / (86400. / M.DTP))
        M.IMES, M.IDIA, M.IANO = caldat(M.JDIA)
        M.DIAH[it] = M.IDIA; M.MESH[it] = M.IMES; M.ANOH[it] = M.IANO
        M.HORAH[it] = IT % 24 + M.HORAINI
        M.P[:] = M.rainfile.read(IT)
        M.PM[:] = M.PM + M.P
        M.TONTEM[:] = M.TA
        hy.leclima(M)
        sed.sed_bacia(M, calib_dsup=M.DSUP_SC[:, it, :], calib_tks=M.TKS_SC)
        sed_calib_rede(M)
        if M.SED_INDEX == 1:
            e6 = R(10.) ** R(6.)
            for kss in range(M.AUXSS):
                idx = M.ISSOBS[kss] - 1
                M.CSSCAL[kss][:, it] = (M.CSJ2[idx, 1] + M.CSJ2[idx, 2]) * e6 * M.QJ2[idx] * R(0.0864)


def sed_calib_rede(M):
    """SED_CALIB_REDE.f90 (usa as grandezas hidraulicas gravadas em *_SC)."""
    R = M.R
    NC = M.NC
    M.CSM1[:] = M.CSM2; M.CSJ1[:] = M.CSJ2; M.CFL1[:] = M.CFL2
    M.CARGM[:] = 0.0; M.CSJ2[:] = 0.0; M.CFL2[:] = 0.0
    Qlim = R(0.01)
    it = M.IT - 1
    act = M.sed_active
    Hm = M.zeros(NC); Sf = M.zeros(NC); Vm = M.zeros(NC); Uat = M.zeros(NC)
    QtFL = M.zeros(NC); HFLSD = M.zeros(NC)

    def upsum(col):
        s = np.zeros(NC, R)
        for ic in range(NC):
            ne = M.MINIMONT[ic, 0]
            if ne > 0:
                s[ic] = col[M.MINIMONT[ic, 1:1 + ne] - 1].sum()
        return s

    if M.IT == 1:
        M.VolTREC1[act] = 0.0; M.VFL1[act] = 0.0; M.QM1[:NC][act] = 0.0; M.QJ1[act] = 0.0
    else:
        M.VolTREC1[act] = M.VOLTREC2_SC[act, it - 1]
        M.VFL1[act] = M.VFL2_SC[act, it - 1]
        M.QM1[:NC][act] = upsum(M.Q_SC[:, it - 1])[act]
        M.QJ1[act] = M.Q_SC[act, it - 1]
    M.QJ2[act] = M.Q_SC[act, it]
    M.QM2[:NC][act] = upsum(M.Q_SC[:, it])[act]
    # (sic) VOLTREC2=VOLTREC2_SC(IC,IT) atribui um escalar ao vetor inteiro;
    # dentro do laco por minibacia, o valor efetivo para IC e o proprio VOLTREC2_SC(IC,IT)
    M.VolTREC2[act] = M.VOLTREC2_SC[act, it]
    M.VFL2[act] = M.VFL2_SC[act, it]
    Sf[act] = M.SfTRECx_SC[act, it]
    Hm[act] = M.HmJx_SC[act, it]
    HFLSD[act] = M.HFLSD_SC[act, it]
    QtFL[act] = (M.VFL2[act] - M.VFL1[act]) / R(M.DTP)
    for cells in M.levels:
        cells = cells[act[cells]]
        if cells.size:
            sed._sed_propag(M, cells, Hm[cells], Vm[cells], Sf[cells], Uat[cells], Qlim,
                            QtFL[cells], HFLSD[cells])


def calibra_sed(M, seed=None):
    """CALIBRA_SED (MOCOM-UA multiobjetivo)."""
    import os
    from .simulation import modelo
    R = M.R
    NS, NPAR, NF, NB = M.NS, M.NPAR, M.NF, M.NB
    NPCOMPLEX = 5
    alloca_calib(M)
    M.CSSCAL = [M.zeros(int(M.NOBSS[k]), M.NT) for k in range(M.AUXSS)]
    fevo = open(os.path.join(M.OUTPUT_DIRECTORY, "EVOLUTION.TXT"), "w")
    fres = open(os.path.join(M.OUTPUT_DIRECTORY, "RESULCALIB.TXT"), "w")
    ran1 = Ran1()
    iseed = [semente() if seed is None else int(seed)]
    for _ in range(iseed[0] % 100):
        ran1(iseed)
    M.PAR[0, 0] = 11.8; M.PAR[1, 0] = 0.56; M.PAR[2, 0] = 1.0
    for i in range(1, NB):
        M.PAR[3 * i, 0] = 11.8; M.PAR[3 * i + 1, 0] = 0.56; M.PAR[3 * i + 2, 0] = 1.00
    for i in range(1, NS):
        for L in range(NPAR):
            M.PPAR[L, i] = ran1(iseed)
            M.PAR[L, i] = M.PMIN[L] + M.PPAR[L, i] * (M.PMAX[L] - M.PMIN[L])
    M.xyz = 1
    M.PARX[:] = M.PAR[:, 0]
    _prepare_run(M)
    modelo(M, writers={}, progress=True)
    KOUNTF = 0
    M.log("Gera populacao inicial:")
    for i in range(NS):
        M.PARX[:] = M.PAR[:, i]
        _prepare_run(M)
        sed_calib_modelo(M)
        # (sic) o Fortran nao chama SED_Fob_Calib nem atribui FO(I,:) aqui
    M.log("AVISO: FO da populacao inicial nao e avaliada no codigo original (FO=0).")
    ISHUFFLE = 0
    IRMAX = 10
    while IRMAX > 1 and ISHUFFLE < M.iMaxGen:
        ISHUFFLE += 1
        M.IRUIM[:] = 0
        M.IPARET[:] = 1
        ISOMA = 0
        IRMAX = 1
        IDOMIN = 0
        while True:                                     # label 3000
            for I in range(NS):
                if M.IPARET[I] == IRMAX:
                    for J in range(NS):
                        if J != I and M.IPARET[J] == IRMAX:
                            dominated = True
                            for K in range(NF):
                                if M.FO[I, K] > M.FO[J, K]:
                                    IDOMIN = 1
                                else:
                                    IDOMIN = 0
                                    dominated = False
                                    break               # GOTO 1000
                            if dominated and IDOMIN == 1:
                                break                   # GOTO 2000
                    if IDOMIN == 1:
                        M.IPARET[I] += 1
                        ISOMA = 1
                    else:
                        M.IPARET[I] = IRMAX
                    IDOMIN = 0
            if ISOMA == 1:
                IRMAX += 1
                ISOMA = 0
                continue
            break
        if IRMAX == 1:
            break
        fevo.write("\n " + str(ISHUFFLE) + "\n")
        for I in range(NS):
            fevo.write("".join(f"{v:10.3f}" for v in M.PAR[:, I]) +
                       "".join(f"{v:12.4f}" for v in M.FO[I]) + f"{M.IPARET[I]:8d}\n")
        NPC = 0
        for I in range(NS):
            if M.IPARET[I] == IRMAX:
                M.IRUIM[NPC] = I
                NPC += 1
        SPARET = R(M.IPARET.sum())
        RMAX = R(IRMAX)
        for I in range(NS):
            M.PROB[I] = (RMAX - M.IPARET[I] + R(1.)) / (NS * (RMAX + R(1.)) - SPARET)
        for I in range(1, NS):
            M.PROB[I] = M.PROB[I - 1] + M.PROB[I]
        if abs(M.PROB[NS - 1] - 1.0) > 0.001:
            raise RuntimeError("ALGO ERRADO PROBABILID")
        M.PROB[NS - 1] = 1.0
        SPAR = np.zeros((NPC, NPCOMPLEX, NPAR), R)
        FPLEX = np.zeros((NPC, NPCOMPLEX, NF), R)
        for K in range(NPC):
            LRUIM = M.IRUIM[K]
            SPAR[K, 0, :] = M.PAR[:, LRUIM]
            for J in range(1, NPCOMPLEX):
                XRAND = ran1(iseed)
                PINI = 0.0
                for I in range(NS):
                    if XRAND >= PINI and XRAND <= M.PROB[I]:
                        SPAR[K, J, :] = M.PAR[:, I]
                        FPLEX[K, J, :] = M.FO[I, :]
                    PINI = M.PROB[I]
        for K in range(NPC):
            LRUIM = M.IRUIM[K]
            IREJECT = 0
            SOMAPAR = SPAR[K, 1:, :].sum(axis=0) / (NPCOMPLEX - 1)
            DIFPAR = SOMAPAR - SPAR[K, 0, :]
            REFLEX = SOMAPAR + DIFPAR
            CONTRA = SOMAPAR - DIFPAR / R(2.)
            for L in range(NPAR):
                if REFLEX[L] < M.PMIN[L] or REFLEX[L] >= M.PMAX[L]:
                    M.log(" reflex rejeitado parametro", L + 1)
                    IREJECT = 1
                    break
            if IREJECT == 0:
                M.PARX[:] = REFLEX
                _prepare_run(M)
                sed_calib_modelo(M)
                sed_fob_calib(M)
                KOUNTF += 1
                M.FO[LRUIM, :] = M.VFO[:NF]
                for J in range(1, NPCOMPLEX):
                    nworse = int((FPLEX[K, J, :] < M.FO[LRUIM, :]).sum())
                    if nworse == NF:
                        M.log("PONTO REJEITADO POR", FPLEX[K, J, :])
                        IREJECT = 1
                        break
            if IREJECT == 0:
                M.log("PONTO ACEITO")
                M.PAR[:, LRUIM] = REFLEX
            else:
                M.PARX[:] = CONTRA
                _prepare_run(M)
                sed_calib_modelo(M)
                sed_fob_calib(M)
                KOUNTF += 1
                M.FO[LRUIM, :] = M.VFO[:NF]
                M.PAR[:, LRUIM] = CONTRA
        M.VMIN[:] = M.FO.min(axis=0)
        M.log(f"{ISHUFFLE:8d}{IRMAX:8d}{NPC:8d}" + "".join(f"{v:12.4f}" for v in M.VMIN))
    fevo.close()
    fres.close()
