"""
Leitura dos arquivos de entrada.

Traducao de: LEFIX.f90, ALLOCA_VARS.f90 (parte), LEVAR.f90, LEUSO.f90,
LECELL.f90, LECLIMED.f90, ARQCLISUB.f90, LEQOBS.f90, LESUBST.f90,
LEFLOODPLN.f90, SED_LEObs.f90 (LEOBSed), SED_PARAM.F90, LeCalib.f90,
SED_LeCalib.f90 (LeCalib_SED), LECHUVA.f90, leprev.f90.
"""
from __future__ import annotations

import os

import numpy as np

from .fortran_io import FortranFile, read_fixed
from . import fmath as fm


def _path(M, name):
    """Resolve nome de arquivo sem diferenciar maiusculas/minusculas
    (o Windows nao diferencia; o Linux sim)."""
    p = os.path.join(M.INPUT_DIRECTORY, name)
    if os.path.exists(p):
        return p
    low = name.lower()
    for f in os.listdir(M.INPUT_DIRECTORY):
        if f.lower() == low:
            return os.path.join(M.INPUT_DIRECTORY, f)
    raise FileNotFoundError(p)


# ---------------------------------------------------------------------------
# LEFIX - infoMGB.sim
# ---------------------------------------------------------------------------
def lefix(M):
    f = FortranFile(_path(M, "infoMGB.sim"))
    A20 = [("A", 20)]
    I10x8 = [("I", 10, 8)]
    f.skip(3)                    # READ 75 x3
    f.skip(1)                    # READ(*) TITULO
    f.skip(1)                    # READ 72 CABE
    M.IDIA, M.IMES, M.IANO, M.HORAINI = f.read_fixed(I10x8, 4)
    f.skip(1)
    f.skip(1)
    nt, dtp, alpha = f.read_fixed([("I", 10), ("F", 10, 1), ("F", 10, 2)])
    M.NT, M.DTP, M.alpha = nt, float(np.float32(dtp)), float(np.float32(alpha))
    f.skip(1)
    f.skip(1)
    (M.NC, M.NU, M.NB, M.NCLI, M.CRU_INDEX, M.RESERV_INDEX,
     M.IND_FLOODROUTING, M.SED_INDEX) = f.read_fixed(I10x8, 8)
    M.log("SED_INDEX", M.SED_INDEX)
    f.skip(1)
    f.skip(1)
    M.ICALIB, M.AUXCALIBSS = f.read_fixed(I10x8, 2)
    f.skip(1)
    f.skip(1)
    M.flag_TC, M.FLAG_SEDVAR = f.read_fixed(I10x8, 2)
    f.skip(1)
    f.skip(1)                    # READ 75 TITULO
    M.ARQCLI = []
    if M.NCLI <= 0:
        M.NCLI = -M.NCLI
        M.flagaclimed = 1        # usa climatologia (CRU)
    else:
        M.flagaclimed = 0
        for _ in range(M.NCLI):
            M.ARQCLI.append(f.read_fixed(A20)[0].strip())
    f.skip(1)
    f.skip(1)
    M.ACLIMED = f.read_list(1, [str])[0]
    f.skip(1)
    f.skip(1)
    M.NOBS, M.ARQOBS = f.read_list(2, [int, str])
    f.skip(1)
    f.skip(1)
    M.AUXSS = f.read_list(1, [int])[0]
    M.NOBSS = np.zeros(M.AUXSS, int)
    M.NUMSS = np.zeros(M.AUXSS, int)
    M.ARQOBSSS = []
    for k in range(M.AUXSS):
        n, nm = f.read_list(2, [int, str])
        M.NOBSS[k] = n
        M.ARQOBSSS.append(nm)
        M.log("NUMBER OF DATA RELATED TO SEDIMENTS AND FILE NAME: ", n, nm)
    f.skip(1)
    f.skip(1)
    M.IQOBS = np.array(f.read_list(M.NOBS, [int]), int)          # codigos 1-based
    f.skip(1)
    f.skip(1)
    M.ISSOBS = [np.array(f.read_list(int(M.NOBSS[k]), [int]), int) for k in range(M.AUXSS)]
    f.skip(1)
    f.skip(1)
    M.NUMHIDG = f.read_list(1, [int])[0]
    f.skip(1)
    f.skip(1)
    for k in range(M.AUXSS):
        M.NUMSS[k] = f.read_list(1, [int])[0]
    f.skip(1)
    f.skip(1)
    M.IHIDG = np.array([f.read_list(1, [int])[0] for _ in range(M.NUMHIDG)], int)
    M.ISS = []
    for k in range(M.AUXSS):
        f.skip(1)                # READ 75 TITULO
        M.ISS.append(np.array(f.read_list(int(M.NUMSS[k]), [int]), int))
    f.skip(1)
    f.skip(1)
    M.NUMSUBST, M.ARQSUBST = f.read_list(2, [int, str])
    f.skip(1)
    f.skip(1)
    M.ISUBST = np.array([f.read_list(1, [int])[0] for _ in range(M.NUMSUBST)], int)
    f.skip(1)
    f.skip(1)
    M.SUBini, M.SUBfim = f.read_list(2, [int, int])


# ---------------------------------------------------------------------------
# ALLOCA_VARS(0) - alocacao das variaveis principais usadas na simulacao
# ---------------------------------------------------------------------------
def alloca_vars(M):
    R = M.R
    NC, NU, NB, NT, NCLI = M.NC, M.NU, M.NB, M.NT, M.NCLI
    z = M.zeros
    M.DIAH = np.zeros(NT, int); M.MESH = np.zeros(NT, int)
    M.ANOH = np.zeros(NT, int); M.HORAH = np.zeros(NT, int)
    M.EVQ = z(NC)
    M.QRIOINI = z(NC, M.NUMUSK + 1)
    M.QCONTORM = z(NC + 1, 25)
    M.QCONTORJ = z(NC + 1, 25)
    M.SI = z(NC, NU)
    M.XYC = z(NCLI, 2)
    M.TAMM = z(NCLI, 12); M.URMM = z(NCLI, 12); M.VVMM = z(NCLI, 12)
    M.PAMM = z(NCLI, 12); M.SOLMM = z(NCLI, 12)
    M.PMB = z(NB); M.KPM = np.zeros(NB, int)
    M.QREF = z(NC)
    M.ICBOM = np.zeros(NC, int)
    M.PUSO = z(NC, NU)
    M.ACEL = z(NC); M.ACUR = z(NC); M.SRIO = z(NC); M.DECL = z(NC)
    M.X = z(NC); M.Y = z(NC)
    M.IBAC = np.zeros(NC, int); M.HCEL = z(NC); M.LCEL = z(NC)
    M.CELJUS = np.zeros(NC, int)
    M.PM = z(NC)
    M.BRIO = z(NC)
    M.IEXUT = np.full(NB, -999999, int)
    M.CI = z(NB); M.CB = z(NB); M.CS = z(NB)
    M.PLAM = z(NB, NU); M.CAP = z(NB, NU); M.WC = z(NB, NU)
    M.WM = z(NB, NU); M.B = z(NB, NU); M.KINS = z(NB, NU); M.KBAS = z(NB, NU)
    M.NSUBT = np.zeros(NC, int)
    M.DT = z(NC)
    M.CEL = z(NC); M.TIND = z(NC)
    M.ALB = z(NU, 12); M.RIAF = z(NU, 12); M.Z = z(NU, 12); M.RS = z(NU, 12)
    M.QESP = z(NB)
    M.QB = z(NT); M.QBI = z(NT); M.QBIS = z(NT)
    M.DSUM = z(NB, NT); M.DINM = z(NB, NT); M.DBAM = z(NB, NT)
    M.SIM = z(NB, NT); M.EM = z(NB, NT); M.WBM = z(NB, NT); M.PM2 = z(NB, NT)
    M.KCB = np.zeros(NB, int)
    M.QR = z(M.NOBS, NT); M.QOBS = z(M.NOBS, NT)
    M.QLIDO = z(max(M.NUMSUBST, 0), NT)
    M.QRG = z(M.NUMHIDG, NT)
    M.QRB = z(NB, NT)
    M.QM1 = z(NC + 1); M.QJ1 = z(NC); M.QM2 = z(NC + 1); M.QJ2 = z(NC)
    M.QCEL1 = z(NC); M.QCEL2 = z(NC)
    M.PMB2 = z(NC + 1); M.PMI2 = z(NC + 1); M.PMS2 = z(NC + 1)
    M.PJB2 = z(NC + 1); M.PJI2 = z(NC + 1); M.PJS2 = z(NC + 1)
    M.VRB = z(NC); M.VRI = z(NC); M.VRS = z(NC)
    M.ET = z(NC, NU); M.CAF = z(NC, NU); M.W = z(NC, NU)
    M.QBAS = z(NC); M.QINT = z(NC); M.QSUP = z(NC)
    M.VBAS = z(NC); M.VINT = z(NC); M.VSUP = z(NC)
    M.TONTEM = z(NC); M.P = z(NC)
    M.TA = z(NC); M.UR = z(NC); M.VV = z(NC); M.SOL = z(NC); M.PA = z(NC)
    M.R2 = z(M.NOBS); M.ERRV = z(M.NOBS); M.R2L = z(M.NOBS)
    M.VFO = z(M.NFP)
    M.OD = np.zeros(NC, int)
    M.sbtFLAG = np.zeros(NC, int)      # alocado e nunca preenchido no Fortran (=0)
    M.ATIVABACIA = np.zeros(NB, int)
    M.RUGMAN = z(NC)
    M.E0agua = z(NC); M.E0TOPO = z(NC); M.E0SUP = z(NC)
    # Campos (NC,NT) do Fortran: QTUDO, QITUDO, QBASTUDO, EVAPTUDO, WTUDO.
    # So sao lidos no modo de previsao; para economizar memoria sao
    # alocados sob demanda (ver simulation.alloc_full_fields).
    M.QTUDO = None; M.QITUDO = None; M.QBASTUDO = None
    M.EVAPTUDO = None; M.WTUDO = None
    # acumuladores escalares usados em CELULA
    M.CALIBFLAG_SET = False


def alloca_varsinerc(M):
    """ALLOCA_VARSINERC(0)."""
    NC = M.NC
    M.HDFLAG = np.zeros(NC, int)
    M.MINIMONT = np.zeros((NC, 10), int)
    M.NPFL = np.zeros(NC, int)
    M.ZFUNDOFL = np.zeros(NC, int)      # INTEGER no Fortran (sic)
    M.ZFL = np.zeros((100, NC))
    M.AFL = np.zeros((100, NC))
    M.HRIO = np.zeros(NC)
    M.HRIO1 = np.zeros(NC); M.HRIO2 = np.zeros(NC)
    M.ZTAB = np.zeros((102, NC)); M.VTAB = np.zeros((102, NC)); M.ATAB = np.zeros((102, NC))
    M.Q2fl = np.zeros(NC); M.Vel2fl = np.zeros(NC)
    M.Vol2 = np.zeros(NC); M.Vol1 = np.zeros(NC)
    M.Area2 = np.zeros(NC)
    M.Hfl = np.zeros(NC); M.Yfl = np.zeros(NC)
    M.jtab = np.zeros(NC, int)
    M.nMan = np.zeros(NC)
    M.HRG = M.zeros(M.NUMHIDG, M.NT)
    M.AFLTUDO = M.zeros(1, M.NT)
    M.YTUDO = None


# ---------------------------------------------------------------------------
# LEVAR - ALBIAF.fix
# ---------------------------------------------------------------------------
def levar(M):
    f = FortranFile(_path(M, "ALBIAF.fix"))
    fmt72 = [("A", 10), ("F", 5, 2, 12)]
    M.AUSO = [""] * M.NU
    for arr in (M.ALB, M.RIAF, M.Z, M.RS):
        f.skip(2)
        for iu in range(M.NU):
            v = f.read_fixed(fmt72)
            M.AUSO[iu] = v[0]
            arr[iu, :] = v[1:]


# ---------------------------------------------------------------------------
# LEUSO - PARUSO.cal
# ---------------------------------------------------------------------------
def leuso(M):
    f = FortranFile(_path(M, "paruso.cal"))
    fmt71 = [("A", 10), ("F", 8, 3, 7)]
    for kb in range(M.NB):
        f.skip(2)
        for iu in range(M.NU):
            v = f.read_fixed(fmt71)
            M.AUSO[iu] = v[0]
            (M.WM[kb, iu], M.B[kb, iu], M.KBAS[kb, iu], M.KINS[kb, iu],
             M.PLAM[kb, iu], M.CAP[kb, iu], M.WC[kb, iu]) = v[1:]
        M.CS[kb] = f.read_fixed([("A", 10), ("F", 8, 1)])[1]
        M.CI[kb] = f.read_fixed([("A", 10), ("F", 8, 1)])[1]
        M.CB[kb] = f.read_fixed([("A", 10), ("F", 8, 1)])[1]
        M.QESP[kb] = f.read_fixed([("A", 10), ("F", 8, 6)])[1]


# ---------------------------------------------------------------------------
# LECELL - MINI.gtp
# ---------------------------------------------------------------------------
def lecell(M):
    f = FortranFile(_path(M, "MINI.gtp"))
    f.skip(1)
    NU = M.NU
    for i in range(M.NC):
        v = f.read_list(17 + NU)
        icell = int(v[1])
        if icell != i + 1:
            raise ValueError("MINI.gtp: minibacias devem estar em ordem 1..NC")
        M.X[i] = float(v[2]); M.Y[i] = float(v[3])
        M.IBAC[i] = int(v[4]) - 1                    # 0-based
        M.ACEL[i] = float(v[5]); M.ACUR[i] = float(v[6])
        M.SRIO[i] = float(v[7]); M.DECL[i] = float(v[8])
        M.LCEL[i] = float(v[9]); M.HCEL[i] = float(v[10])
        M.CELJUS[i] = int(v[11])                     # 1-based; -1 = exutorio
        M.OD[i] = int(v[12]); M.HDFLAG[i] = int(v[13])
        M.BRIO[i] = float(v[14]); M.HRIO[i] = float(v[15])
        M.RUGMAN[i] = float(v[16])
        M.PUSO[i, :] = [float(x) for x in v[17:17 + NU]]
        s = M.PUSO[i, :].sum()
        if s > 100.01 or s < 99.99 or M.PUSO[i].min() < 0 or M.PUSO[i].max() > 100:
            M.log("Error in the % of HRUs in the catchments", i + 1)
        kw = M.IBAC[i]
        if M.IEXUT[kw] < icell:
            M.IEXUT[kw] = icell                       # 1-based (como no Fortran)
    M.QREF[:] = 30 * M.ACUR / 1000.
    M.nMan[:] = M.RUGMAN.astype(np.float64)
    M.HRIO1[:] = M.HRIO
    M.HRIO2[:] = M.HRIO
    M.DECL[:] = M.DECL / 1000.0
    M.LCEL[:] = np.maximum(M.LCEL, M.R(0.001))
    # indice 0-based da minibacia de jusante (-1 = nenhuma)
    M.jus = np.where(M.CELJUS > 0, M.CELJUS - 1, -1)


# ---------------------------------------------------------------------------
# LECLIMED - medias.cli
# ---------------------------------------------------------------------------
def leclimed(M):
    f = FortranFile(_path(M, "medias.cli"))
    for ip in range(M.NCLI):
        v = f.read_list(2, [float, float])
        M.XYC[ip, :] = v
    fmt713 = [("A", 20), ("F", 10, 2, 12)]
    for arr in (M.TAMM, M.URMM, M.SOLMM, M.VVMM, M.PAMM):
        f.skip(2)
        for ip in range(M.NCLI):
            arr[ip, :] = f.read_fixed(fmt713)[1:]


# ---------------------------------------------------------------------------
# ARQCLISUB - dados climaticos diarios
# ---------------------------------------------------------------------------
def arqclisub(M):
    NB = M.NB
    M.KCB[0] = M.IEXUT[0]
    for kb in range(1, NB):
        M.KCB[kb] = M.IEXUT[kb] - M.IEXUT[kb - 1]
    shape = (M.NCLI, M.NT)
    if M.flagaclimed == 1:
        # Fortran preenche TD..PAD com -9999 (NCLI x NT). Nao sao usados.
        M.TD = M.UD = M.VD = M.SOLD = M.PAD = None
    else:
        M.TD = M.zeros(*shape); M.UD = M.zeros(*shape); M.VD = M.zeros(*shape)
        M.SOLD = M.zeros(*shape); M.PAD = M.zeros(*shape)
        fmt = [("F", 10, 2, 5)]
        for kli, nome in enumerate(M.ARQCLI):
            f = FortranFile(_path(M, nome))
            f.skip(1)
            for it in range(M.NT):
                v = f.read_fixed(fmt)
                M.TD[kli, it], M.UD[kli, it], M.VD[kli, it], M.SOLD[kli, it], M.PAD[kli, it] = v


# ---------------------------------------------------------------------------
# LEQOBS / LESUBST / LEOBSed
# ---------------------------------------------------------------------------
def _read_series(path, nseries, nt, dtype):
    out = np.zeros((nseries, nt), dtype=dtype)
    f = FortranFile(path)
    f.skip(1)
    for it in range(nt):
        v = f.read_list(3 + nseries)
        out[:, it] = [float(x) for x in v[3:3 + nseries]]
    return out


def leqobs(M):
    M.QOBS[:, :] = _read_series(_path(M, M.ARQOBS), M.NOBS, M.NT, M.R)


def lesubst(M):
    M.QLIDO = _read_series(_path(M, M.ARQSUBST), M.NUMSUBST, M.NT, M.R)


def leobsed(M):
    """LEOBSed: le dados observados de sedimentos (SSOBS)."""
    M.SSOBS = []
    for kss in range(M.AUXSS):
        M.SSOBS.append(_read_series(_path(M, M.ARQOBSSS[kss]), int(M.NOBSS[kss]), M.NT, M.R))


# ---------------------------------------------------------------------------
# LEFLOODPLN - COTA_AREA.FLP (somente modelo inercial)
# ---------------------------------------------------------------------------
def lefloodpln(M):
    f = FortranFile(_path(M, "COTA_AREA.FLP"))
    f.skip(1)
    M.NPFL[:] = 0
    while not f.eof():
        rec = f.next_record()
        if rec.strip() == "":
            continue
        v = rec.split()
        ic = int(v[0]) - 1
        M.ZFUNDOFL[ic] = int(float(v[1])) if "." not in v[1] else int(float(v[1]))
        k = M.NPFL[ic]
        M.ZFL[k, ic] = float(v[2])
        M.AFL[k, ic] = float(v[3])
        M.NPFL[ic] += 1


# ---------------------------------------------------------------------------
# Chuva: CHUVABIN.PBI (acesso direto, RECL=NC palavras de 4 bytes)
# ---------------------------------------------------------------------------
class RainFile:
    def __init__(self, path, nc):
        self.mm = np.memmap(path, dtype="<f4", mode="r")
        self.nc = nc
        self.nrec = self.mm.size // nc

    def read(self, rec1: int) -> np.ndarray:
        """LECHUVA: READ(FILPLU,REC=ITCHUVA). rec1 e 1-based."""
        i0 = (rec1 - 1) * self.nc
        return np.array(self.mm[i0:i0 + self.nc])


# ---------------------------------------------------------------------------
# SED_PARAM - parametros da MUSLE e arquivos do pre-processamento
# ---------------------------------------------------------------------------
def sed_param(M):
    NB, NU, NC = M.NB, M.NU, M.NC
    R = M.R
    if M.xyz == 1:
        M.Kusle = M.zeros(NB, NU); M.Cusle = M.zeros(NB, NU); M.Pusle = M.zeros(NB, NU)
        M.Rgros = M.zeros(NB, NU); M.Ksdr = M.zeros(NB, NU)
        M.Mareia = M.zeros(NB, NU); M.Msilte = M.zeros(NB, NU); M.Margila = M.zeros(NB, NU)
        M.Morg = M.zeros(NB, NU); M.Mrocha = M.zeros(NB, NU)
        M.LSAcu = M.zeros(NC, NU); M.NPXU = np.zeros((NC, NU), int); M.SDR = M.zeros(NC, NU + 1)
        M.PSC = np.zeros((NC, 4)); M.SLC = np.zeros(NC); M.QSC = np.zeros(NC)
        M.QSAUX = np.zeros((NC, NU + 2))                     # REAL(8)
        M.SLU = M.zeros(NC, NU)
        M.QSU = M.zeros(NC, NU); M.PSU = M.zeros(NC, NU); M.SLXU = M.zeros(NU)
    M.Kusle[:] = 0; M.Cusle[:] = 0; M.Pusle[:] = 0; M.Rgros[:] = 0; M.Ksdr[:] = 0
    M.PSC[:] = 0; M.SLC[:] = 0; M.QSC[:] = 0; M.QSAUX[:] = 0; M.SLU[:] = 0
    M.QSU[:] = 0; M.PSU[:] = 0; M.SLXU[:] = 0
    M.log("DETERMINANDO OS PARAMETROS DA MUSLE...")
    fp = FortranFile(_path(M, "PARUSO_MUSLE.TXT"))
    ft = FortranFile(_path(M, "PARTEXT_MUSLE.TXT"))
    fmt73 = [("A", 10), ("F", 10, 5, 5)]
    fmt74 = [("A", 10), ("F", 10, 3, 5)]
    M.Nuso = [""] * NU
    f32 = R
    for ib in range(NB):
        fp.skip(2); ft.skip(2)
        for iu in range(NU):
            v = fp.read_fixed(fmt73)
            M.Nuso[iu] = v[0]
            M.Kusle[ib, iu], M.Cusle[ib, iu], M.Pusle[ib, iu], M.Rgros[ib, iu], M.Ksdr[ib, iu] = v[1:]
            v = ft.read_fixed(fmt74)
            M.Mareia[ib, iu], M.Msilte[ib, iu], M.Margila[ib, iu], M.Morg[ib, iu], M.Mrocha[ib, iu] = v[1:]
            if M.Kusle[ib, iu] == -1.:
                np.seterr(invalid='ignore', divide='ignore')   # Fortran gera NaN p/ URH sem textura (agua)
                ma = M.Mareia[ib, iu]; ms = M.Msilte[ib, iu]; mg = M.Margila[ib, iu]; mo = M.Morg[ib, iu]
                one = f32(1.)
                Fagros = f32(0.2) + f32(0.3) * fm.exp(f32(-0.0256) * ma * (one - ms / f32(100.)))
                Fargila = fm.pw(ms / (mg + ms), f32(0.3))
                Forg = one - (f32(0.25) * mo / (mo + fm.exp(f32(3.72) - f32(2.95) * mo)))
                xa = one - ma / f32(100.)
                Fareia = one - (f32(0.7) * xa / (xa + fm.exp(f32(-5.51) + f32(22.9) * xa)))
                M.Kusle[ib, iu] = Fagros * Fargila * Forg * Fareia
            if M.Rgros[ib, iu] == -1:
                M.Rgros[ib, iu] = fm.exp(f32(-0.053) * M.Mrocha[ib, iu])
    M.log("LENDO ARQUIVOS DO PRE_MUSLE...")
    fl = FortranFile(_path(M, "SED_LSm.txt"))
    fh = FortranFile(_path(M, "SED_HRU.txt"))
    fs = FortranFile(_path(M, "SED_SDR.txt"))
    for ic in range(NC):
        M.LSAcu[ic, :] = fl.read_fixed([("F", 12, 4, NU)])
        M.NPXU[ic, :] = fh.read_fixed([("I", 8, NU)])
        M.SDR[ic, :] = fs.read_fixed([("F", 12, 4, NU + 1)])
    M.log("CALCULO DOS PARAMETROS DA MUSLE ENCERRADO!!!")


# ---------------------------------------------------------------------------
# LeCalib (ParCalib.aut) e LeCalib_SED (ParCalibSED.txt)
# ---------------------------------------------------------------------------
def lecalib(M):
    f = FortranFile(_path(M, "ParCalib.aut"))
    f.skip(2)
    M.NS = f.read_list(1, [int])[0]
    f.skip(1)
    M.NF = f.read_list(1, [int])[0]
    f.skip(4)
    M.iFO = np.array(f.read_list(M.NF, [int]), int)
    f.skip(1)
    M.iMaxGen = f.read_list(1, [int])[0]
    f.skip(1)
    M.NCONGEL = f.read_list(1, [int])[0]
    f.skip(1)
    M.IBCONGEL = np.array(f.read_list(M.NCONGEL, [int]), int) if M.NCONGEL > 0 else np.zeros(0, int)
    if M.NCONGEL == 0:
        f.skip(1)
    M.p_Calib = np.zeros((M.NU + 3, 7), int)
    f.skip(5)
    for i in range(M.NU):
        v = f.read_fixed([("A", 10), ("I", 6, 7)])
        M.p_Calib[i, :] = v[1:]
    for i in range(M.NU, M.NU + 3):
        v = f.read_fixed([("A", 10), ("I", 6)])
        M.p_Calib[i, 0] = v[1]
    M.NPAR = int((M.p_Calib == 1).sum())
    f.skip(2)
    M.PMIN = M.zeros(M.NPAR); M.PMAX = M.zeros(M.NPAR)
    for i in range(M.NPAR):
        v = f.read_list(3)
        M.PMIN[i], M.PMAX[i] = float(v[0]), float(v[1])
    M.CalibFLAG = M.zeros(M.NOBS, M.NF)
    f.skip(1)
    for i in range(M.NOBS):
        M.CalibFLAG[i, :] = [float(x) for x in f.read_list(M.NF)]


def lecalib_sed(M):
    f = FortranFile(_path(M, "ParCalibSED.txt"))
    f.skip(1)
    M.NS = f.read_list(1, [int])[0]
    f.skip(1)
    M.NF = f.read_list(1, [int])[0]
    f.skip(1)
    M.iMaxGen = f.read_list(1, [int])[0]
    f.skip(1)
    M.NPAR = f.read_list(1, [int])[0]
    f.skip(1)
    M.PMIN = M.zeros(M.NPAR); M.PMAX = M.zeros(M.NPAR)
    for i in range(M.NPAR):
        M.PMIN[i], M.PMAX[i] = [float(x) for x in f.read_list(2)]
    f.skip(1)
    M.CalibFLAG = M.zeros(M.NBP, M.NF)
    M.CalibFLAG[0, :] = [float(x) for x in f.read_list(M.NF)]
    kb = int(M.NOBSS[M.AUXCALIBSS - 1])
    for i in range(1, kb):
        M.CalibFLAG[i, :] = M.CalibFLAG[0, :]


# ---------------------------------------------------------------------------
# PAR_MUSLE_SED.txt (novo; substitui os valores escritos em SED_BACIA.F90)
# ---------------------------------------------------------------------------
def le_par_musle_sed(M):
    """Le alfa, beta e AuxTKS da MUSLE por sub-bacia, se o arquivo existir.

    Formato (texto livre; linhas iniciadas por '!' sao comentarios):
        ! IB    ALF    BET   AuxTKS
           1  11.80   0.56   1.00
           ...  (uma linha por sub-bacia, IB = 1..NB)
    Sem o arquivo, SED_BACIA usa os valores padrao (sediment.ALF_PADRAO etc.).
    """
    try:
        p = _path(M, "PAR_MUSLE_SED.txt")
    except FileNotFoundError:
        M.PAR_MUSLE_SED = None
        M.log("PAR_MUSLE_SED.txt ausente: MUSLE com valores padrao (ALF=11.8, BET=0.56, AuxTKS=1.0)")
        return
    alf = np.full(M.NB, np.nan); bet = np.full(M.NB, np.nan); atks = np.full(M.NB, np.nan)
    with open(p, encoding="latin-1") as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("!"):
                continue
            t = ln.split()
            ib = int(t[0]) - 1
            if not 0 <= ib < M.NB:
                raise ValueError(f"PAR_MUSLE_SED.txt: sub-bacia {ib + 1} fora de 1..{M.NB}")
            alf[ib], bet[ib], atks[ib] = float(t[1]), float(t[2]), float(t[3])
    faltam = np.nonzero(np.isnan(alf))[0] + 1
    if faltam.size:
        raise ValueError(f"PAR_MUSLE_SED.txt: faltam as sub-bacias {faltam.tolist()}")
    M.PAR_MUSLE_SED = (alf, bet, atks)
    M.log("MUSLE: alfa, beta e AuxTKS lidos de PAR_MUSLE_SED.txt")
