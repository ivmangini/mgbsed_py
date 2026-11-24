"""
Hidrologia das minibacias (balanco vertical) e pre-calculos.

Traducao de: PARCEL.f90, INTECLIM.f90, REGION.f90, PARCUNGE.f90, CondInic.f90,
leclima.F90, RADIACAO.f90, EVAPO.f90, SOLO.f90, CELULA.f90.

Estrategia de vetorizacao
-------------------------
No Fortran, CELULA percorre minibacias (IC) e URHs (IU) em lacos aninhados,
chamando RADIACAO, EVAPO e SOLO para cada par (IC,IU). Esses calculos sao
independentes entre pares e foram vetorizados em matrizes (NC,NU). As
acumulacoes por minibacia (VBAS, VINT, VSUP, DCAP, medias) sao feitas num laco
sobre IU, na mesma ordem do Fortran, para preservar a ordem das somas em
precisao simples.

Detalhes de fidelidade reproduzidos
-----------------------------------
* CondInic: TK e INTEGER no Fortran (truncamento de CB*3600 e CI*TIND).
* EVAPO: a variavel de modulo F4 nao e reinicializada; quando W <= WPM o
  Fortran usa RSX = 999999999*F4, com o F4 do ultimo par (IC,IU) que o
  atribuiu (inclusive de passos de tempo anteriores). Reproduzido com
  "forward fill" na ordem (IC, IU).
* SOLO: o laco DO ... ENDDO de ajuste do balanco foi vetorizado com mascaras,
  repetindo a mesma sequencia de ramos.
"""
from __future__ import annotations

import numpy as np

from .fortran_io import caldat, julday
from . import fmath as fm


# ---------------------------------------------------------------------------
# PARCEL + INTECLIM
# ---------------------------------------------------------------------------
def parcel(M):
    R = M.R
    XLADO = M.LCEL                         # km
    DH = M.HCEL * M.LCEL                   # m
    M.TIND[:] = R(3600.) * fm.pw((R(0.868) * fm.pw(XLADO, R(3.0))) / DH, R(0.385))
    inteclim(M)
    # Blocos "if(1==0)" do Fortran (REGION e leitura de GEOR.RG) estao
    # desativados no original e nao sao executados aqui (ver region()).


def region(ARX):
    """REGION.f90 (nao usada na simulacao: bloco desativado em PARCEL)."""
    QMESP = 30.
    QRX = QMESP * ARX / 1000.
    BRX = 0.8054 * ARX ** 0.5289
    HRX = 0.9 * ARX ** 0.19
    return QRX, BRX, HRX


def inteclim(M):
    R = M.R
    for ic in range(M.NC):
        dmin = R(1000000.0)
        ipmin = 0
        for ip in range(M.NCLI):
            dist = fm.pw(fm.pw(M.X[ic] - M.XYC[ip, 0], R(2.)) + fm.pw(M.Y[ic] - M.XYC[ip, 1], R(2.)), R(0.5))
            if dist < dmin:
                dmin = dist
                ipmin = ip
        M.ICBOM[ic] = ipmin


# ---------------------------------------------------------------------------
# PARCUNGE - parametros de Muskingum-Cunge
# ---------------------------------------------------------------------------
def parcunge(M):
    R = M.R
    M.ICODMUSK = np.zeros(M.NC, int)
    for ic in range(M.NC):
        if M.OD[ic] > 1:
            COEF = R(1.67) * fm.pw(M.DECL[ic], R(0.3)) / fm.pw(M.RUGMAN[ic], R(0.6))
            M.CEL[ic] = COEF * fm.pw(M.QREF[ic] / M.BRIO[ic], R(0.4))
            AT = R(M.DTP)
            DTCAL = AT / R(24.)
            QRX = M.QREF[ic]
            DECLIV = M.DECL[ic]
            BRX = M.BRIO[ic]
            CLX = M.CEL[ic]
            DX = R(0.5) * CLX * DTCAL * (R(1) + fm.pw(R(1) + R(1.5) * (QRX / BRX) /
                                                     (fm.pw(CLX, R(2)) * DECLIV * DTCAL), R(0.5)))
            ALONG = M.SRIO[ic] * R(1000.)
            x = float(ALONG / DX)
            NTRECH = int(np.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)   # NINT
            NTRECH = max(NTRECH, 1)
            if NTRECH > M.NUMUSK:
                raise RuntimeError(f"ntrech number greater then {M.NUMUSK} in parcunge (IC={ic + 1})")
            M.DT[ic] = DTCAL
            M.NSUBT[ic] = NTRECH
        else:
            M.DT[ic] = 0.0
            M.NSUBT[ic] = 0


# ---------------------------------------------------------------------------
# CONDINIC - condicoes iniciais
# ---------------------------------------------------------------------------
def condinic(M):
    R = M.R
    M.PM[:] = 0; M.PM2[:] = 0; M.SIM[:] = 0; M.EM[:] = 0; M.WBM[:] = 0
    M.DSUM[:] = 0; M.DINM[:] = 0; M.DBAM[:] = 0
    M.SI[:] = 0
    ib = M.IBAC
    M.W[:, :] = M.WM[ib, :] * R(0.40)
    M.VBAS[:] = M.QESP[ib] * M.ACEL * M.CB[ib] * R(3600.)
    M.VINT[:] = 0.0
    M.VSUP[:] = 0.0
    # TK e INTEGER no Fortran -> truncamento
    TK = np.trunc(M.CB[ib] * R(3600.)).astype(np.int64)
    M.QBAS[:] = M.VBAS / TK
    TK = np.trunc(M.CI[ib] * M.TIND).astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        M.QINT[:] = M.VINT / TK
    TK = np.trunc(M.CS[ib] * M.TIND).astype(np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        M.QSUP[:] = M.VSUP / TK
    M.TA[:] = 20.0
    M.QM2[:] = 0; M.PMB2[:] = 0; M.PMI2[:] = 0; M.PMS2[:] = 0
    M.QM1[:] = 0; M.QJ1[:] = 0; M.QJ2[:] = 0; M.QCEL1[:] = 0; M.QCEL2[:] = 0
    M.QRIOINI[:] = 0; M.QCONTORM[:] = 0; M.QCONTORJ[:] = 0
    subst = {int(c): i for i, c in enumerate(M.ISUBST)}
    for ic in range(M.NC):
        M.QCEL2[ic] = M.QBAS[ic] + M.QINT[ic] + M.QSUP[ic]
        M.QCEL1[ic] = M.QCEL2[ic]
        jus = M.jus[ic]
        if M.NSUBT[ic] > 0:
            M.QM2[ic] = M.QM2[ic] + M.QCEL2[ic]
            M.QJ2[ic] = M.QM2[ic]
            if jus >= 0:
                M.QM2[jus] = M.QM2[jus] + M.QJ2[ic]
        else:
            M.QJ2[ic] = M.QCEL2[ic]
            M.QM2[ic] = M.QCEL2[ic]
            if jus >= 0:
                M.QM2[jus] = M.QM2[jus] + M.QCEL2[ic]
        if (ic + 1) in subst:
            # (sic) no Fortran usa JUS mesmo se JUS<=0
            M.QM2[jus] = M.QLIDO[subst[ic + 1], 0]
    M.QJ1[:] = M.QJ2
    M.QM1[:] = M.QM2
    M.QRIOINI[:M.NC, :] = M.QM2[:M.NC, None]
    M.QCONTORM[:M.NC, :] = M.QM2[:M.NC, None]
    M.QCONTORJ[:M.NC, :] = M.QJ2[:M.NC, None]
    if M.hdFLAG0 > 0:
        from .inertial import flood_ini
        flood_ini(M)


# ---------------------------------------------------------------------------
# LECLIMA - clima do intervalo
# ---------------------------------------------------------------------------
def leclima(M):
    mm, _, _ = caldat(M.JDIA)
    k = M.ICBOM
    m = mm - 1
    it = M.IT - 1
    if M.flagaclimed == 1:
        M.TA[:] = M.TAMM[k, m]
        M.UR[:] = M.URMM[k, m]
        M.VV[:] = M.VVMM[k, m]
        M.SOL[:] = M.SOLMM[k, m]
        M.PA[:] = M.PAMM[k, m]
    else:
        td = M.TD[k, it]; ud = M.UD[k, it]; vd = M.VD[k, it]; sd = M.SOLD[k, it]; pd = M.PAD[k, it]
        M.TA[:] = np.where(td < -200.0, M.TAMM[k, m], td)
        M.UR[:] = np.where(ud < 0.0, M.URMM[k, m], ud)
        M.VV[:] = np.where(vd < 0.0, M.VVMM[k, m], vd)
        M.SOL[:] = np.where(sd < 0.0, M.SOLMM[k, m], sd)
        M.PA[:] = np.where(pd < 0.0, M.PAMM[k, m], pd)


# ---------------------------------------------------------------------------
# RADIACAO + EVAPO + SOLO (vetorizados em (NC,NU))
# ---------------------------------------------------------------------------
def _sind(x):
    return np.sin(np.deg2rad(x))


def _cosd(x):
    return np.cos(np.deg2rad(x))


def _tand(x):
    return np.tan(np.deg2rad(x))


def radiacao(M, T2, T1, URX, GLAT, SOLX, ALBX):
    """Retorna (RLX, CLATE, ES, ED, STO, SSUP). Entradas com broadcast."""
    R = M.R
    PI = R(M.PI)
    CLATE = (R(2.501) - R(0.002361) * T2)
    ES = R(0.6108) * fm.exp((R(17.27) * T2) / (R(237.3) + T2))
    ED = ES * URX / R(100.0)
    RDIA = R(M.JDIA)
    SDECL = R(0.4093) * fm.sin((R(2) * PI / R(365.0)) * RDIA - R(1.405))
    SANG = fm.arccos(-fm.tand(GLAT) * fm.tan(SDECL))
    HIM = R(24) * SANG / PI
    DR = R(1.0) + R(0.033) * fm.cos(R(2) * PI * RDIA / R(365.0))
    STO = R(15.392) * DR * (SANG * fm.sind(GLAT) * fm.sin(SDECL) +
                            fm.cosd(GLAT) * fm.cos(SDECL) * fm.sin(SANG))
    STO = STO / R(1000.0)
    STO = STO * R(M.MESP) * CLATE
    if M.flagaclimed == 1:
        if M.CRU_INDEX == 1:
            SSUP = STO * (R(0.25) + R(0.5) * SOLX / R(100.))
        else:
            SSUP = STO * (R(0.25) + R(0.5) * SOLX / HIM)
    else:
        SSUP = STO * (R(0.25) + R(0.5) * SOLX / HIM)
    SN = SSUP * (R(1) - ALBX)
    SLONG = R(M.STEBOL) * fm.pw(T2 + R(273.2), R(4.0))
    SLONG = SLONG * (R(0.34) - R(0.14) * fm.pw(ED, R(0.5)))
    SLONG = SLONG * (R(1.0) * SSUP / STO + R(0.0))
    GS = R(0.38) * (T2 - T1)
    RLX = SN - SLONG - GS
    return RLX, CLATE, ES, ED, STO, SSUP


def _stale_forward_fill(assigned_mask, values, carry):
    """Para cada posicao (ordem achatada), devolve o ultimo valor atribuido
    estritamente antes dela (ou ``carry`` se nenhum). Retorna tambem o novo
    carry (ultimo valor atribuido no vetor)."""
    flat_m = assigned_mask.ravel()
    flat_v = values.ravel()
    n = flat_m.size
    idx = np.where(flat_m, np.arange(n), -1)
    last = np.maximum.accumulate(idx)
    # ultimo atribuido estritamente antes de k
    prev = np.empty(n, dtype=np.int64)
    prev[0] = -1
    prev[1:] = last[:-1]
    out = np.where(prev >= 0, flat_v[np.maximum(prev, 0)], carry)
    new_carry = flat_v[last[-1]] if last[-1] >= 0 else carry
    return out.reshape(assigned_mask.shape), new_carry


def evapo(M, valid, ES, ED, PAX, CLATE, TAX, RLX, ZX, VVX, WMX, WX, RSX, RIAFX, SIX, PX):
    """EVAPO vetorizado. Retorna dict com PX, SIX, EIX, EIXP, REIX, EX, ETX."""
    R = M.R
    D = ES - ED
    GAMA = R(0.0016286) * PAX / CLATE
    DEL = (R(4098.0) * ES) / fm.pw(R(237.3) + TAX, R(2.))
    MSPA = R(3.486) * (PAX / (R(275.0) + TAX))
    RLX = RLX / (R(3600.0) * R(24.0))
    RUG = R(0.1) * ZX
    VVX = VVX + R(0.1)
    VVX2 = VVX / ((fm.log(R(10.0) / RUG)) / fm.log(R(2.0) / RUG))
    ZX = np.minimum(ZX, R(12.5))
    RA = ((fm.log((R(10.) - ZX * R(0.67)) / (ZX * R(0.123)))) *
          (fm.log((R(10.) - ZX * R(0.67)) / (ZX * R(0.0123))))) / (R(0.41) * R(0.41) * VVX)
    RAEP = (R(4.72) * fm.pw(fm.log(R(2.) / R(0.00137)), R(2.))) / (R(1.) + R(0.536) * VVX2)
    WL = R(0.5) * WMX
    WPM = R(0.1) * WMX
    # --- resistencia superficial com F4 "stale" ---------------------------
    soil = valid & ~(WMX < R(0.1))
    assigned = soil & (WX > WPM)
    with np.errstate(divide="ignore", invalid="ignore"):
        f4frac = (WX - WPM) / (WL - WPM)
        f4new = np.where(WX < WL, R(1.) / f4frac, R(1.0)).astype(R)
    f4stale, M.F4 = _stale_forward_fill(assigned, f4new, M.F4)
    f4stale = f4stale.astype(R)
    RSXo = np.where(WMX < R(0.1), R(0.0),
                    np.where(WX <= WPM, f4stale * R(999999999.0), f4new * RSX)).astype(R)
    SIL = R(0.2) * RIAFX
    SF = np.minimum(SIX + PX, SIL)
    PX = PX - (SF - SIX)
    fac = ((R(1000.0) * R(M.DTP)) / (CLATE * R(M.MESP)))
    EIX = ((DEL * RLX + MSPA * R(M.CSPA) * D / RA) / (DEL + GAMA * (R(1.0) + R(0.0) / RA))) * fac
    EIXP = ((DEL * RLX + MSPA * R(M.CSPA) * D / RAEP) / (DEL + GAMA * (R(1.0) + R(0.0) / RAEP))) * fac
    REIX = np.minimum(EIX, SF)
    SIX = SF - REIX
    with np.errstate(divide="ignore", invalid="ignore"):
        FDE = (EIX - REIX) / EIX
    EX = ((DEL * RLX + MSPA * R(M.CSPA) * D / RA) / (DEL + GAMA * (R(1.0) + RSXo / RA))) * fac
    EX = FDE * EX
    return PX, SIX, EIX, EIXP, REIX, EX


def solo(M, PX, EX, WX, WMX, BX, KIX, KBX, XL, CAPX, WCX, DTP):
    """SOLO vetorizado (modelo ARNO/VIC simplificado). Entradas 1D (pares ativos)."""
    R = M.R
    one = R(1.)
    WB = WCX * WMX
    WZ = WB
    WC = WCX * WMX
    with np.errstate(invalid="ignore", divide="ignore"):
        VTEMP = fm.pw(one - WX / WMX, one / (BX + one)) - (PX / ((BX + one) * WMX))
        d_sat = np.maximum(PX - (WMX - WX), R(0.0))
        d_uns = np.maximum(PX - (WMX - WX) + WMX * fm.pw(VTEMP, BX + one), R(0.0))
    DSUP = np.where(PX > 0.0, np.where(VTEMP <= 0.0, d_sat, d_uns), R(0.0)).astype(R)
    with np.errstate(invalid="ignore", divide="ignore"):
        DBAS = np.where(WX >= WB, KBX * (WX - WB) / (WMX - WB), R(0.0)).astype(R)
        CAPX = np.maximum(CAPX * (WC - WX) / WC, R(0.0)).astype(R)
        dint = KIX * fm.pw((WX - WZ) / (WMX - WZ), R(3.) + R(2.) / XL)
    DINT = np.where(WX <= WZ, R(0.0), dint).astype(R)
    DINT = np.minimum(DINT, WX / R(10.0))
    DINT = np.maximum(DINT, R(0.0))
    WXOLD = WX.copy()
    # Fortran: CAPX=CAPX*DTP/86400. (avaliado da esquerda para a direita)
    CAPX = (CAPX * R(DTP)) / R(86400.)
    DINT = (DINT * R(DTP)) / R(86400.)
    DBAS = (DBAS * R(DTP)) / R(86400.)
    EX = EX.copy()
    WXout = WX.copy()
    active = np.ones(WX.shape, bool)
    for _ in range(1000):
        if not active.any():
            break
        a = active
        wx = (WXOLD[a] + PX[a]) + ((((CAPX[a] - EX[a]) - DSUP[a]) - DINT[a]) - DBAS[a])
        wmx = WMX[a]
        c = (wx < wmx + R(0.01)) & (wx > wmx)
        wx = np.where(c, wmx, wx)
        c = (wx > R(-0.01)) & (wx < R(0.0))
        wx = np.where(c, R(0.0), wx)
        over = wx > wmx
        under = (~over) & (wx < 0.0)
        normal = ~(over | under)
        ia = np.nonzero(a)[0]
        # excesso: escoa superficialmente e repete
        io = ia[over]
        DSUP[io] = DSUP[io] + wx[over] - wmx[over]
        # deficit: retira de EX, DBAS, DINT, DSUP (nesta ordem)
        iu = ia[under]
        wxu = wx[under]
        keep = np.zeros(iu.size, bool)
        dry = np.zeros(iu.size, bool)
        ex = EX[iu] - (-wxu)
        EX[iu] = ex
        k1 = ex >= 0.0
        keep |= k1
        r = ~k1                      # residual = ex (<0)
        if r.any():
            ir = iu[r]
            wres = ex[r]
            EX[ir] = 0.0
            db = DBAS[ir] - (-wres)
            DBAS[ir] = db
            k2 = db >= 0.0
            kk = np.zeros(iu.size, bool); kk[np.nonzero(r)[0][k2]] = True
            keep |= kk
            r2 = ~k2
            if r2.any():
                ir2 = ir[r2]
                wres2 = db[r2]
                DBAS[ir2] = 0.0
                di = DINT[ir2] - (-wres2)
                DINT[ir2] = di
                k3 = di >= 0.0
                pos2 = np.nonzero(r)[0][r2]
                kk = np.zeros(iu.size, bool); kk[pos2[k3]] = True
                keep |= kk
                r3 = ~k3
                if r3.any():
                    ir3 = ir2[r3]
                    wres3 = di[r3]
                    DINT[ir3] = 0.0
                    ds = DSUP[ir3] - (-wres3)
                    DSUP[ir3] = ds
                    k4 = ds >= 0.0
                    pos3 = pos2[r3]
                    kk = np.zeros(iu.size, bool); kk[pos3[k4]] = True
                    keep |= kk
                    r4 = ~k4
                    if r4.any():
                        ir4 = ir3[r4]
                        DSUP[ir4] = 0.0
                        WXout[ir4] = R(0.001)          # 'Soil went dry'
                        kk = np.zeros(iu.size, bool); kk[pos3[r4]] = True
                        dry |= kk
        # atualiza mascaras
        inormal = ia[normal]
        WXout[inormal] = wx[normal]
        new_active = np.zeros_like(active)
        new_active[io] = True
        new_active[iu[keep]] = True
        active = new_active
    return WXout, EX, DSUP, DINT, DBAS, CAPX


# ---------------------------------------------------------------------------
# CELULA
# ---------------------------------------------------------------------------
def celula(M):
    R = M.R
    NC, NU = M.NC, M.NU
    DTP = R(M.DTP)
    it = M.IT - 1
    M.JDIA = M.JDIA - julday(1, 1, M.IANO) + 1          # dia do ano
    ib = M.IBAC
    imes = M.IMES - 1
    act_cell = (M.sbtFLAG != 1)
    valid = (M.PUSO >= R(0.0001)) & act_cell[:, None]
    M.EVQ[act_cell] = 0.0

    # cargas (NC,NU)
    PX = np.broadcast_to(M.P[:, None], (NC, NU)).astype(R)
    WX = M.W.copy()
    SIX = M.SI.copy()
    WMX = M.WM[ib, :]
    PAX = M.PA[:, None]; TAX = M.TA[:, None]; VVX = M.VV[:, None]
    URX = M.UR[:, None]; SOLX = M.SOL[:, None]
    XL = M.PLAM[ib, :]
    T1 = M.TONTEM[:, None]; T2 = M.TA[:, None]
    RIAFX = M.RIAF[:, imes][None, :]; ALBX = M.ALB[:, imes][None, :]
    ZX = np.broadcast_to(M.Z[:, imes][None, :], (NC, NU)).astype(R)
    RSX = np.broadcast_to(M.RS[:, imes][None, :], (NC, NU)).astype(R)
    BX = M.B[ib, :]; KIX = M.KINS[ib, :]; KBX = M.KBAS[ib, :]
    CAPX = M.CAP[ib, :].copy(); WCX = M.WC[ib, :]

    ta = M.TA[act_cell]
    if (ta < -50.0).any() or (ta > 50.0).any():
        raise RuntimeError("ERROR AIR TEMP SUB CELULA")
    if (M.VV[act_cell] < 0).any() or (M.VV[act_cell] > 100.).any():
        raise RuntimeError("ERROR WIND SUB CELULA")
    if (M.UR[act_cell] < 0).any() or (M.UR[act_cell] > 100.01).any():
        raise RuntimeError("ERROR UMIDITY SUB CELULA")

    GLAT = M.Y[:, None]
    RLX, CLATE, ES, ED, STO, SSUP = radiacao(M, T2, T1, URX, GLAT, SOLX, ALBX)
    RLX = np.broadcast_to(RLX, (NC, NU)).astype(R)
    PX, SIX_n, EIX, EIXP, REIX, EX = evapo(M, valid, ES, ED, PAX, CLATE, TAX, RLX, ZX, VVX,
                                           WMX, WX, RSX, RIAFX, SIX, PX)

    soilm = valid & (WMX > R(0.001))
    water = valid & ~(WMX > R(0.001))
    DSUP = M.zeros(NC, NU); DINT = M.zeros(NC, NU); DBAS = M.zeros(NC, NU)
    CAPo = M.zeros(NC, NU)
    ETX = M.zeros(NC, NU)
    WXn = WX.copy()
    if soilm.any():
        ii = soilm
        wx_s, ex_s, dsup_s, dint_s, dbas_s, cap_s = solo(
            M, PX[ii], EX[ii], WX[ii], WMX[ii], BX[ii], KIX[ii], KBX[ii], XL[ii],
            CAPX[ii], WCX[ii], M.DTP)
        WXn[ii] = wx_s
        DSUP[ii] = dsup_s; DINT[ii] = dint_s; DBAS[ii] = dbas_s; CAPo[ii] = cap_s
        ETX[ii] = ex_s + REIX[ii]
    if water.any():
        DSUP[water] = PX[water]
        ETX[water] = EIX[water]
    # EVQ: ultimo HRU "agua" da minibacia sobrescreve (ordem IU)
    for iu in range(NU):
        w = water[:, iu]
        if w.any():
            if True:  # HDFLAG == 0 (Muskingum-Cunge)
                hd0 = w & (M.HDFLAG == 0)
                M.EVQ[hd0] = (EIXP[hd0, iu] * R(1000.) * M.ACEL[hd0] * (M.PUSO[hd0, iu] / R(100.))) / DTP
                hd1 = w & (M.HDFLAG != 0)
                M.EVQ[hd1] = EIXP[hd1, iu]

    M.ET[valid] = ETX[valid]
    M.W[valid] = WXn[valid]
    M.SI[valid] = SIX_n[valid]
    M.CAF[valid] = CAPo[valid]
    M.QSAUX[:, 2:][valid] = DSUP[valid]           # DSUP (mm) para o modulo de sedimentos

    # acumulacoes por minibacia, na ordem de IU (como no Fortran)
    WXB = M.zeros(NC); SIB = M.zeros(NC); ETB = M.zeros(NC); DCAP = M.zeros(NC)
    hd0c = (M.HDFLAG == 0)
    area_eff = np.where(hd0c, M.ACEL, np.maximum((M.ACEL - M.Area2.astype(R)), R(0.0))).astype(R)
    VBAS = M.VBAS; VINT = M.VINT; VSUP = M.VSUP
    for iu in range(NU):
        v = valid[:, iu]
        if not v.any():
            continue
        PPU = M.PUSO[:, iu] / R(100.)
        WXB = np.where(v, WXB + WXn[:, iu] * PPU, WXB)
        SIB = np.where(v, SIB + (M.P - PX[:, iu]) * PPU, SIB)
        ETB = np.where(v, ETB + ETX[:, iu] * PPU, ETB)
        dsup = ((DSUP[:, iu] * R(1000.)) * area_eff * (M.PUSO[:, iu] / R(100.))) / DTP
        dint = ((DINT[:, iu] * R(1000.)) * area_eff * (M.PUSO[:, iu] / R(100.))) / DTP
        dbas = ((DBAS[:, iu] * R(1000.)) * area_eff * (M.PUSO[:, iu] / R(100.))) / DTP
        DCAP = np.where(v, DCAP + (M.PUSO[:, iu] / R(100.)) * (CAPo[:, iu] * R(1000.) * M.ACEL) / DTP, DCAP)
        VBAS[v] = (VBAS + dbas * DTP)[v]
        VINT[v] = (VINT + dint * DTP)[v]
        VSUP[v] = (VSUP + dsup * DTP)[v]

    c = act_cell
    VBAS[c] = np.maximum(VBAS[c] - DCAP[c] * DTP, R(0.0))
    TKB = M.CB[ib] * R(3600.)
    TKI = M.CI[ib] * M.TIND
    TKS = M.CS[ib] * M.TIND
    QBAS = VBAS / TKB
    QINT = VINT / TKI
    QSUP = VSUP / TKS
    M.QSAUX[c, 0] = TKS[c]
    M.QSAUX[c, 1] = VSUP[c]
    VBX = VBAS - QBAS * DTP
    VIX = VINT - QINT * DTP
    VSX = VSUP - QSUP * DTP
    # superficial
    neg_s = np.where(hd0c, VSX < 0.0, (VSUP < 0.0) | (VSX < 0.0))
    QSUP = np.where(neg_s, VSUP / DTP, QSUP)
    VSUPn = np.where(neg_s, R(0.0), VSX)
    neg_i = VIX < 0.0
    QINT = np.where(neg_i, VINT / DTP, QINT)
    VINTn = np.where(neg_i, R(0.0), VIX)
    neg_b = VBX < 0.0
    QBAS = np.where(neg_b, VBAS / DTP, QBAS)
    VBASn = np.where(neg_b, R(0.0), VBX)
    M.QSUP[c] = QSUP[c]; M.QINT[c] = QINT[c]; M.QBAS[c] = QBAS[c]
    M.VSUP[c] = VSUPn[c]; M.VINT[c] = VINTn[c]; M.VBAS[c] = VBASn[c]

    if M.EVAPTUDO is not None:
        M.EVAPTUDO[c, it] = ETB[c]
    if M.WTUDO is not None:
        M.WTUDO[c, it] = M.VINT[c]

    # medias por sub-bacia (BACIAS.HIG). Acumulacao sequencial em IC.
    kcb = M.KCB[ib].astype(R)
    DSB = (M.QSUP * R(3600.) * R(24.)) / (R(1000.) * M.ACEL)
    DIB = (M.QINT * R(3600.) * R(24.)) / (R(1000.) * M.ACEL)
    DBB = (M.QBAS * R(3600.) * R(24.)) / (R(1000.) * M.ACEL)
    for arr, val in ((M.PM2, M.P / kcb), (M.SIM, SIB / kcb), (M.EM, ETB / kcb),
                     (M.WBM, WXB / kcb), (M.DSUM, DSB / kcb), (M.DINM, DIB / kcb),
                     (M.DBAM, DBB / kcb)):
        col = arr[:, it]
        np.add.at(col, ib[c], val[c])
        arr[:, it] = col

    # E0 (ultimo IU processado da minibacia - EIXP/STO/SSUP sao globais)
    last_iu = np.where(valid.any(axis=1), NU - 1 - np.argmax(valid[:, ::-1], axis=1), -1)
    rows = np.arange(NC)
    has = last_iu >= 0
    li = np.maximum(last_iu, 0)
    STOb = np.broadcast_to(STO, (NC, NU)); SSUPb = np.broadcast_to(SSUP, (NC, NU))
    CLb = np.broadcast_to(CLATE, (NC, NU))
    M.E0agua[has] = EIXP[rows, li][has]
    M.E0TOPO[has] = (R(1000.) * STOb[rows, li] / (R(M.MESP) * CLb[rows, li]))[has]
    M.E0SUP[has] = (R(1000.) * SSUPb[rows, li] / (R(M.MESP) * CLb[rows, li]))[has]
