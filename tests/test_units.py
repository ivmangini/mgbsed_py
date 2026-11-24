"""
Testes unitarios das partes vetorizadas contra traducoes escalares literais
(linha a linha) das rotinas Fortran.

Executar:  python -m pytest -q tests/   (ou  python tests/test_units.py)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mgbsed.fortran_io import julday, caldat, fmt_F, fmt_ES, fmt_A, read_fixed  # noqa: E402
from mgbsed import hydrology as hy  # noqa: E402
from mgbsed import routing as rt  # noqa: E402
from mgbsed import fmath as fm  # noqa: E402


# ---------------------------------------------------------------------------
def test_calendar_roundtrip():
    for jd in range(julday(1, 1, 1950), julday(12, 31, 2030)):
        mm, dd, yy = caldat(jd)
        assert julday(mm, dd, yy) == jd
    assert caldat(julday(1, 1, 1974)) == (1, 1, 1974)
    assert julday(2, 29, 2000) - julday(2, 28, 2000) == 1


def test_formats():
    assert fmt_F(10.786, 15, 3) == "         10.786"
    assert fmt_F(-0.0001, 15, 3) == "         -0.000"
    assert fmt_F(1e13, 15, 3) == "*" * 15
    assert fmt_F(float("nan"), 15, 3).strip() == "NaN"
    assert fmt_ES(12908.25, 16, 7) == "   1.2908250E+04"
    assert fmt_ES(0.0, 16, 7) == "   0.0000000E+00"
    assert fmt_A("IT", 10) == "        IT"
    # Fw.d sem ponto decimal -> casas implicitas
    assert read_fixed("       150", [("F", 10, 3)]) == [0.15]
    assert read_fixed("     150.0", [("F", 10, 3)]) == [150.0]
    assert read_fixed("   1   2\tx", [("I", 4, 3)], nitems=2) == [1, 2]


# ---------------------------------------------------------------------------
class _MS:
    R = np.float32


def solo_scalar(PX, EX, WX, WMX, BX, KIX, KBX, XL, CAPX, WCX, DTP):
    """Traducao literal de SOLO.f90 (escalar, float32)."""
    f = np.float32
    WB = WCX * WMX; WZ = WB; WC = WCX * WMX
    if PX > 0.0:
        VTEMP = fm.pw(f(1.) - WX / WMX, f(1.) / (BX + f(1.))) - (PX / ((BX + f(1.)) * WMX))
        if VTEMP <= 0.0:
            DSUP = max(PX - (WMX - WX), f(0.0))
        else:
            DSUP = max(PX - (WMX - WX) + WMX * fm.pw(VTEMP, BX + f(1.)), f(0.0))
    else:
        DSUP = f(0.0)
    DBAS = KBX * (WX - WB) / (WMX - WB) if WX >= WB else f(0.0)
    CAPX = max(CAPX * (WC - WX) / WC, f(0.0))
    DINT = f(0.0) if WX <= WZ else KIX * fm.pw((WX - WZ) / (WMX - WZ), f(3.) + f(2.) / XL)
    DINT = min(DINT, WX / f(10.0)); DINT = max(DINT, f(0.0))
    WXOLD = WX
    CAPX = CAPX * f(DTP) / f(86400.); DINT = DINT * f(DTP) / f(86400.); DBAS = DBAS * f(DTP) / f(86400.)
    while True:
        WX = WX + PX + ((((CAPX - EX) - DSUP) - DINT) - DBAS)
        if WMX < WX < WMX + f(0.01):
            WX = WMX
        if f(-0.01) < WX < f(0.0):
            WX = f(0.0)
        if WX > WMX:
            DSUP = DSUP + WX - WMX
            WX = WXOLD
        elif WX < 0.0:
            EX = EX - (-WX)
            if EX >= 0.0:
                WX = WXOLD; continue
            WX = EX; EX = f(0.0)
            DBAS = DBAS - (-WX)
            if DBAS >= 0.0:
                WX = WXOLD; continue
            WX = DBAS; DBAS = f(0.0)
            DINT = DINT - (-WX)
            if DINT >= 0.0:
                WX = WXOLD; continue
            WX = DINT; DINT = f(0.0)
            DSUP = DSUP - (-WX)
            if DSUP >= 0.0:
                WX = WXOLD; continue
            DSUP = f(0.0); WX = f(0.001)
            break
        else:
            break
    return WX, EX, DSUP, DINT, DBAS, CAPX


def test_solo_vectorized_matches_scalar():
    rng = np.random.default_rng(1)
    n = 4000
    f = np.float32
    WMX = rng.choice([50., 150., 550.], n).astype(f)
    WX = (rng.uniform(0, 1.0, n) * WMX).astype(f)
    WX[:200] = WMX[:200]                     # saturado
    WX[200:400] = 0.0005                     # quase seco
    PX = np.where(rng.uniform(size=n) < 0.5, 0, rng.exponential(20, n)).astype(f)
    EX = rng.uniform(0, 8, n).astype(f)
    BX = np.full(n, 0.08, f); KIX = rng.uniform(1, 20, n).astype(f); KBX = rng.uniform(0.1, 2, n).astype(f)
    XL = np.full(n, 0.38, f); CAPX = rng.uniform(0, 1, n).astype(f); WCX = np.full(n, 0.1, f)
    out = hy.solo(_MS, PX, EX, WX, WMX, BX, KIX, KBX, XL, CAPX, WCX, 86400.)
    for i in range(n):
        ref = solo_scalar(PX[i], EX[i], WX[i], WMX[i], BX[i], KIX[i], KBX[i], XL[i], CAPX[i], WCX[i], 86400.)
        for a, b in zip(out, ref):
            assert a[i] == b or (np.isnan(a[i]) and np.isnan(b)), (i, [o[i] for o in out], ref)


# ---------------------------------------------------------------------------
def test_stale_forward_fill():
    rng = np.random.default_rng(2)
    m = rng.uniform(size=(50, 15)) < 0.4
    v = rng.uniform(size=(50, 15)).astype(np.float32)
    out, carry = hy._stale_forward_fill(m, v, 7.0)
    last = 7.0
    for i in range(50):
        for j in range(15):
            assert out[i, j] == np.float32(last)
            if m[i, j]:
                last = v[i, j]
    assert carry == np.float32(last)


# ---------------------------------------------------------------------------
def test_musk_linear_equals_recurrence():
    """Matriz linear de REDE == recorrencia de MUSK quando nao ha valores negativos;
    e o fallback exato e usado quando ha."""
    class M:
        R = np.float64
        NC = 3
        NSUBT = np.array([1, 3, 4])
        DT = np.full(3, 3600.)
        DTP = 86400.
        jus = np.array([2, 2, -1])
        def zeros(self, *s, dtype=None):
            return np.zeros(s, dtype=dtype or np.float64)
    m = M()
    m.MC_C1 = np.array([0.6, 0.5, 0.45]); m.MC_C2 = np.array([0.2, -0.1, 0.25]); m.MC_C3 = np.array([0.2, 0.6, 0.3])
    m.EVQ = np.zeros(3)
    m.NTMUSK = 24
    # monta G como em prepare_routing
    T = 24; nmax = 4
    G = np.zeros((3, nmax * T, nmax + 1 + T))
    for ic in range(3):
        n = m.NSUBT[ic]
        QC = np.zeros((n + 1, T + 1, nmax + 1 + T))
        for r in range(n + 1):
            QC[r, 0, r] = 1
        for t in range(1, T + 1):
            QC[0, t, nmax + t] = 1
        for t in range(T):
            for r in range(n):
                QC[r + 1, t + 1] = m.MC_C1[ic] * QC[r, t] + m.MC_C2[ic] * QC[r, t + 1] + m.MC_C3[ic] * QC[r + 1, t]
        for r in range(1, n + 1):
            G[ic, (r - 1) * T:r * T] = QC[r, 1:]
    rng = np.random.default_rng(3)
    for ic in range(3):
        n = m.NSUBT[ic]
        qcm = rng.uniform(10, 20, T + 1)
        qr = rng.uniform(10, 20, n + 1)
        QC, cj = rt._musk_exact(m, ic, qcm, qr)
        x = np.zeros(nmax + 1 + T); x[:n + 1] = qr; x[nmax + 1:] = qcm[1:]
        grid = (G[ic] @ x).reshape(nmax, T)
        if (grid[:n] >= 0).all():
            assert np.allclose(grid[n - 1], cj[1:], rtol=1e-12)
            assert np.allclose(grid[:n, T - 1], QC[1:, T], rtol=1e-12)


if __name__ == "__main__":
    for k, fn in list(globals().items()):
        if k.startswith("test_"):
            fn()
            print("OK", k)
