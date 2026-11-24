"""
Funcoes matematicas com arredondamento controlado.

O numpy calcula pow/exp/log em float32 com rotinas SIMD que podem errar 1 ulp
(cerca de 20% dos casos para pow). O Fortran (Intel) usa rotinas de precisao
simples com erro tipico < 0,5 ulp. Para aproximar o comportamento do Fortran,
quando os argumentos sao float32 as funcoes sao avaliadas em float64 e
arredondadas uma unica vez para float32 (resultado "corretamente arredondado"
na quase totalidade dos casos). Para float64, as funcoes do numpy sao usadas
diretamente.

Expoentes inteiros pequenos (2, 3, 4) sao calculados por multiplicacao, como o
compilador Intel faz para X**2.0, X**3.0 e X**4.0.
"""
from __future__ import annotations

import numpy as np

_F32 = np.dtype(np.float32)


def _dt(*xs):
    return np.result_type(*xs)


def pw(x, y):
    dt = _dt(x, y)
    if np.isscalar(y) or (isinstance(y, np.ndarray) and y.ndim == 0):
        yv = float(y)
        if yv == 2.0:
            return (x * x).astype(dt) if isinstance(x, np.ndarray) else dt.type(x * x)
        if yv == 3.0:
            r = x * x * x
            return r.astype(dt) if isinstance(r, np.ndarray) else dt.type(r)
        if yv == 4.0:
            t = x * x
            r = t * t
            return r.astype(dt) if isinstance(r, np.ndarray) else dt.type(r)
        if yv == 0.0:
            return np.ones_like(np.asarray(x), dtype=dt) if np.ndim(x) else dt.type(1.0)
    if dt == _F32:
        with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
            r = np.power(np.asarray(x, np.float64), np.asarray(y, np.float64))
        return r.astype(np.float32) if np.ndim(r) else np.float32(r)
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        return np.power(x, y)


def _u(fn, x):
    if _dt(x) == _F32:
        with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
            r = fn(np.asarray(x, np.float64))
        return r.astype(np.float32) if np.ndim(r) else np.float32(r)
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        return fn(x)


def exp(x):
    return _u(np.exp, x)


def log(x):
    return _u(np.log, x)


def log10(x):
    return _u(np.log10, x)


def sin(x):
    return _u(np.sin, x)


def cos(x):
    return _u(np.cos, x)


def tan(x):
    return _u(np.tan, x)


def arccos(x):
    return _u(np.arccos, x)


def sqrt(x):
    return _u(np.sqrt, x)


def sind(x):
    return _u(lambda v: np.sin(np.deg2rad(v)), x)


def cosd(x):
    return _u(lambda v: np.cos(np.deg2rad(v)), x)


def tand(x):
    return _u(lambda v: np.tan(np.deg2rad(v)), x)
