"""
Estado global do modelo.

No Fortran, as variaveis vivem em modulos (VARS_MAIN, VARS_CALIB, VARS_INERC,
SED_VARS) com atributo SAVE. Aqui todas ficam como atributos de um unico objeto
``Model``; as sub-rotinas sao funcoes que recebem esse objeto (``M``).

Convencoes
----------
* Indices de minibacia (IC), sub-bacia (IB) e URH (IU) sao 0-based nos arrays.
  Os codigos lidos/gravados em arquivo continuam 1-based (como no Fortran).
* ``M.R`` e o tipo REAL padrao. Por padrao np.float32, igual ao REAL(4) do
  Fortran original (reprodutibilidade numerica). Pode ser np.float64.
  Variaveis declaradas REAL(8) no Fortran sao sempre np.float64.
"""
from __future__ import annotations

import numpy as np


class Model:
    # Parametros (VARS_MAIN)
    NUP = 15
    NBP = 505
    NCLIP = 501
    NUMUSK = 150
    NFP = 3
    NTMUP = 20
    NMAXHID = 800

    MESP = 1000.0           # massa especifica da agua (kg/m3)
    PI = 3.141592           # (sic) valor usado no Fortran
    STEBOL = 4.903E-9       # constante de Stefan-Boltzmann
    CSPA = 0.001013         # calor especifico do ar umido MJ/kg/C

    # VARS_INERC
    g = 9.81
    dxart = 5.0

    def __init__(self, input_dir: str, output_dir: str, real=np.float32,
                 write_big_outputs: bool = True, verbose: bool = True):
        self.INPUT_DIRECTORY = input_dir
        self.OUTPUT_DIRECTORY = output_dir
        self.R = real
        self.write_big_outputs = write_big_outputs
        self.verbose = verbose
        # escalares globais que o Fortran mantem entre chamadas
        self.F4 = 0.0       # EVAPO: fator F4 (estado "stale" entre URHs, ver hydrology.evapo)
        self.IT = 0
        self.xyz = 0

    # utilitarios ---------------------------------------------------------
    def zeros(self, *shape, dtype=None):
        return np.zeros(shape, dtype=self.R if dtype is None else dtype)

    def full(self, shape, val, dtype=None):
        return np.full(shape, val, dtype=self.R if dtype is None else dtype)

    def log(self, *msg):
        txt = " ".join(str(m) for m in msg)
        if self.verbose:
            print(txt, flush=True)
        if getattr(self, "_logf", None) is not None:
            self._logf.write(" " + txt + "\n")
            self._logf.flush()
