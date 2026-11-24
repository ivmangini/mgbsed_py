# -*- coding: utf-8 -*-
"""
Roda o MGB-SED (versao Python) no Spyder ou em qualquer IDE.

1. Instale o pacote uma vez:  pip install .   (na pasta do repositorio)
2. Ajuste as pastas abaixo e rode (F5 no Spyder).
"""
import time

import numpy as np

from mgbsed.main import setup, run

PASTA_INPUT = r"C:\caminho\para\input"
PASTA_OUTPUT = r"C:\caminho\para\output"     # criada se nao existir

PRECISAO = np.float32            # np.float32 = igual ao Fortran | np.float64 = precisao dupla
GRAVAR_ARQUIVOS_GRANDES = True   # False: nao grava MUSLE_*, QSUP, QCEL, CONC_*, DEPdin, VAZAO_MINI (~5 GB)
NT_TESTE = None                  # ex.: 400 para um teste rapido; None = periodo completo

t0 = time.time()
M = setup(PASTA_INPUT, PASTA_OUTPUT, real=PRECISAO,
          write_big_outputs=GRAVAR_ARQUIVOS_GRANDES, nt_override=NT_TESTE)
run(M)
print(f"\nConcluido em {(time.time() - t0) / 60:.1f} min. Saidas em: {PASTA_OUTPUT}")

# Estado do modelo ao final (exemplos):
#   M.QRG       vazao nos postos de VAZAO.HIG (n_postos, NT), m3/s
#   M.SEDstore  lista; SEDstore[0] = QSS nos postos (n_postos, NT), t/dia
#   M.CSJ2      concentracao areia/silte/argila no ultimo dia (NC, 3), t/m3
