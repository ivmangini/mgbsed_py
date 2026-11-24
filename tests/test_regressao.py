"""
Teste de regressao: simulacao do rio Ivai por 60 dias (01/01 a 01/03/1974).

Dados (tests/dados_teste/):
   entrada/                 entradas do Ivai (topologia, solos, chuva de 400 dias e vazao
                            observada da ANA) com parametros GENERICOS, nao calibrados
                            (PARUSO.cal, PARUSO_MUSLE.txt, PARTEXT_MUSLE.txt; MUSLE com os
                            valores padrao, sem PAR_MUSLE_SED.txt) e observacoes de
                            sedimento SINTETICAS (QSS_T5.txt, CSS_T5.txt; uma a cada 5 dias).
                            infoMGB.sim com NT = 60.
   ref_fortran_mc/          saidas do executavel Fortran original (MGB-SED v1.0, Intel
   ref_fortran_inercial/    Fortran ifx 2026.1.1, Debug x64) com essas mesmas entradas
   ref_fortran_css/         e os coeficientes padrao da MUSLE no SED_BACIA.F90, rodado
   ref_fortran_400dias_mc/  por 400 dias: 60 primeiras linhas das saidas diarias, e as
   ref_fortran_400dias_css/ metricas e os acumulados de fim de rodada dos 400 dias.

Os 400 dias sao necessarios porque o Fortran so calcula as metricas de sedimento a
partir do dia 365 (SED_Fob_Calib, iniaju = 365).

O modelo depende so do passado: as 60 primeiras linhas de uma rodada longa sao
identicas as de uma rodada de 60 dias.

Criterios (por arquivo, sobre todos os valores):
    * mesma forma (linhas x colunas) e nenhum NaN a mais ou a menos;
    * correlacao de Pearson >= CORR_MIN;
    * fracao de valores com |Py - Ref| > 1 % * max(|Ref|, 1) <= FRAC_MAX;
    * erro absoluto maximo <= ABS_MAX[arquivo] (unidade do arquivo).

As tolerancias vem das diferencas medidas entre Python e Fortran (arredondamento
em float32 entre Intel Fortran e NumPy), com margem de ~10x para outras
plataformas. Ver VERIFICACAO.md, secao "Teste automatico".
"""
from __future__ import annotations

import gzip
import io
import os
import shutil

import numpy as np
import pytest

from mgbsed.main import run, setup

AQUI = os.path.dirname(os.path.abspath(__file__))
DADOS = os.path.join(AQUI, "dados_teste")

# arquivo: (linhas de cabecalho, colunas iniciais de data/indice a ignorar, erro abs. maximo)
ARQUIVOS = {
    "VAZAO.HIG": (1, 1, 0.005),          # m3/s
    "NOSOLO.txt": (0, 1, 0.0005),        # mm, m3/s
    "BACIAS.HIG": (0, 0, 0.0001),        # mm
    "QPROP.HIG": (0, 1, 0.005),          # m3/s (so MC; no inercial o Fortran grava memoria nao inicializada)
    "VAZAO_MINI.HIG": (1, 1, 0.01),      # m3/s
    "QSUP.txt": (1, 4, 0.005),           # m3/s
    "QCEL.txt": (1, 4, 0.005),           # m3/s
    "MUSLE_TOTAL.txt": (1, 4, 0.5),      # t/dia (medido: 0,034)
    "MUSLE_areia.txt": (1, 4, 0.5),
    "MUSLE_silte.txt": (1, 4, 0.5),
    "MUSLE_argila.txt": (1, 4, 0.5),
    "CONC_areia.txt": (1, 3, 50.0),      # mg/L (medido: 0,05 MC / 0,40 inercial; margem para trocas de ramo)
    "CONC_silte.txt": (1, 3, 5.0),
    "CONC_argila.txt": (1, 3, 5.0),
    "DEPdin_RIO.txt": (0, 3, 5.0),       # t (deposito acumulado de areia; medido: 0,07)
    "MINI_areia.txt": (1, 1, 0.05),      # mg/L nos postos
    "MINI_silte.txt": (1, 1, 0.05),
    "MINI_argila.txt": (1, 1, 0.05),
    "MINI_QSS.txt": (1, 1, 0.5),         # t/dia nos postos
}
CORR_MIN = 0.99999
FRAC_MAX = 1e-3


def _ler(linhas, cab, pular):
    linhas = [ln for ln in linhas[cab:] if ln.strip()]
    vals = []
    for ln in linhas:
        lin = []
        for t in ln.split()[pular:]:
            try:
                lin.append(float(t))
            except ValueError:
                lin.append(np.nan)
        vals.append(lin)
    return np.array(vals, dtype=float)


def _ler_ref(modo, nome, cab, pular):
    with gzip.open(os.path.join(DADOS, f"ref_fortran_{modo}", nome + ".gz"), "rt", encoding="latin-1") as f:
        return _ler(f.read().splitlines(), cab, pular)


def _ler_py(pasta, nome, cab, pular, nlin):
    with open(os.path.join(pasta, nome), encoding="latin-1") as f:
        linhas = f.read().splitlines()
    return _ler(linhas[:cab + nlin], cab, pular)


def _rodar(tmp, inercial, sedvar_css=False, nt=60):
    ent = os.path.join(tmp, "entrada")
    shutil.copytree(os.path.join(DADOS, "entrada"), ent)
    for nome in os.listdir(ent):                      # COTA_AREA.FLP.gz, CHUVABIN.PBI.gz
        if nome.endswith(".gz"):
            with gzip.open(os.path.join(ent, nome), "rb") as g, \
                    open(os.path.join(ent, nome[:-3]), "wb") as f:
                shutil.copyfileobj(g, f)
    if inercial:
        p = os.path.join(ent, "infoMGB.sim")
        s = open(p, "rb").read()
        a = b"      1816        15        12        -9         0         0         0         1"
        assert s.count(a) == 1
        open(p, "wb").write(s.replace(a, a[:-20] + b"         1         1"))
    if sedvar_css:
        # FLAG_SEDVAR = 1: compara CSS (mg/L) em vez de QSS e grava MINI_CSS.txt
        p = os.path.join(ent, "infoMGB.sim")
        s = open(p, "rb").read()
        a = b"         1         2\t"
        assert s.count(a) == 1
        s = s.replace(a, b"         1         1\t").replace(b"9 QSS_T5.txt", b"9 CSS_T5.txt")
        open(p, "wb").write(s)
    if nt != 60:
        p = os.path.join(ent, "infoMGB.sim")
        s = open(p, "rb").read()
        a = b"        60    86400."
        assert s.count(a) == 1
        open(p, "wb").write(s.replace(a, b"%10d    86400." % nt))
    saida = os.path.join(tmp, "saida")
    M = setup(ent, saida, real=np.float32, verbose=False)
    assert (M.IND_FLOODROUTING == 1) == inercial
    run(M, progress=False)
    return saida, M


@pytest.fixture(scope="module", params=["mc", "inercial"])
def rodada(request, tmp_path_factory):
    modo = request.param
    tmp = str(tmp_path_factory.mktemp(modo))
    saida, M = _rodar(tmp, inercial=(modo == "inercial"))
    return modo, saida, M


def comparar(F, P):
    """Estatisticas de comparacao entre referencia F e Python P."""
    assert F.shape == P.shape, f"forma diferente: Fortran {F.shape}, Python {P.shape}"
    nan_f, nan_p = np.isnan(F), np.isnan(P)
    m = ~(nan_f | nan_p)
    f, p = F[m], P[m]
    d = np.abs(p - f)
    if f.std() > 0 and p.std() > 0:
        corr = float(np.corrcoef(f, p)[0, 1])
    else:
        corr = 1.0 if np.array_equal(f, p) else 0.0
    return dict(
        nan_dif=int((nan_f != nan_p).sum()),
        max_abs=float(d.max()) if d.size else 0.0,
        frac_1pct=float((d > 0.01 * np.maximum(np.abs(f), 1.0)).mean()) if d.size else 0.0,
        frac_ident=float((d == 0).mean()) if d.size else 1.0,
        corr=corr,
    )


@pytest.mark.parametrize("nome", list(ARQUIVOS))
def test_saida_igual_fortran(rodada, nome):
    modo, saida, M = rodada
    if modo == "inercial" and nome == "QPROP.HIG":
        pytest.skip("QPROP.HIG nao tem significado no modo inercial (memoria nao inicializada no Fortran)")
    cab, pular, abs_max = ARQUIVOS[nome]
    F = _ler_ref(modo, nome, cab, pular)
    P = _ler_py(saida, nome, cab, pular, F.shape[0])
    st = comparar(F, P)
    msg = f"{modo}/{nome}: {st}"
    assert st["nan_dif"] == 0, msg
    assert st["corr"] >= CORR_MIN, msg
    assert st["frac_1pct"] <= FRAC_MAX, msg
    assert st["max_abs"] <= abs_max, msg


def test_balanco_massa_sedimentos(rodada):
    """Entrada - saida de sedimentos na rede fecha (MODELO.f90, verificacao final)."""
    modo, saida, M = rodada
    bal = M.mass_balance
    erro = 100 * (bal["IN"] - bal["OUT"]) / bal["IN"]
    # medido em 60 dias: |erro| <= 0,004 % (areia, inercial); serie completa: <= 0,003 %
    assert np.all(np.abs(erro) < 0.01), erro


def test_arquivos_gerados(rodada):
    modo, saida, M = rodada
    for nome in ("AJUSTE.fob", "AJUSTESED.fob", "AJUSTE_HIDRO_PY.fob", "SEDRIO_SAIDA.txt",
                 "SEDRIO_DEPOSITO.txt", "SEDRIO_EROSAO.txt", "BAL_PLAN.txt", "LOG_MGB.txt"):
        assert os.path.getsize(os.path.join(saida, nome)) > 0, nome


def test_modo_css(tmp_path):
    """FLAG_SEDVAR = 1 (CSS): MINI_CSS.txt igual ao do Fortran e igual a silte + argila.

    O FLAG_SEDVAR so muda a variavel gravada e comparada (MODELO.f90, SIMULA.f90):
    CSS = (CSJ2(silte) + CSJ2(argila)) * 1e6.
    """
    saida, M = _rodar(str(tmp_path), inercial=False, sedvar_css=True)
    assert M.FLAG_SEDVAR == 1
    css = _ler_py(saida, "MINI_CSS.txt", 1, 1, 60)
    F = _ler_ref("css", "MINI_CSS.txt", 1, 1)
    st = comparar(F, css)
    assert st["nan_dif"] == 0 and st["corr"] >= CORR_MIN and st["max_abs"] <= 0.05, st
    soma = _ler_ref("mc", "MINI_silte.txt", 1, 1) + _ler_ref("mc", "MINI_argila.txt", 1, 1)
    # silte e argila gravados com 3 casas: a soma pode diferir de ate ~0,002 + arredondamento
    assert np.abs(css - soma).max() <= 0.01
    assert not os.path.exists(os.path.join(saida, "MINI_QSS.txt"))


def _ler_txt(p):
    vals = []
    for ln in open(p, encoding="latin-1"):
        for t in ln.split():
            try:
                vals.append(float(t))
            except ValueError:
                vals.append(np.nan)
    return np.array(vals)


@pytest.mark.parametrize("modo", ["mc", "css"])
def test_metricas_400_dias(tmp_path, modo):
    """Rodada de 400 dias: AJUSTE.fob e AJUSTESED.fob identicos ao Fortran (texto gravado)
    e acumulados de fim de rodada (SEDRIO_*) dentro de 0,05 t."""
    saida, M = _rodar(str(tmp_path), inercial=False, sedvar_css=(modo == "css"), nt=400)
    ref = os.path.join(DADOS, f"ref_fortran_400dias_{modo}")
    for nome in ("AJUSTE.fob", "AJUSTESED.fob"):
        f = open(os.path.join(ref, nome), encoding="latin-1").read().replace("\r", "")
        p = open(os.path.join(saida, nome), encoding="latin-1").read().replace("\r", "")
        assert f.split() == p.split(), nome
    for nome in ("SEDRIO_SAIDA.txt", "SEDRIO_DEPOSITO.txt", "SEDRIO_EROSAO.txt"):
        F = _ler_txt(os.path.join(ref, nome)); P = _ler_txt(os.path.join(saida, nome))
        assert F.shape == P.shape, nome
        assert np.array_equal(np.isnan(F), np.isnan(P)), nome
        m = ~np.isnan(F)
        assert np.abs(F[m] - P[m]).max() <= 0.05, nome
