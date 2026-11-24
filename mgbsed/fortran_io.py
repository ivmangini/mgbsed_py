"""
Utilitarios de entrada/saida que reproduzem a semantica do Fortran (Intel Fortran)
usada no MGB-SED v1.0.

Conteudo
--------
* Leitura de registros com formato fixo (A, I, F, X) -> ``read_fixed``
* Leitura "list-directed" (``READ(u,*)``) -> ``read_list``
* Escrita formatada (Fw.d, Iw, ESw.d, Aw) -> ``fmt_F``, ``fmt_I``, ``fmt_ES``, ``fmt_A``
* Calendario juliano (Numerical Recipes): ``julday`` e ``caldat``

Observacoes de fidelidade
-------------------------
* Arquivos sao lidos como latin-1, pois o Fortran le bytes: nomes com acento
  (ex.: "Maringa" com a-agudo em 1 byte) precisam ocupar a mesma coluna.
* Campo Fw.d sem ponto decimal: o Fortran aplica d casas implicitas
  (ex.: "  150" em F8.3 -> 0.150). Isso e reproduzido aqui.
* Campos em branco valem 0 (BLANK='NULL', padrao).
* Registros mais curtos que o formato sao completados com brancos (PAD='YES').
"""
from __future__ import annotations

import math
import re
from typing import Iterable, List, Sequence

import numpy as np

ENCODING = "latin-1"


# ----------------------------------------------------------------------------
# Leitura
# ----------------------------------------------------------------------------
class FortranFile:
    """Leitor sequencial de registros (linhas) de um arquivo texto."""

    def __init__(self, path: str):
        with open(path, "r", encoding=ENCODING, newline="") as f:
            txt = f.read()
        # Fortran (Windows) trata CRLF como fim de registro
        self.lines: List[str] = txt.replace("\r\n", "\n").replace("\r", "\n").split("\n")
        if self.lines and self.lines[-1] == "":
            self.lines.pop()
        self.pos = 0
        self.path = path

    def next_record(self) -> str:
        if self.pos >= len(self.lines):
            raise EOFError(f"Fim de arquivo em {self.path} (registro {self.pos + 1})")
        rec = self.lines[self.pos]
        self.pos += 1
        return rec

    def skip(self, n: int = 1) -> None:
        for _ in range(n):
            self.next_record()

    def eof(self) -> bool:
        return self.pos >= len(self.lines)

    # READ(u,fmt)
    def read_fixed(self, spec: Sequence, nitems: int | None = None) -> list:
        return read_fixed(self.next_record(), spec, nitems)

    # READ(u,*)
    def read_list(self, n: int | None = None, types: Sequence | None = None) -> list:
        """Leitura list-directed: consome registros ate obter ``n`` itens.

        Se ``n`` for None, le apenas o primeiro registro (todos os itens).
        O restante do ultimo registro consumido e descartado (semantica Fortran).
        """
        items: List[str] = []
        if n is None:
            items = _split_list_directed(self.next_record())
        else:
            while len(items) < n:
                toks = _split_list_directed(self.next_record())
                items.extend(toks)
            items = items[:n]
        if types is None:
            return items
        out = []
        for i, t in enumerate(items):
            tp = types[i] if i < len(types) else types[-1]
            out.append(_conv(t, tp))
        return out


def _conv(tok: str, tp):
    if tp is str:
        return tok.strip("'\"")
    if tp is int:
        return int(float(tok)) if re.search(r"[.eEdD]", tok) else int(tok)
    return float(tok.replace("D", "E").replace("d", "e"))


def _split_list_directed(rec: str) -> List[str]:
    # separadores: brancos, virgula, tab. '/' encerra a leitura do registro.
    rec = rec.split("/")[0]
    toks = re.split(r"[\s,]+", rec.strip())
    out = []
    for t in toks:
        if t == "":
            continue
        m = re.match(r"^(\d+)\*(.*)$", t)  # repeticao r*c
        if m:
            out.extend([m.group(2)] * int(m.group(1)))
        else:
            out.append(t)
    return out


def _read_F(field: str, d: int) -> float:
    s = re.sub(r"\s", "", field)  # BLANK='NULL': brancos ignorados
    if s == "" or s in "+-":
        return 0.0
    s = s.replace("D", "E").replace("d", "e")
    if "." in s:
        return float(s)
    # sem ponto decimal: d casas implicitas na mantissa
    m = re.match(r"^([+-]?\d+)([eE][+-]?\d+|[+-]\d+)?$", s)
    if not m:
        return float(s)
    mant = int(m.group(1)) / (10.0 ** d)
    if m.group(2):
        e = m.group(2)
        e = e[1:] if e[0] in "eE" else e
        mant *= 10.0 ** int(e)
    return mant


def _read_I(field: str) -> int:
    s = re.sub(r"\s", "", field)
    if s == "" or s in "+-":
        return 0
    return int(s)


def expand_spec(spec: Sequence) -> list:
    """Expande ('F',8,3,n) em n descritores. Formatos: ('A',w), ('I',w),
    ('F',w,d), ('X',n). Um quarto elemento opcional e a repeticao."""
    out = []
    for item in spec:
        kind = item[0]
        if kind == "A":
            rep = item[2] if len(item) > 2 else 1
            out += [("A", item[1])] * rep
        elif kind == "I":
            rep = item[2] if len(item) > 2 else 1
            out += [("I", item[1])] * rep
        elif kind == "F":
            rep = item[3] if len(item) > 3 else 1
            out += [("F", item[1], item[2])] * rep
        elif kind == "X":
            out.append(("X", item[1]))
        else:
            raise ValueError(item)
    return out


def read_fixed(rec: str, spec: Sequence, nitems: int | None = None) -> list:
    """Le um registro segundo um formato fixo (lista de descritores).

    ``nitems``: numero de itens da lista de E/S. Como no Fortran, o formato e
    processado apenas ate esgotar a lista (campos restantes nao sao lidos).
    """
    desc = expand_spec(spec)
    if nitems is not None:
        out_desc, k = [], 0
        for d in desc:
            if k >= nitems:
                break
            out_desc.append(d)
            if d[0] != "X":
                k += 1
        desc = out_desc
    width = sum(d[1] for d in desc)
    rec = rec.ljust(width)  # PAD='YES'
    vals = []
    p = 0
    for d in desc:
        w = d[1]
        field = rec[p:p + w]
        p += w
        if d[0] == "A":
            vals.append(field)
        elif d[0] == "I":
            vals.append(_read_I(field))
        elif d[0] == "F":
            vals.append(_read_F(field, d[2]))
        # X: ignora
    return vals


# ----------------------------------------------------------------------------
# Escrita formatada
# ----------------------------------------------------------------------------
def fmt_A(s: str, w: int) -> str:
    """Saida Aw: texto menor que w e alinhado a direita; maior, truncado a esquerda."""
    s = str(s)
    if len(s) >= w:
        return s[:w]
    return s.rjust(w)


def fmt_I(v, w: int) -> str:
    s = str(int(v))
    return s.rjust(w) if len(s) <= w else "*" * w


def fmt_F(v, w: int, d: int) -> str:
    v = float(v)
    if math.isnan(v):
        return "NaN".rjust(w)
    if math.isinf(v):
        s = "Infinity" if w >= 9 else "Inf"
        s = ("-" if v < 0 else "") + s
        return s.rjust(w) if len(s) <= w else "*" * w
    s = f"{v:.{d}f}"
    if s.startswith("-") and float(s) == 0.0:
        s = s  # Intel mantem o sinal: -0.000
    if len(s) > w:
        # Fortran omite o zero inicial opcional (0.123 -> .123) antes de estourar
        if s.startswith("0.") and len(s) - 1 <= w:
            s = s[1:]
        elif s.startswith("-0.") and len(s) - 1 <= w:
            s = "-" + s[2:]
        else:
            return "*" * w
    return s.rjust(w)


def fmt_ES(v, w: int, d: int) -> str:
    v = float(v)
    if math.isnan(v):
        return "NaN".rjust(w)
    if math.isinf(v):
        return (("-" if v < 0 else "") + "Infinity").rjust(w)
    s = f"{v:.{d}E}"
    mant, exp = s.split("E")
    e = int(exp)
    if abs(e) < 100:
        s = f"{mant}E{'+' if e >= 0 else '-'}{abs(e):02d}"
    else:
        s = f"{mant}{'+' if e >= 0 else '-'}{abs(e):03d}"
    return s.rjust(w) if len(s) <= w else "*" * w


# Versoes vetorizadas (rapidas) para matrizes grandes -------------------------
def format_F_matrix(arr: np.ndarray, w: int, d: int) -> List[str]:
    """Formata cada linha de ``arr`` (2D) como concatenacao de campos Fw.d.

    Usa formatacao em lote do numpy e corrige NaN/Inf/estouro de campo.
    """
    arr = np.asarray(arr, dtype=np.float64)
    if arr.ndim == 1:
        arr = arr[None, :]
    fmt = f"%{w}.{d}f"
    out = []
    # limite de valor que cabe no campo
    big = 10.0 ** (w - d - 2) - 0.5 * 10.0 ** (-d)
    special = (~np.isfinite(arr)) | (np.abs(arr) >= big)
    for i in range(arr.shape[0]):
        row = arr[i]
        if special[i].any():
            out.append("".join(fmt_F(x, w, d) for x in row))
        else:
            out.append((fmt * row.size) % tuple(row.tolist()))
    return out


# ----------------------------------------------------------------------------
# Calendario (Numerical Recipes) - julday.f90 / caldat.f90
# ----------------------------------------------------------------------------
def _fint(x: float) -> int:
    """INT() do Fortran: truncamento em direcao a zero."""
    return int(x)


def julday(mm: int, iday: int, iyyy: int) -> int:
    IGREG = 15 + 31 * (10 + 12 * 1582)
    jy = iyyy
    if jy == 0:
        raise ValueError("julday: there is no year zero")
    if jy < 0:
        jy += 1
    if mm > 2:
        jm = mm + 1
    else:
        jy -= 1
        jm = mm + 13
    jul = _fint(365.25 * jy) + _fint(30.6001 * jm) + iday + 1720995
    if iday + 31 * (mm + 12 * iyyy) >= IGREG:
        ja = _fint(0.01 * jy)
        jul = jul + 2 - ja + _fint(0.25 * ja)
    return jul


def caldat(julian: int):
    """Retorna (mm, id, iyyy) - mesma ordem de argumentos de saida do Fortran."""
    IGREG = 2299161
    if julian >= IGREG:
        jalpha = _fint(((julian - 1867216) - 0.25) / 36524.25)
        ja = julian + 1 + jalpha - _fint(0.25 * jalpha)
    else:
        ja = julian
    jb = ja + 1524
    jc = _fint(6680. + ((jb - 2439870) - 122.1) / 365.25)
    jd = 365 * jc + _fint(0.25 * jc)
    je = _fint((jb - jd) / 30.6001)
    iday = jb - jd - _fint(30.6001 * je)
    mm = je - 1
    if mm > 12:
        mm -= 12
    iyyy = jc - 4715
    if mm > 2:
        iyyy -= 1
    if iyyy <= 0:
        iyyy -= 1
    return mm, iday, iyyy
