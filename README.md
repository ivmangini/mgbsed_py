# MGB-SED v1.0 em Python

[![testes](https://github.com/ivmangini/mgbsed_py/actions/workflows/testes.yml/badge.svg)](https://github.com/ivmangini/mgbsed_py/actions/workflows/testes.yml)

Tradução para Python (NumPy) do modelo hidrossedimentológico **MGB-SED v1.0** (MGB-IPH + módulo de sedimentos, IPH/UFRGS), escrito originalmente em Fortran. A tradução foi verificada contra as saídas do Fortran original na bacia do rio Ivaí (PR), 1974–2020: 17.167 dias, 1.816 minibacias, 12 sub-bacias e 15 URHs.

- **Fiel ao Fortran:** mesmos arquivos de entrada e de saída, mesmos formatos e precisão simples (`float32`), como o `REAL` do Fortran.
- **Verificado:** na série completa, a QSS nos postos não tem nenhum valor com diferença acima de 1 % e o `AJUSTESED.fob` sai idêntico (`VERIFICACAO.md`).
- **Testado a cada mudança:** o GitHub Actions roda os testes unitários e a regressão contra o Fortran (60 dias, Muskingum-Cunge e inercial) em Windows e Linux, Python 3.10 a 3.12.

## Instalação

Requisitos: Python ≥ 3.10 e NumPy ≥ 1.26.

```bash
git clone https://github.com/ivmangini/mgbsed_py.git
cd mgbsed_py
pip install .            # ou: pip install ".[test]" para rodar os testes
```

## Uso

Igual ao Fortran: uma pasta de entrada com os arquivos do MGB-SED (`infoMGB.sim`, `MINI.gtp`, `CHUVABIN.PBI`, ...) e uma pasta de saída.

```bash
mgbsed --input ./input --output ./output
#   --float64          REAL em float64 (padrão: float32, igual ao Fortran)
#   --no-big-outputs   não grava os arquivos de ~470 MB (MUSLE_*, QSUP, QCEL, CONC_*, DEPdin_RIO, VAZAO_MINI)
```

Dentro do Python (Spyder, Jupyter), com o estado completo do modelo disponível ao final:

```python
import numpy as np
from mgbsed.main import setup, run
M = setup("input", "output", real=np.float32)   # lê entradas e faz os pré-cálculos
run(M)                                          # executa o modo definido por ICALIB no infoMGB.sim
M.QRG          # vazão nos postos de VAZAO.HIG (n_postos, NT)
M.SEDstore[0]  # QSS nos postos (t/dia)
```

Um script pronto para o Spyder está em `examples/rodar_spyder.py`.

### Parâmetros da MUSLE por sub-bacia (`PAR_MUSLE_SED.txt`)

No Fortran original, α, β e o fator do tempo de retardo (`AuxTKS`) de cada sub-bacia são escritos no código-fonte (`SED_BACIA.F90`). Aqui eles são lidos do arquivo opcional `PAR_MUSLE_SED.txt` na pasta de entrada:

```
! IB    ALF    BET  AuxTKS
   1  11.80   0.56   1.00
   2  11.80   0.56   1.00
 ...  (uma linha por sub-bacia)
```

Sem esse arquivo, todas as sub-bacias usam os valores padrão da MUSLE (α = 11,8; β = 0,56; AuxTKS = 1,0), os mesmos do bloco "VALORES PADRÃO" comentado no `SED_BACIA.F90`. O `LOG_MGB.txt` informa qual dos dois foi usado.

Tempo da série completa do Ivaí em 1 núcleo: ~27 min com Muskingum-Cunge e ~94 min com o inercial (o Fortran leva 6,3 min com Muskingum-Cunge).

## Testes

```bash
pip install ".[test]"
pytest
```

| Teste | O que verifica |
|---|---|
| `tests/test_units.py` | Calendário, formatos Fortran, `SOLO` vetorizado = tradução escalar linha a linha (bit a bit), estado de `F4` em `EVAPO`, Muskingum-Cunge matricial = recorrência original |
| `tests/test_regressao.py` | Simulação de 60 dias do rio Ivaí, com Muskingum-Cunge e com o inercial: 19 arquivos diários comparados com as saídas do Fortran original, balanço de massa de sedimentos, modo de saída CSS (`FLAG_SEDVAR=1`) e, numa rodada de 400 dias, métricas de ajuste (`AJUSTE.fob`, `AJUSTESED.fob`) idênticas às do Fortran |

Os dados de teste (`tests/dados_teste/`, ~8 MB) usam a topologia, os solos, a chuva e a vazão observada (ANA) do rio Ivaí a partir de 01/01/1974, mas **com parâmetros genéricos, não calibrados** (`PARUSO.cal`, `PARUSO_MUSLE.txt`, `PARTEXT_MUSLE.txt`; MUSLE com os valores padrão) e **observações de sedimento sintéticas** (uma a cada 5 dias por posto). A referência são as saídas do executável Fortran original rodado com essas mesmas entradas por 400 dias:

- `ref_fortran_mc/`, `ref_fortran_inercial/`, `ref_fortran_css/`: 60 primeiras linhas dos arquivos diários (o modelo só depende do passado, então são as mesmas de uma rodada de 60 dias);
- `ref_fortran_400dias_mc/`, `ref_fortran_400dias_css/`: `AJUSTE.fob`, `AJUSTESED.fob` e `SEDRIO_*` dos 400 dias. São 400 dias porque o Fortran só calcula as métricas de sedimento a partir do dia 365.

O teste exige `AJUSTE.fob` e `AJUSTESED.fob` idênticos ao Fortran. Nos 400 dias completos, todos os arquivos diários também passam nos mesmos critérios (ver `VERIFICACAO.md`).

Para conferir que o teste realmente detecta erros, foram introduzidas de propósito três alterações pequenas, uma de cada vez: coeficiente do fluxo de calor no solo 0,38 → 0,37 (9 testes falharam); constante de Muskingum-Cunge 3,125 → 3,0 (13); β padrão da MUSLE 0,56 → 0,57 (27).

A verificação da série completa (47 anos) não roda no CI, porque exige os dados completos. Os scripts estão em `verification/` (`compare.py`, `plots.py`, `metricas_qss.py`), e os resultados em `VERIFICACAO.md` e `verification/resultados/`.

## Estrutura e correspondência Fortran → Python

| Python | Rotinas Fortran |
|---|---|
| `main.py` | `1main.f90` |
| `state.py` | módulos `VARS_MAIN`, `VARS_CALIB`, `VARS_INERC`, `SED_VARS` (atributos do objeto `Model`) |
| `fortran_io.py` | semântica de E/S Fortran (formatos fixos, list-directed, `Fw.d`, `ESw.d`), `julday.f90`, `caldat.f90` |
| `io_input.py` | `LEFIX`, `ALLOCA_VARS`, `ALLOCA_VARSINERC`, `LEVAR`, `LEUSO`, `LECELL`, `LECLIMED`, `ARQCLISUB`, `LEQOBS`, `LESUBST`, `LEFLOODPLN`, `LEOBSed`, `SED_PARAM`, `LeCalib`, `LeCalib_SED`, `LECHUVA`, `LEPREV` |
| `hydrology.py` | `PARCEL`, `INTECLIM`, `REGION`, `PARCUNGE`, `CondInic`, `LECLIMA`, `CELULA`, `RADIACAO`, `EVAPO`, `SOLO` |
| `routing.py` | `REDE`, `MUSK` |
| `inertial.py` | `flood_TOPO`, `flood_TAB`, `flood_ini`, `flood_inercial`, `flood_timestep`, `flood_discharge`, `flood_continuity`, `flood_write`, `hunt`, `FINT`, ramo inercial de `SED_REDE` |
| `sediment.py` | `ALLOCA_VARS_SED`, `SED_INICIAL`, `FALLVELOCITY`, `SED_BACIA`, `SED_REDE`, `SED_PROPAG`, `SED_CT`, `NEWTON_H` |
| `objective.py` | `FOBJ`, `SED_Fob_Calib` |
| `simulation.py` | `MODELO`, `SIMULA` |
| `calibration.py` | `CALIBRA` (vazia no original), `CalibParam`, `SEMENTE`, `ran1`, `ALLOCA_CALIB`, `CALIBRA_SED`, `CalibParam_SED`, `SED_CALIB_MODELO`, `SED_CALIB_BACIA`, `SED_CALIB_REDE` |
| `forecast.py` | `PREVISAO`, `ATUALIZA` |
| `fmath.py` | (novo) pow/exp/log em float32 com arredondamento único, para aproximar a libm do Intel |

Todos os 71 arquivos `.f90/.F90` foram traduzidos. `TEMPO`/`SECNDS` (cronômetro) viraram `time.time()`.

### Convenções
* Índices de minibacia/sub-bacia/URH são 0-based nos arrays; códigos gravados em arquivo continuam 1-based.
* Precisão: por padrão `float32` em tudo que é `REAL` no Fortran e `float64` no que é `REAL(8)`
  (`PSC`, `QSAUX`, `CARGM`, `DEPTREC`, `EROSTREC`, `iENTRA`, `iSAI`, variáveis do modelo inercial).
  As promoções de tipo do NumPy (float32 ⊕ float64 → float64) reproduzem as do Fortran.
* Formatos de saída idênticos (larguras, cabeçalhos, `CRLF`), para que pós-processadores existentes continuem funcionando.

### Vetorização (o que mudou na forma, não no cálculo)
* `CELULA/RADIACAO/EVAPO/SOLO`: vetorizados em matrizes (NC, NU). Somas por minibacia feitas num laço sobre URH, na ordem do Fortran.
* `SOLO`: o laço `DO ... ENDDO` de ajuste do balanço virou um laço com máscaras que repete a mesma sequência de ramos (teste unitário compara com a tradução escalar linha a linha: resultados bit a bit iguais).
* `REDE/MUSK`: minibacias agrupadas por nível topológico (85 níveis no Ivaí). Sem o `MAX(·,0)`, Muskingum-Cunge é linear; a grade completa é obtida por matriz pré-calculada e, se qualquer ponto der negativo, a minibacia é recalculada pela recorrência original. O resultado é o da recorrência original.
* `SED_REDE/SED_PROPAG`: a única dependência entre minibacias é `CARGM(JUS)`; processamento por nível topológico preserva essa dependência exatamente.

## Comportamentos herdados do código original (reproduzidos, não corrigidos)

A tradução é fiel: os itens abaixo estão no Fortran e foram mantidos para que as saídas coincidam.
Eles são relevantes para a interpretação dos resultados e para a etapa de IA.

**Afetam resultados físicos**
1. **`SED_CT`, Meyer-Peter & Müller (`flag_TC=1`, usado no Ivaí)**: `Kr=26/(D90**(1/6))` e `(G/9.82)**(1/3)` usam divisão inteira (`1/6=0`, `1/3=0`), então `Kr=26` sempre e o denominador vale 1. Quando a tensão está abaixo da crítica a base de `**1.5` é negativa → `CTS=NaN`. Em `SED_PROPAG`, `CSJ2 >= NaN` é falso e o código entra no ramo de **erosão**, onde `EROSaux=NaN`, `NaN <= CSdep` é falso, e **todo o depósito de areia `CSdep` é remobilizado de uma vez**. Ou seja: abaixo do limiar de movimento, a areia depositada é integralmente erodida. O Fortran original faz o mesmo: as séries de areia coincidem (`VERIFICACAO.md`).
2. **`SED_CT`** (outras opções): `(2/3)`, `(4/3)`, `(1/3)` também são divisões inteiras em Rottner, Schoklitsch e Ackers-White.
3. **`SED_BACIA`**: no Fortran, α, β da MUSLE e o fator de `TKS` por sub-bacia estão **fixos no código-fonte** quando `ICALIB=0` (sub-bacias 1–12), e `ParCalibSED.txt` não é usado na simulação. **Esta é a única diferença de interface da versão Python:** os valores vêm de `PAR_MUSLE_SED.txt` ou, sem ele, dos valores padrão (ver "Uso").
4. **`SED_PROPAG`**: o termo de variação do leito (`EQ3` da areia) é calculado e depois zerado; `HRIO1/HRIO2` são atualizados mas não influenciam nada.
5. **`EVAPO`**: quando `W ≤ WPM`, `RSX = 999999999·F4` usa o `F4` "velho" (variável de módulo da última URH que o atribuiu, inclusive de outro dia).
6. **`CondInic`**: `TK` é `INTEGER` (truncamento de `CB·3600` e `CI·TIND`) na condição inicial.
7. **`SED_BACIA` x `SED_CALIB_BACIA`**: a vazão de pico usa `/360` na simulação e `/3.6` na calibração.

**Afetam arquivos de saída/diagnóstico**
8. **`AJUSTE.fob`**: `SED_Fob_Calib` sobrescreve `R2`, `R2L`, `ERRV` (vazão) com −100/−100/1000 antes da gravação — confirmado no `AJUSTE.fob` original. A versão Python reproduz isso e grava em separado **`AJUSTE_HIDRO_PY.fob`** (arquivo novo) com os valores de NSE, NSE-log e erro de volume de vazão calculados por `FOBJ`.
9. **`QCEL.txt`** usa `QCEL2` do passo anterior (`SED_BACIA` roda antes de `REDE`).
10. **`DEPdin_RIO.txt`** contém o depósito **acumulado** de areia `CSdep` (t), não um valor diário; não tem cabeçalho.
11. `SED_Fob_Calib`: `auxcont` não é zerado entre postos (curva de permanência acumulada); `FOBJ` começa no dia 360 e `SED_Fob_Calib` no 365.
12. A unidade 25001 (NOSOLO da minibacia 2) nunca é aberta: o Fortran grava `fort.25001`; o Python também.
13. `AREA_ACU.txt` (3,4 GB) não é lido por nenhuma rotina.
13b. Modo inercial: `QPROP.HIG` usa `PJB2/PJI2/PJS2`, que só são calculados em `REDE` (não chamada nesse modo). No Fortran saem valores de memória não inicializada; no Python, zeros. O arquivo não tem significado nesse modo.

**Modos não usados na configuração do Ivaí**
14. `CALIBRA` (calibração hidrológica) é uma sub-rotina vazia.
15. `CALIBRA_SED`: a função objetivo da população inicial nunca é avaliada (`FO=0`), então nenhum ponto é dominado e o algoritmo termina na 1ª iteração. Além disso `SED_CALIB_REDE` lê `Q_SC`, `HmJx_SC` etc., que só são preenchidos no modo inercial, e usa `MINIMONT` (só existe no modo inercial).
16. `PREVISAO` usa caminhos fixos `c:\MGB-hora\...` e nomes `ARQHID` nunca lidos; foram trocados pelas pastas de entrada/saída (documentado em `forecast.py`). A previsão só dispara quando a hora de `Previsao.meta` coincide com `MOD(IT,24)+HORAINI` (lógica pensada para passo horário; com passo diário a "hora" é um contador sem significado).


## Limitações

- **Não verificados contra o Fortran** (não havia saídas de referência): calibração automática de sedimentos (`ICALIB=3`), previsão (`ICALIB=2`), opções de capacidade de transporte 2 a 6 e comparação por QS total (`FLAG_SEDVAR` 3). No modo CSS (`FLAG_SEDVAR` 1), o teste confere que `MINI_CSS.txt` = silte + argila do Fortran; as métricas do `AJUSTESED.fob` nesse modo não foram comparadas com uma rodada do Fortran. Esses caminhos foram traduzidos linha a linha; alguns foram só executados por poucos dias, sem erro.
- **Equivalência não é bit a bit:** as diferenças vêm do arredondamento em `float32` entre o Intel Fortran e o NumPy (ver `VERIFICACAO.md`).
- **Correção física:** a tradução reproduz o modelo como está, inclusive os comportamentos listados acima. Equivalência com o Fortran não significa correção física.

## Como citar

Cite este software (`CITATION.cff`; botão "Cite this repository" no GitHub) e os trabalhos que descrevem o modelo:

- Buarque, D. C. (2015). *Simulação da geração e do transporte de sedimentos em grandes bacias: estudo de caso do Rio Madeira* [Tese de doutorado, Instituto de Pesquisas Hidráulicas, Universidade Federal do Rio Grande do Sul].
- Collischonn, W., Allasia, D. G., Silva, B. C., & Tucci, C. E. M. (2007). The MGB-IPH model for large-scale rainfall–runoff modelling. *Hydrological Sciences Journal*, 52(5), 878–895. https://doi.org/10.1623/hysj.52.5.878

## Créditos

Código Fortran original: MGB-SED v1.0, Instituto de Pesquisas Hidráulicas (IPH/UFRGS); o código-fonte registra alterações de Hugo Fagundes no módulo de sedimentos. Plugin do MGB-SED para QGIS: [HGE-IPH/mgb-sed-plugin](https://github.com/HGE-IPH/mgb-sed-plugin) (licença MIT).

Tradução para Python e verificação: Ivan Mangini.

## Licença

MIT (`LICENSE`).
