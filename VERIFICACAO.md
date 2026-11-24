# Verificação da tradução Python × Fortran original (MGB-SED v1.0, rio Ivaí)

## Diagnóstico curto
A versão Python reproduz as saídas do Fortran na simulação completa (01/01/1974 a 31/12/2020,
17.167 passos diários, 1.816 minibacias, Muskingum-Cunge, MUSLE + transporte com Meyer-Peter & Müller).
As diferenças estão no nível do arredondamento em precisão simples, não há deriva ao longo dos 47 anos,
e os arquivos de ajuste (`AJUSTE.fob`, `AJUSTESED.fob`) saem idênticos.

**Sobre os dados:** a verificação da série completa foi feita com a configuração calibrada do rio Ivaí usada na dissertação do autor. Os parâmetros calibrados e as séries simuladas não são publicados; este relatório traz apenas estatísticas de diferença entre Python e Fortran. O teste automático do repositório usa parâmetros genéricos e uma referência própria do Fortran (seção "Teste automático"). O script `verification/plots.py` gera as figuras Fortran × Python a partir de quaisquer duas rodadas.

## Método
1. **Referência**: saídas do Fortran fornecidas (`output_mgb-sed.rar`, `output_mgb-sed2.rar`), mesmas entradas (`input_mgb-sed.rar`).
2. **Execução Python**: `python -m mgbsed.main` com REAL = float32 (padrão), sem nenhum ajuste a partir das saídas.
3. **Comparação** (`verification/compare.py`): arquivo a arquivo, valor a valor, na precisão em que o Fortran grava (ex.: F15.3).
4. **Testes unitários** (`tests/test_units.py`): partes vetorizadas comparadas com traduções escalares linha a linha das rotinas Fortran.

Métricas da tabela: `idênticos` = % de valores iguais na precisão gravada; `erro rel. p99` = percentil 99 de |Py−F|/max(|F|, 10⁻³);
`>1%` = % de valores com |Py−F| > 1%·max(|F|,1); `viés` = diferença percentual das somas.

## Resultados (série completa)

| Arquivo | n valores | idênticos (%) | erro rel. p99 | >1% (%) | erro abs. máx. | correlação | viés (%) |
|---|---|---|---|---|---|---|---|
| VAZAO.HIG | 154,503 | 94.1 | 4.8e-06 | 0.0000 | 0.002 | 1.00000000 | +0.00000 |
| NOSOLO.txt | 137,336 | 98.8 | 8.8e-07 | 0.0000 | 0.0001 | 1.00000000 | +0.00000 |
| QPROP.HIG | 51,501 | 100.0 | 0.0e+00 | 0.0000 | 0.001 | 1.00000000 | +0.00003 |
| BACIAS.HIG | 120,169 | 99.8 | 0.0e+00 | 0.0000 | 1e-05 | 1.00000000 | +0.00000 |
| VAZAO_MINI.HIG | 31,175,272 | 99.6 | 0.0e+00 | 0.0000 | 0.003 | 1.00000000 | +0.00000 |
| QSUP.txt | 31,175,272 | 100.0 | 0.0e+00 | 0.0000 | 0.001 | 1.00000000 | +0.00000 |
| QCEL.txt | 31,175,272 | 100.0 | 0.0e+00 | 0.0000 | 0.001 | 1.00000000 | -0.00000 |
| MUSLE_TOTAL.txt | 31,175,272 | 97.9 | 1.9e-04 | 0.0014 | 0.214 | 1.00000000 | +0.00000 |
| MUSLE_areia.txt | 31,175,272 | 99.1 | 0.0e+00 | 0.0009 | 0.116 | 1.00000000 | +0.00001 |
| MUSLE_silte.txt | 31,175,272 | 99.4 | 0.0e+00 | 0.0005 | 0.052 | 1.00000000 | +0.00000 |
| MUSLE_argila.txt | 31,175,272 | 99.2 | 0.0e+00 | 0.0005 | 0.047 | 1.00000000 | +0.00000 |
| CONC_areia.txt | 31,175,272 | 90.7 | 4.1e-04 | 0.0127 | 602.5 | 0.99999999 | +0.00003 |
| CONC_silte.txt | 31,175,272 | 89.6 | 4.3e-04 | 0.0098 | 134.8 | 0.99999975 | +0.00000 |
| CONC_argila.txt | 31,175,272 | 87.7 | 4.6e-04 | 0.0120 | 199.4 | 0.99999974 | -0.00000 |
| DEPdin_RIO.txt | 31,175,272 | 89.9 | 1.0e-05 | 0.0007 | 9.739 | 1.00000000 | +0.00001 |
| MINI_areia.txt | 154,503 | 99.6 | 0.0e+00 | 0.0000 | 4 | 1.00000000 | +0.00000 |
| MINI_silte.txt | 154,503 | 97.8 | 3.8e-05 | 0.0000 | 0.001 | 1.00000000 | +0.00000 |
| MINI_argila.txt | 154,503 | 97.1 | 3.8e-05 | 0.0000 | 0.002 | 1.00000000 | -0.00000 |
| MINI_QSS.txt | 154,503 | 64.0 | 2.5e-05 | 0.0000 | 0.11 | 1.00000000 | +0.00000 |
| SEDRIO_SAIDA.txt | 5,448 | 16.7 | 1.3e-06 | 0.0000 | 2 | 1.00000000 | -0.00000 |
| SEDRIO_DEPOSITO.txt | 5,448 | 64.7 | 5.0e-04 | 0.0551 | 1.924 | 1.00000000 | +0.00000 |
| SEDRIO_EROSAO.txt | 5,448 | 80.9 | 1.2e-04 | 0.0184 | 1.921 | 1.00000000 | +0.00000 |
| BAL_PLAN.txt | 16,344 | 100.0 | 0.0e+00 | 0.0000 | 0 | — | +— |
| AJUSTESED.fob | 27 | 100.0 | 0.0e+00 | 0.0000 | 0 | 1.00000000 | -0.00000 |
| AJUSTE.fob | 66 | 100.0 | 0.0e+00 | 0.0000 | 0 | 1.00000000 | +0.00000 |
Outros conferidos: precipitação anual média por sub-bacia (log) idêntica nas 12 sub-bacias; balanço de massa de sedimentos
fecha com erro de −2×10⁻⁵ % (areia) e −4×10⁻⁶ % (silte, argila).

O erro relativo diário entre Python e Fortran, no máximo entre os postos, fica entre 10⁻⁶ e 10⁻⁴ durante todo o período, sem tendência (sem deriva em 47 anos).

## Interpretação das diferenças

**Verificado**
* Vazões, escoamento superficial e balanço vertical: > 94% dos valores idênticos; maior diferença = 1 a 3 unidades do último dígito gravado.
* QSS nos postos (`MINI_QSS.txt`, variável usada na avaliação): nenhuma diferença acima de 1%; métricas KGE/NSE/erro de volume idênticas.

**Inferência (fundamentada, não provada bit a bit)**
* As diferenças residuais vêm da ordem e da implementação de operações em float32 (pow/exp da libm do Intel × NumPy, ordem de somas).
  Elas são ampliadas onde a fórmula é mal condicionada: `DSUP` no `SOLO` é a diferença de termos ~10² mm que resulta em ~10⁻² mm.
* Os poucos casos isolados com diferença grande (até 0,013% dos valores em `CONC_*`, ex. 602 mg/L × 0) são **trocas de ramo por limiar**
  no `SED_PROPAG`: quando `CTS` da areia vira NaN (tensão abaixo da crítica, item 1 do README), o depósito `CSdep` inteiro é liberado;
  uma diferença de 1 ulp perto do limiar antecipa ou atrasa essa liberação em 1 dia numa minibacia. O efeito não se propaga
  (correlação > 0,9999997, viés de soma ≈ 0, postos sem diferença > 1%).
* O fato de as séries de areia coincidirem em 90,8% dos valores, inclusive nos eventos de liberação, indica que o Fortran original
  também trata comparações com NaN como falsas (o mesmo comportamento da tradução).

**Pendências (não verificáveis com os dados disponíveis)**
* Calibração de sedimentos (`ICALIB=3`) e previsão (`ICALIB=2`): traduzidos linha a linha,
  mas sem saída de referência do Fortran. Para verificá-los seria preciso rodar o executável Fortran nesses modos.
  (O modo inercial foi verificado; ver seção seguinte.)
* Opções de capacidade de transporte 2–6 (`flag_TC`): traduzidas, não exercitadas.
* Bit a bit: não é alcançável sem o mesmo compilador/libm; também não é necessário para os usos previstos.

## Modo inercial (`INERTIAL = 1`)

Mesmas entradas, só com `INERTIAL = 1` no `infoMGB.sim`; saídas do Fortran fornecidas por Ivan (03/11/2026). Série completa.

| Arquivo | n valores | idênticos (%) | erro rel. p99 | >1% (%) | erro abs. máx. | correlação | viés (%) |
|---|---|---|---|---|---|---|---|
| VAZAO.HIG | 154,503 | 98.2 | 1.5e-06 | 0.0000 | 0.001 | 1.00000000 | +0.00000 |
| NOSOLO.txt | 137,336 | 98.8 | 8.8e-07 | 0.0000 | 0.0001 | 1.00000000 | +0.00000 |
| QPROP.HIG | — | — | — | — | — | — | ver nota |
| BACIAS.HIG | 120,169 | 99.8 | 0.0e+00 | 0.0000 | 1e-05 | 1.00000000 | +0.00000 |
| VAZAO_MINI.HIG | 31,175,272 | 99.8 | 0.0e+00 | 0.0000 | 0.002 | 1.00000000 | +0.00000 |
| QSUP.txt | 31,175,272 | 100.0 | 0.0e+00 | 0.0000 | 0.001 | 1.00000000 | +0.00000 |
| QCEL.txt | 31,175,272 | 100.0 | 0.0e+00 | 0.0000 | 0.001 | 1.00000000 | +0.00000 |
| MUSLE_TOTAL.txt | 31,175,272 | 97.9 | 1.9e-04 | 0.0014 | 0.214 | 1.00000000 | +0.00000 |
| MUSLE_areia.txt | 31,175,272 | 99.1 | 0.0e+00 | 0.0009 | 0.116 | 1.00000000 | +0.00001 |
| MUSLE_silte.txt | 31,175,272 | 99.4 | 0.0e+00 | 0.0005 | 0.052 | 1.00000000 | +0.00000 |
| MUSLE_argila.txt | 31,175,272 | 99.2 | 0.0e+00 | 0.0005 | 0.047 | 1.00000000 | +0.00000 |
| CONC_areia.txt | 31,175,272 | 92.8 | 2.7e-04 | 0.0088 | 3360 | 0.99999998 | +0.00007 |
| CONC_silte.txt | 31,175,272 | 89.6 | 4.2e-04 | 0.0096 | 92.58 | 0.99999992 | +0.00001 |
| CONC_argila.txt | 31,175,272 | 87.8 | 4.5e-04 | 0.0117 | 138.3 | 0.99999991 | +0.00001 |
| DEPdin_RIO.txt | 31,175,272 | 83.1 | 1.1e-05 | 0.0005 | 7.022 | 1.00000000 | -0.00000 |
| MINI_areia.txt | 154,503 | 99.7 | 0.0e+00 | 0.0000 | 1 | 1.00000000 | +0.00000 |
| MINI_silte.txt | 154,503 | 98.0 | 3.7e-05 | 0.0000 | 0.001 | 1.00000000 | -0.00000 |
| MINI_argila.txt | 154,503 | 97.6 | 3.6e-05 | 0.0000 | 0.002 | 1.00000000 | -0.00001 |
| MINI_QSS.txt | 154,503 | 64.9 | 2.7e-05 | 0.0000 | 0.078 | 1.00000000 | +0.00000 |
| SEDRIO_SAIDA.txt | 5,448 | 17.1 | 1.1e-06 | 0.0000 | 6.18 | 1.00000000 | -0.00000 |
| SEDRIO_DEPOSITO.txt | 5,448 | 57.7 | 4.4e-04 | 0.0367 | 6.18 | 1.00000000 | -0.00000 |
| SEDRIO_EROSAO.txt | 5,448 | 73.1 | 2.7e-05 | 0.0000 | 6.18 | 1.00000000 | -0.00001 |
| BAL_PLAN.txt | 16,344 | 89.7 | 1.2e-06 | 0.0000 | 0.038 | 1.00000000 | -0.00001 |
| AJUSTESED.fob | 27 | 100.0 | 0.0e+00 | 0.0000 | 0 | 1.00000000 | -0.00000 |
| AJUSTE.fob | 66 | 100.0 | 0.0e+00 | 0.0000 | 0 | 1.00000000 | +0.00000 |

**Verificado:** vazão (98,2% idênticos em `VAZAO.HIG`), QSS nos postos (nenhum valor > 1%), trocas com a planície
(`BAL_PLAN.txt`: erro máx. 0,038 em valores da ordem de 10⁴), `AJUSTESED.fob` e `AJUSTE.fob` idênticos.

**Nota sobre `QPROP.HIG` (não é erro de tradução):** no modo inercial `REDE` não é chamada, então `PJB2/PJI2/PJS2` nunca
recebem valor. No Fortran ficam com o conteúdo não inicializado da memória (o arquivo do Fortran tem vazões de
−5×10⁷ m³/s, com as três colunas em proporção exata 1:2:3 e 1.282 linhas com campos colados); no Python começam em zero e o arquivo sai com zeros.
O arquivo não tem significado físico nesse modo em nenhuma das versões.

## Desempenho
Python (NumPy, 1 núcleo): ~28 min (Muskingum-Cunge) e ~94 min (inercial, ~250 sub-passos por dia) para a série completa, incluindo ~5 GB de arquivos texto. O Fortran levou 6,3 min.
A opção `--no-big-outputs` evita gravar os ~5 GB de texto (tempo não medido).

## Teste automático (no repositório)

`tests/test_regressao.py` compara o Python diretamente com o executável Fortran original. Os dados (`tests/dados_teste/`) usam a topologia, os solos, a chuva e a vazão observada (ANA) do Ivaí a partir de 01/01/1974, com **parâmetros genéricos, não calibrados**, MUSLE com os valores padrão (α = 11,8; β = 0,56; AuxTKS = 1,0) e **observações de sedimento sintéticas**. O Fortran (MGB-SED v1.0, Intel Fortran ifx 2026.1.1, Debug x64, com o bloco "VALORES PADRÃO" ativo no `SED_BACIA.F90`) rodou 400 dias nas três configurações: Muskingum-Cunge, inercial e modo CSS.

O teste roda o Python por 60 dias e compara 19 arquivos diários com as 60 primeiras linhas do Fortran (critérios por arquivo: nenhum NaN a mais ou a menos; correlação ≥ 0,99999; no máximo 0,1 % dos valores com diferença > 1 %·max(|Ref|, 1); erro absoluto máximo dentro de uma tolerância declarada). Também roda 400 dias (Muskingum-Cunge e CSS) e exige `AJUSTE.fob` e `AJUSTESED.fob` idênticos ao Fortran e `SEDRIO_*` com diferença ≤ 0,05 t.

Comparação dos 400 dias completos, Python × Fortran, com os mesmos critérios:

| Configuração | Arquivos diários | Valores idênticos (na precisão gravada) | Maior diferença absoluta | Valores com diferença > 1 % | `AJUSTE.fob` e `AJUSTESED.fob` |
|---|---|---|---|---|---|
| Muskingum-Cunge | 19 | 95,4–100 % | 0,005 (DEPdin_RIO, t) | 0 | idênticos |
| Inercial | 18 (sem `QPROP`) | 97,7–100 % | 0,024 mg/L (CONC_areia) | 0 | idênticos |
| CSS (`FLAG_SEDVAR=1`) | 19 (com `MINI_CSS`) | 95,4–100 % | 0,005 (DEPdin_RIO, t) | 0 | idênticos |

Os acumulados de fim de rodada (`SEDRIO_*`) diferem em no máximo 0,006 t e o `BAL_PLAN.txt` em no máximo 0,0001.

**Sensibilidade do teste:** três alterações pequenas introduzidas de propósito, uma de cada vez, fizeram o teste falhar: coeficiente do fluxo de calor no solo de 0,38 para 0,37 (9 testes falharam); constante de Muskingum-Cunge de 3,125 para 3,0 (13); β padrão da MUSLE de 0,56 para 0,57 (27).

**Limites:** 400 dias não exercitam o comportamento de longo prazo, coberto pela verificação da série completa acima. Com menos de 2 observações de sedimento a partir do dia 365, o Fortran compilado em Debug aborta em `SED_Fob_Calib` (índice 0 na curva de permanência) e o Python grava métricas nulas com um aviso.

## Windows × Linux

A rodada completa feita no Windows (PC do autor, Python 3.12) e a feita no Linux (Python 3.13, NumPy 2.3) geram `VAZAO.HIG`, `MINI_silte.txt` e `MINI_argila.txt` idênticos na precisão gravada (100 % dos valores iguais). As diferenças em relação ao Fortran, portanto, não dependem do sistema operacional.

## Modo CSS (`FLAG_SEDVAR = 1`)

No Fortran, o `FLAG_SEDVAR` só escolhe a variável gravada e comparada com as observações (`MODELO.f90`, `SIMULA.f90`): CSS = (silte + argila)·10⁶ mg/L, sem efeito na simulação. Na rodada completa do Windows, `MINI_CSS.txt` difere de `MINI_silte.txt` + `MINI_argila.txt` em no máximo 0,001 mg/L (arredondamento da gravação). O teste automático compara o `MINI_CSS.txt` do Python diretamente com o do Fortran (60 dias, diferença máxima 0,001 mg/L) e exige o `AJUSTESED.fob` do modo CSS idêntico numa rodada de 400 dias. O `AJUSTESED.fob` do modo CSS da rodada original do autor (24/11/2026) coincide, nos 9 postos e com 3 casas decimais (KGE, NSE e erro de volume), com as métricas calculadas sobre a saída do Python.
