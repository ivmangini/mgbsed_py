# Changelog

## 1.0.0 (2026-11-24)
Primeira versão publicada.

- Tradução para Python (NumPy) de todas as 71 rotinas Fortran do MGB-SED v1.0.
- Verificação contra as saídas do Fortran original, rio Ivaí, 1974–2020: propagação Muskingum-Cunge e inercial (ver `VERIFICACAO.md`).
- Testes automáticos: unitários (partes vetorizadas × tradução escalar) e regressão em 60 dias, nos dois modos de propagação.
- Comportamentos do código original reproduzidos de propósito e documentados no README (seção "Comportamentos herdados").
- Diferença de interface: α, β e AuxTKS da MUSLE por sub-bacia são lidos de `PAR_MUSLE_SED.txt` (opcional) em vez de escritos no código; sem o arquivo, valores padrão (11,8; 0,56; 1,0).
- Teste de regressão contra o executável Fortran original, com parâmetros genéricos e observações sintéticas (`tests/dados_teste/`): 60 dias de saídas diárias nos modos Muskingum-Cunge, inercial e CSS; métricas de ajuste idênticas numa rodada de 400 dias.
