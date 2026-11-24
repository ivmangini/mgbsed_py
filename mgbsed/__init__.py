"""MGB-SED v1.0 - traducao para Python do modelo hidrologico-sedimentologico
MGB-IPH + modulo de sedimentos (versao Fortran de jul/2019 - nov/2020).

Modulos:
    fortran_io  - E/S com semantica Fortran, julday/caldat
    state       - objeto Model (variaveis dos modulos Fortran)
    io_input    - leitura de arquivos de entrada
    hydrology   - PARCEL, PARCUNGE, CONDINIC, LECLIMA, CELULA (RADIACAO, EVAPO, SOLO)
    routing     - REDE + MUSK (Muskingum-Cunge)
    inertial    - modelo inercial (flood_*)
    sediment    - SED_INICIAL, SED_BACIA, SED_REDE, SED_PROPAG, SED_CT
    objective   - FOBJ, SED_Fob_Calib
    simulation  - MODELO, SIMULA
    calibration - CALIBRA (vazia no original), CALIBRA_SED (MOCOM-UA)
    forecast    - PREVISAO, ATUALIZA
    main        - programa principal (1main.f90)
"""
__version__ = "1.0.0"
