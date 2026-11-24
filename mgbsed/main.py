"""
Programa principal (traducao de 1main.f90 - PROGRAM MGB_Inercial).

Uso:
    python -m mgbsed.main --input ./Input --output ./Output [--float64] [--no-big-outputs]
"""
from __future__ import annotations

import argparse
import os
import time

import numpy as np

from . import io_input as io
from . import hydrology as hy
from . import routing as rt
from . import sediment as sed
from .fortran_io import julday
from .state import Model


def setup(input_dir, output_dir, real=np.float32, write_big_outputs=True, verbose=True,
          nt_override=None):
    """Leitura das entradas e pre-calculos (1main.f90 ate a selecao do modo).

    ``nt_override``: limita o numero de passos de tempo (testes rapidos).
    """
    os.makedirs(output_dir, exist_ok=True)
    M = Model(input_dir, output_dir, real=real, write_big_outputs=write_big_outputs, verbose=verbose)
    M._logf = open(os.path.join(output_dir, "LOG_MGB.txt"), "w", newline="\r\n", encoding="latin-1")
    M.log("Modelo Hidrologico de Grandes Bacias: MGB-IPH  (versao Python do MGB-SED v1.0)")
    M.log("Reading infoMGB.sim file...")
    io.lefix(M)
    if nt_override is not None:
        M.NT = int(nt_override)
    M.log("Allocating variables...")
    io.alloca_vars(M)
    io.alloca_varsinerc(M)
    M.IDINI = julday(M.IMES, M.IDIA, M.IANO)
    M.log("Reading fixed parameters file...")
    io.levar(M)
    M.log("Reading calibrated parameters file...")
    io.leuso(M)
    M.log("Reading model topology (MINI.GTP) file...")
    io.lecell(M)
    M.log("Reading observed flows, substituted flows, and climate files...")
    io.leclimed(M)
    M.log("CLIMED OK")
    io.arqclisub(M)
    M.log("CLISUB OK")
    io.leqobs(M)
    M.log("QOBS OK")
    if M.NUMSUBST > 0:
        io.lesubst(M)
    M.log("QSUBS OK")
    if M.IND_FLOODROUTING == 1:
        M.HDFLAG[:] = 1
        io.lefloodpln(M)
    else:
        M.HDFLAG[:] = 0
    M.hdFLAG0 = int(M.HDFLAG.sum())
    M.log("Calculating river geometry informations...")
    hy.parcel(M)
    M.log("PARCEL Routine OK")
    hy.parcunge(M)
    M.log("PARCUNGE Routine OK")
    if M.HDFLAG[0] == 1:
        from . import inertial as inr
        inr.flood_topo(M)
        M.log("FLOOD TOPOLOGY OK")
        inr.flood_tab(M)
        M.log("FLOOD TAB OK")
    rt.topo_levels(M)
    if M.hdFLAG0 == 0:
        rt.prepare_routing(M)
    # QSAUX e usado por CELULA mesmo sem sedimentos
    M.QSAUX = np.zeros((M.NC, M.NU + 2))
    if M.SED_INDEX == 1:
        sed.alloca_vars_sed(M)
        io.leobsed(M)
        M.log("SED OBS OK")
        M.xyz = 1
        if M.ICALIB == 0:
            io.le_par_musle_sed(M)
            io.sed_param(M)
            M.log("Reading fixed MUSLE parameters file...")
    M.log("Opening input/output files...")
    M.rainfile = io.RainFile(io._path(M, "CHUVAbin.pbi"), M.NC)
    M.log("Input data reading process ended")
    M.log(f"{'Number of catchments:':>30}{M.NC:8d}")
    M.log(f"{'Number of subbasins:':>30}{M.NB:8d}")
    M.log(f"{'Number of HRUs:':>30}{M.NU:8d}")
    M.log(f"{'Start date (dd mm yyyy hh):':>30}{M.IDIA:5d}{M.IMES:5d}{M.IANO:5d}{M.HORAINI:5d}")
    M.log(f"{'Time step (s):':>30}{M.DTP:16.2f}")
    M.log(f"{'Number of time steps:':>30}{M.NT:8d}")
    M.log(f"Number of observation gauges (discharge): {M.NOBS}")
    M.log(f"Number of substitution gauges (discharge): {M.NUMSUBST}")
    for k in range(M.AUXSS):
        M.log(f"Number of observation gauges (sediment): {M.NOBSS[k]}")
    M.log(f"{'Number of climate gauges:':>30}{M.NCLI:8d}")
    M.log("ICALIB (simulation=0; autocalibration=1; forecasting=2; Sed autocalibration=3)")
    M.log(f"{'ICALIB :':>30}{M.ICALIB:8d}")
    M.log(f"{'Reservoir Flag (0=no; 1=yes):':>40}{M.RESERV_INDEX:8d}")
    M.log(f"{'Inertial routing Flag (0=no; 1=yes):':>40}{M.IND_FLOODROUTING:8d}")
    M.log(f"{'Sediment module Flag (0=no; 1=yes):':>40}{M.SED_INDEX:8d}")
    M.log(f"{'Alpha for Inertial routing:':>40}{M.alpha:10.2f}")
    M.log(f"Precisao REAL: {np.dtype(real).name}")
    return M


def run(M, progress=True):
    t0 = time.time()
    if M.ICALIB == 0:
        M.log("SIMULATION MODE")
        from .simulation import simula
        simula(M, progress=progress)
    elif M.ICALIB == 1:
        M.log("AUTOMATIC CALIBRATION MODE - HIDROLOGICAL MODEL")
        from . import calibration as cal
        io.lecalib(M)
        cal.calibra(M)
    elif M.ICALIB == 2:
        M.log("FORECASTING MODE")
        from .forecast import previsao
        previsao(M)
    elif M.ICALIB == 3:
        M.log("AUTOMATIC CALIBRATION MODE - SEDIMENT MODEL")
        from . import calibration as cal
        io.lecalib_sed(M)
        cal.calibra_sed(M)
    else:
        raise RuntimeError(" ERROR: UNKNOWM ICALIB!!!")
    M.log("MGB ENDED SUCCESSFULLY")
    M.log("TOTAL TIME:", (time.time() - t0) / 60., "MINUTES")
    M._logf.close()
    M._logf = None
    return M


def cli():
    ap = argparse.ArgumentParser(description="MGB-SED v1.0 (traducao Python)")
    ap.add_argument("--input", default="./Input")
    ap.add_argument("--output", default="./Output")
    ap.add_argument("--float64", action="store_true", help="usa REAL=float64 (padrao: float32, como no Fortran)")
    ap.add_argument("--no-big-outputs", action="store_true",
                    help="nao grava MUSLE_*, QSUP, QCEL, CONC_*, DEPdin_RIO, VAZAO_MINI (arquivos de ~470 MB)")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    M = setup(a.input, a.output, real=np.float64 if a.float64 else np.float32,
              write_big_outputs=not a.no_big_outputs, verbose=not a.quiet)
    run(M, progress=not a.quiet)


if __name__ == "__main__":
    cli()
