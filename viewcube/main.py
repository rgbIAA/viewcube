"""Módulo principal de ViewCube para visualización de datos astronómicos."""
############################################################################
#                            VIEW-CUBE/RSS                                 #
#                               PYTHON 3                                   #
#                                                                          #
# RGB@IAA ---> Last Change: 2024/11/30                                     #
############################################################################

"""
 Author: Ruben Garcia-Benito (RGB)
"""
import argparse
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.rcsetup

from viewcube import config as vc
from viewcube import version
from viewcube.cubeviewer_old import CubeViewer
from viewcube.viewrss import RSSViewer

VERSION = "0.0.7"


def parse_arguments() -> argparse.Namespace:
    """Configura y analiza los argumentos de línea de comandos."""
    list_backends = matplotlib.rcsetup.interactive_bk
    slist_backends = " | ".join(list_backends)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Argumentos principales
    parser.add_argument("name", type=str, help="FITS file", nargs="*")

    # Configuración de visualización
    parser.add_argument(
        "-b",
        type=str,
        help=f"Matplotlib backend. Available: {slist_backends}",
        default=matplotlib.rcParams["backend"]
    )
    parser.add_argument(
        "-y",
        type=str,
        help="Plot style (comma-separated): 'dark_background,seaborn-ticks'"
    )

    # Parámetros de procesamiento
    parser.add_argument("--data", type=int, help="DATA extension")
    parser.add_argument("--error", type=int, help="ERROR extension")
    parser.add_argument("--flag", type=int, help="FLAG/MASK extension")
    parser.add_argument("--header", type=int, default=0, help="HEADER extension")
    parser.add_argument(
        "-a",
        type=float,
        help="Angle to rotate position table (RSS only)"
    )
    parser.add_argument(
        "-p",
        type=str,
        help="External position table for RSS Viewer"
    )
    parser.add_argument("-s", type=int, help="Spectral dimension")

    # Opciones de comportamiento
    parser.add_argument("-v", help="Print version", action="store_true")
    parser.add_argument(
        "--config-file",
        action="store_true",
        dest="configFile",
        help="Write config file"
    )
    parser.add_argument(
        "--fc",
        type=str,
        default="1.0",
        help="Multiplicative factor for comparison file"
    )
    parser.add_argument(
        "--fo",
        type=str,
        default="1.0",
        help="Multiplicative factor for original file"
    )
    parser.add_argument(
        "-i",
        action="store_true",
        help="Conversion from IVAR to error"
    )
    parser.add_argument(
        "-k",
        action="store_false",
        help="Use X,Y instead of sky coords for computing fiber distance"
    )
    parser.add_argument(
        "-m",
        action="store_false",
        help="Do NOT use masked arrays for flagged values"
    )
    parser.add_argument(
        "-e",
        action="store_true",
        help="Position table is in RSS and -p indicates the extension"
    )
    parser.add_argument("-w", help="HDU number extension for the wavelength array")
    parser.add_argument("-c", type=str, help="FITS file for comparison")

    return parser.parse_args()


def configure_backend(backend: str, available_backends: list) -> None:
    """Configura el backend de Matplotlib."""
    if backend not in available_backends:
        print(f'*** Backend "{backend}" NOT available ***')
        print(f">>> Available backends: {' | '.join(available_backends)}")
        sys.exit(1)
    matplotlib.use(backend)


def handle_config_and_version(args: argparse.Namespace) -> None:
    """Maneja las opciones de configuración y versión."""
    if args.configFile:
        cfgfile = vc.viewcuberc if not os.path.exists(vc.configfile) else vc.configfile
        vc.WriteConfigFile(filerc=cfgfile)
        sys.exit()

    if args.v:
        print(f'ViewCube version: {version.__version__}')
        sys.exit()


def update_config_from_args(args: argparse.Namespace) -> None:
    """Actualiza el diccionario de configuración global desde los argumentos."""
    config_mappings = {
        "specaxis": args.s,
        "exdata": args.data,
        "exerror": args.error,
        "exflag": args.flag,
        "exhdr": args.header,
        "masked": args.m,
        "skycoord": args.k,
        "sensf": args.f
    }

    for key, value in config_mappings.items():
        if value is not None:
            vc.defaultDictParams[key] = value

    if args.w is not None:
        try:
            vc.defaultDictParams["exwave"] = int(args.w)
        except ValueError as e:
            print(f"Error converting wavelength extension: {e}")
            vc.defaultDictParams["exwave"] = args.w


def setup_viewer(args: argparse.Namespace, config_params: dict):
    """Configura y devuelve el visor apropiado según los argumentos."""
    fc = float(args.fc) if args.fc else 1.0
    fo = float(args.fo) if args.fo else 1.0

    if args.p is None:
        config_params.pop("angle", None)
        config_params.pop("skycoord", None)
        return CubeViewer(
            args.name[0],
            fitscom=args.c,
            fc=fc,
            fo=fo,
            ivar=args.i,
            **config_params
        )

    excluded_keys = [
        "exdata", "exwave", "exhdr", "specaxis", "exerror", "exflag",
        "c", "mval", "cflag", "lflag", "lcom", "ccom", "lspec", "cspec",
        "dsoni", "ref_mode", "soni_start"
    ]

    for key in excluded_keys:
        config_params.pop(key, None)

    if args.a is not None:
        config_params["angle"] = args.a

    return RSSViewer(
        args.name[0],
        args.p,
        fitscom=args.c,
        fc=fc,
        fo=fo,
        ivar=args.i,
        extension=args.e,
        **config_params
    )


def main() -> None:
    """Función principal de ejecución de ViewCube."""
    args = parse_arguments()

    if not args.name:
        sys.exit("Error: Se requiere al menos un archivo FITS como argumento")

    handle_config_and_version(args)
    configure_backend(args.b, matplotlib.rcsetup.interactive_bk)

    if args.y:
        plt.style.use(args.y.split(','))

    vc.GetConfig(vc.defaultDictParams)
    update_config_from_args(args)
    viewer = setup_viewer(args, vc.defaultDictParams.copy())
    viewer.run()


if __name__ == "__main__":
    main()
