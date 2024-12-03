############################################################################
#                            VIEW-CUBE/RSS                                 #
#                               PYTHON 3                                   #
#                                                                          #
# RGB@IAA ---> Last Change: 2024/11/30                                     #
############################################################################
#
#
#
################################ VERSION ###################################
VERSION = "0.0.7"                                                          #
############################################################################
#
"""
 Author: Ruben Garcia-Benito (RGB)
 """
import viewcube.version as version
import viewcube.config as vc
import matplotlib.rcsetup
import matplotlib
import argparse
import sys
import os

def main():

    # List of matplotlib backends
    list_backends = matplotlib.rcsetup.interactive_bk
    slist_backends = " | ".join(list_backends)

    # Parse options
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )  # --> Escribir defaults en help
    parser.add_argument("--data", type=int, help="DATA extension")
    parser.add_argument("--error", type=int, help="ERROR extension")
    parser.add_argument("--flag", type=int, help="FLAG/MASK extension")
    parser.add_argument("--header", type=int, default=0, help="HEADER extension")
    parser.add_argument("-a", type=float, help="Angle to rotate the position table (only RSS)")
    parser.add_argument(
        "-b",
        type=str,
        help="Matplotlib backend. Use 'TkAgg' if using PyRAF for interactive fitting. Available backends: %s"
        % slist_backends,
        default=matplotlib.rcParams["backend"],
    )
    parser.add_argument("-c", type=str, help="FITS file for comparison")
    parser.add_argument(
        "-e",
        action="store_true",
        help="Position table is in RSS and -p indicates the extension (string or int)",
    )
    parser.add_argument(
        "-f", action="store_false", help="Do NOT apply sensitivity curve (if HDU is available)"
    )
    parser.add_argument("-fo", type=str, default="1.0", help="Multiplicative factor for original file")
    parser.add_argument(
        "-fc", type=str, default="1.0", help="Multiplicative factor for comparison file"
    )
    parser.add_argument("-i", action="store_true", help="Conversion from IVAR to error")
    parser.add_argument(
        "-k", action="store_false", help="Use X,Y instead of sky coords for computing fiber distance"
    )
    parser.add_argument("-m", action="store_false", help="Do NOT use masked arrays for flagged values")
    parser.add_argument("-p", type=str, help="External position table for RSS Viewer")
    parser.add_argument("-s", type=int, help="Spectral dimension")
    parser.add_argument(
        "-y", type=str, help="Plot style, separated by comma: 'dark_background, seaborn-ticks'"
    )
    parser.add_argument("-v", help="Print version", action="store_true")
    parser.add_argument("-w", help="HDU number extension for the wavelength array")
    parser.add_argument(
        "--config-file", action="store_true", dest="configFile", help="Write config file"
    )
    parser.add_argument("name", type=str, help="FITS file", nargs="*")
    args = parser.parse_args()
    
    # Config File
    if args.configFile:
        cfgfile = vc.viewcuberc
        if not os.path.exists(vc.configfile):
            cfgfile = vc.configfile
        vc.WriteConfigFile(filerc=cfgfile)
        sys.exit()

    # Print version
    if args.v:
        print ('ViewCube version: %s' % version.__version__)
        sys.exit()
    
    # Exception arguments
    if args.name is None or len(args.name) < 1:
        sys.exit(parser.print_usage())
    
    # Matplotlib Backend
    if args.b not in list_backends:
        print('*** Backend "%s" NOT available ***' % args.b)
        print(">>> Available backends: %s" % slist_backends)
        sys.exit()
    else:
        matplotlib.use(args.b)

    # Multiplicative factor
    fc = float(eval(args.fc))
    fo = float(eval(args.fo))
    
    # Read Config File
    vc.GetConfig(vc.defaultDictParams)
    if args.w != None:
        try:
            vc.defaultDictParams["exwave"] = int(args.w)
        except:
            vc.defaultDictParams["exwave"] = args.w
    
    # LoadFits options
    vc.defaultDictParams["specaxis"] = args.s
    vc.defaultDictParams["exdata"] = args.data
    vc.defaultDictParams["exerror"] = args.error
    vc.defaultDictParams["exflag"] = args.flag
    vc.defaultDictParams["exhdr"] = args.header
    vc.defaultDictParams["masked"] = args.m
    vc.defaultDictParams["skycoord"] = args.k
    vc.defaultDictParams["sensf"] = args.f
    
    if args.y is not None:
        import matplotlib.pyplot as plt
    
        styles = args.y.split(",")
        plt.style.use(styles)
    
    if args.p is None:
        from viewcube.cubeviewer import CubeViewer
    
        lkey = ["angle", "skycoord"]
        for key in lkey:
            if key in vc.defaultDictParams:
                del vc.defaultDictParams[key]
        CV = CubeViewer(args.name[0], fitscom=args.c, fc=fc, fo=fo, ivar=args.i, **vc.defaultDictParams)
    
    else:
        from viewcube.viewrss import RSSViewer
    
        lkey = [
            "exdata",
            "exwave",
            "exhdr",
            "specaxis",
            "exerror",
            "exflag",
            "c",
            "mval",
            "cflag",
            "lflag",
            "lcom",
            "ccom",
            "lspec",
            "cspec",
            "dsoni",
            "ref_mode",
            "soni_start",
        ]
        for key in lkey:
            if key in vc.defaultDictParams:
                del vc.defaultDictParams[key]
        if args.a is not None:
            vc.defaultDictParams["angle"] = args.a
        CV = RSSViewer(
            args.name[0],
            args.p,
            fitscom=args.c,
            fc=fc,
            fo=fo,
            ivar=args.i,
            extension=args.e,
            **vc.defaultDictParams
        )
