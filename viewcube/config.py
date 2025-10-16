############################################################################
#                              VIEWCUBE CONFIG                             #
#                                 PYTHON 3                                 #
#                                                                          #
# RGB@IAA ---> Last Change: 2024/11/30                                     #
############################################################################
#
#
#
################################ VERSION ###################################
VERSION = "0.1.0"                                                          #
############################################################################
#

from .utils import LoadFits, check_list
from collections import OrderedDict
from .version import __version__
import os, sys, ast
import inspect

if sys.version_info >= (3, 0):
    import configparser as ConfigParser
else:
    import ConfigParser

viewcuberc = "viewcuberc"
configfile = os.path.join(os.environ.get('HOME', ''), '.' + viewcuberc)

# --------------------------------------------------------------------------
defaultParams = [
    ["exdata", None, "Extension Data"],
    ["exhdr", None, "Extension Header"],
    ["exwave", None, "External Wavelength"],
    ["exflag", None, "FLAG/MASK Extension"],
    ["exerror", None, "ERROR Extension"],
    ["specaxis", None, "Spectral Dimension"],
    ["dfilter", None, "Filters directory"],
    ["dsoni", None, "Sonification directory"],
    ["ref_mode", "crpix", "Reference (central) pixel [crpix | max]"],
    ["soni_start", False, "Activate sonification mode (import libraries and check database)"],
    ["norm", "sqrt", "Normalization for colorbar"],
    ["default_filter", None, "Default filter"],
    # ['fcom',               1.0,            'Multiplicative factor'],
    ["mval", 0.0, "Masked values"],
    ["palpha", 0.95, "Alpha channel value"],
    ["plw", 0.1, "Line width"],
    ["plc", "k", "Color"],
    ["clw", 1, "Line width selection"],
    ["cc", "r", "Color selection"],
    ["cf", False, "True or False"],
    ["ca", 0.8, "Alpha selection"],
    ["slw", 2, "Width"],
    ["sf", False, "Save Fits"],
    ["sa", 0.9, "Alpha"],
    ["cspec", "#1f77b4", "Spectra color"],
    ["lspec", 1, "Spectra linewidth"],
    ["ccom", "#ff7f0e", "Comparison spectra color"],
    ["lcom", 1, "Comparision spectra linewidth"],
    ["cflag", "r", "Flag color"],
    ["lflag", 1, "Flag linewidth"],
    ["colorbar", True, "Colorbar"],
    ["fits", False, "Save fits"],
    ["txt", True, "Save ASCII"],
    ["integrated", True, "Save integrated"],
    ["individual", False, "Save individual spectra"],
    ["wlim", None, "Wavelength limits"],
    ["flim", None, "Flux limits"],
    ["iclm", True, "Colorbar limits"],
    ["fp", 1.2, "Factor"],
    ["fig_spaxel_size", (7.1, 6), "Size spaxel viewer"],
    ["fig_spectra_size", (8, 5), "Size spectra viewer"],
    ["fig_window_manager", (5, 5), "Size window manager"],
    ["c", 299792.458, "Speed of light (km/s)"],
    ["cfilter", False, "Center filter in wavelength range"],
    ["remove_cont", False, "Remove continuum"],
    ["angle", None, "Angle for rotating the position table (only RSS)"],
    ["skycoord", True, "Use skycoord to compute fibre distance; x,y otherwise"],
    ["masked", True, "Use masked arrays for flux (flag = mask)"],
    ["vflag", 0, 'Flags with values larger than "vflag" are considered flagged'],
]
# --------------------------------------------------------------------------


# ---------------------------------------------------------------------------
def WriteConfigFile(idict=defaultParams, cdict=None, filerc=viewcuberc):
    ff = open(filerc, "w")
    ff.write("[VIEWCUBE] #%% Version: %s\n" % __version__)
    if isinstance(idict, (tuple, list)):
        for key, val, com in idict:
            if key != "version":
                ff.write("#%% %s\n" % (com))
                ff.write("#%-18s   : %s\n\n" % (key, val))
    else:
        for key in list(idict.keys()):
            if key != "version":
                if cdict != None:
                    ff.write("#%% %s\n" % (cdict[key]))
                if defaultParams[key] != idict[key]:
                    ff.write("%-18s   : %s\n\n" % (key, idict[key]))
                else:
                    ff.write("#%-18s   : %s\n\n" % (key, idict[key]))
    ff.close()

    # Make sure the file is writable
    #os.chmod(filerc, 0o666)

    print('*** ViewCube Config File "%s" created! ***' % (filerc))


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
def get_func_default_parameters(func, remove=None):
    signature = inspect.signature(func)
    ndict = {}
    for par in signature.parameters:
        # Get only default parameters
        if not isinstance(signature.parameters[par].default, inspect._empty):
            ndict[signature.parameters[par].name] = signature.parameters[par].default
    if remove is not None:
        for key in check_list(remove):
            ndict.pop(key, None)
    return ndict


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
def GetDefaultParams(idict=defaultParams, value=True):
    id = 0 if value else 1
    ndict = OrderedDict()
    if isinstance(idict, dict):
        for key in list(idict.keys()):
            ndict[key] = idict[key][id]
    else:
        for item in idict:
            ndict[item[0]] = item[id + 1]
    return ndict


# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
def GetConfig(dictParams, filerc=viewcuberc):
    rcname = os.path.join(os.curdir, filerc)
    if not os.path.exists(rcname):
        rcname = os.path.join(os.environ["HOME"], "." + filerc)
        if not os.path.exists(rcname):
            print('*** ViewCube configuration file "%s" NOT found!! ***' % filerc)
            print('*** Run "ViewCube --config-file" to create one! ***')
            sys.exit()
    confp = ConfigParser.ConfigParser()
    confp.read_file((open(rcname)))
    # if len(confp.items('QBICK')) == 0:
    # print '*** QBICK configuration file "%s" EMPTY!! Set default variables!! ***' % filerc
    def_pars = get_func_default_parameters(LoadFits, remove="name")
    for key, val in confp.items("VIEWCUBE"):
        if key in list(dictParams.keys()) or key in def_pars:
            dictParams[key] = ast.literal_eval(val)
        else:
            sys.exit('*** Key "%s" not found!!! ***' % key)


# ---------------------------------------------------------------------------


defaultDictParams = GetDefaultParams(defaultParams)
