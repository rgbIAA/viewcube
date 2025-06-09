############################################################################
#                               VIEWRSS                                    #
#                               PYTHON 3                                   #
#                                                                          #
# RGB@IAA ---> Last Change: 2024/10/10                                     #
############################################################################
#
#
#
################################ VERSION ###################################
VERSION = "0.1.5"                                                          #
############################################################################
#
"""
 Author: Ruben Garcia-Benito (RGB)
 """
from .utils import ofits, LoadFits, lsfiles, ckfiles, save_spec, get_min_max, convert2iraf_spec
from .cubeviewer_old import GetSpaxelLimits, GetLambdaLimits, GetIdFilter, GetFluxLimits
from matplotlib.collections import PatchCollection, PolyCollection
from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.widgets import RectangleSelector
from astropy.coordinates import SkyCoord
import matplotlib.transforms as mtrans
from .rgbmpl import rnorm, IntColorMap
import matplotlib, fnmatch, sys, os
from scipy.spatial import distance
import matplotlib.cbook as cbook
import astropy.io.fits as pyfits
from matplotlib import rcParams
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy import units as u
import numpy as np
import argparse
import math

# Check module "pyraf"
try:
    from pyraf.iraf import splot

    PYRAF = True
    # Change before matplotlib.pyplot is called by any program
    # matplotlib.use('TkAgg')
except:
    PYRAF = False

# Check module "pyspeckit" (after pyraf)
try:
    import pyspeckit

    PYSPEC = True
except:
    PYSPEC = False


def hexagon(x, y, scale=1.0, m_sqrt3=None, extent=None, patch=True, angle=None, **kwargs):
    if m_sqrt3 is None:
        m_sqrt3 = math.sqrt(3.0)

    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

    x, y = cbook.delete_masked_points(x, y)

    x = np.array(x, float)
    y = np.array(y, float)

    sx = 2 / m_sqrt3 * scale * 0.99
    sy = scale * 0.99

    if extent is not None:
        xmin, xmax, ymin, ymax = extent
    else:
        xmin, xmax = (np.amin(x - sx), np.amax(x + sx)) if len(x) else (0, 1)
        ymin, ymax = (np.amin(y - sy), np.amax(y + sy)) if len(y) else (0, 1)

        # to avoid issues with singular data, expand the min/max pairs
        xmin, xmax = mtrans.nonsingular(xmin, xmax, expander=0.1)
        ymin, ymax = mtrans.nonsingular(ymin, ymax, expander=0.1)

    padding = 1.0e-9 * (xmax - xmin)
    xmin -= padding
    xmax += padding

    n = len(x)
    polygon = np.zeros((6, 2), float)

    S = 1.0 / m_sqrt3
    mx = my = 0.99 * scale
    polygon[:, 0] = mx * np.array([-0.5 * S, 0.5 * S, 1.0 * S, 0.5 * S, -0.5 * S, -1.0 * S])
    polygon[:, 1] = my * np.array([0.5, 0.5, 0.0, -0.5, -0.5, 0.0])

    if angle is not None:
        polygon[:, 0], polygon[:, 1] = rotatePosTable(polygon[:, 0], polygon[:, 1], angle=angle)

    if patch:
        polygons = [Polygon(polygon + np.array([ix, iy]), **kwargs) for (ix, iy) in zip(x, y)]

    else:
        offsets = np.zeros((n, 2), float)
        offsets[:, 0] = x
        offsets[:, 1] = y

        polygons = PolyCollection(
            [polygon],
            offsets=offsets,
            offset_position="data",
            transOffset=mtrans.IdentityTransform(),
            **kwargs
        )

    return polygons


def get_radius(x, y):
    points = list(zip(x, y))
    dist = np.unique(distance.cdist(points, points, "euclidean").flatten())
    radius = dist[dist > 0.0][0]
    return radius


def readMegaraPosTable(
    fits, extname="FIBERS", ext=None, table=True, lkeys=None, angle=None, ref_x=0.0, ref_y=0.0
):
    from collections import OrderedDict

    if lkeys is None:
        lkeys = ["X", "Y", "D", "R", "A", "W1", "W2", "T"]

    if ext is not None:
        extname = None
    hdr = pyfits.getheader(fits, extname=extname, ext=ext)
    df = OrderedDict()

    for key in hdr:
        if key.startswith("FIB") and "_" in key:
            fid = int(key.split("_")[0].replace("FIB", ""))
            if not fid in df:
                df[fid] = {}
            ikey = key.split("_")[-1]
            if ikey in lkeys:
                df[fid][ikey] = hdr[key]

    df = OrderedDict(sorted(df.items()))

    if table:
        from astropy.table import Table

        t = Table(list(df.values()))

        t["ID"] = list(df.keys())
        lkeys = ["ID"] + [key for key in lkeys if key in t.colnames]
        df = t[lkeys]

    if angle is not None:
        df["X"], df["Y"] = rotatePosTable(df["X"], df["Y"], angle, ref_x=ref_x, ref_y=ref_y)

    return df


def readLIFUPosTable(fits, extname="FIBTABLE", skycoord=True, angle=None, ref_x=0.0, ref_y=0.0):
    t = Table.read(fits, hdu=extname)
    if skycoord:
        try:
            co = SkyCoord(t["FIBRERA"] * u.deg, t["FIBREDEC"] * u.deg)
        except:
            # Some tables come with the units already set
            co = SkyCoord(t["FIBRERA"], t["FIBREDEC"])
        idc = np.argmin(t["XPOSITION"] ** 2 + t["YPOSITION"] ** 2)
        ra, dec = co.ra, co.dec
        ra, dec = co.spherical_offsets_to(co[idc])
        t["X"] = ra.to(u.arcsec).value
        t["Y"] = dec.to(u.arcsec).value
    else:
        t["X"] = t["XPOSITION"]
        t["Y"] = t["YPOSITION"]

    if angle is not None:
        t["X"], t["Y"] = rotatePosTable(t["X"], t["Y"], angle, ref_x=ref_x, ref_y=ref_y)

    return t


def rotatePosTable(x, y, angle, ref_x=0.0, ref_y=0.0):
    nx = (x - ref_x) * np.cos(angle / 180.0 * np.pi) - (y - ref_y) * np.sin(angle / 180.0 * np.pi)
    ny = (x - ref_x) * np.sin(angle / 180.0 * np.pi) + (y - ref_y) * np.cos(angle / 180.0 * np.pi)
    return nx, ny


class RSSViewer:
    def __init__(
        self,
        name_fits,
        ptable,
        fitscom=None,
        dfilter="filters/",
        norm="sqrt",
        default_filter="Halpha_KPNO-NOAO",
        fo=1.0,
        fc=1.0,
        extension=False,
        palpha=0.95,
        plw=0.1,
        plc="k",
        clw=1,
        cc="r",
        cf=False,
        ca=0.8,
        slw=2,
        sf=False,
        sa=0.9,
        colorbar=False,
        fits=False,
        txt=True,
        integrated=True,
        mval=0.0,
        individual=False,
        wlim=None,
        flim=None,
        iclm=True,
        fig_spaxel_size=(7, 6),
        fig_spectra_size=(8, 5),
        fig_window_manager=(5, 5),
        fp=1.2,
        ft="C",
        hex_scale=None,
        extent=None,
        cfilter=False,
        remove_cont=False,
        angle=None,
        skycoord=True,
        masked=True,
        vflag=0,
        c=299792.458,
        **kwargs
    ):
        """
        # -----------------------------------------------------------------------------
              USO:
                      rssv = RSSViewer('ngc.fits','ngc.pt.txt')
        # -----------------------------------------------------------------------------
              name_fits --> Name of the fits file
              ptbale --> Name of the position table
              fitscom = None --> Name of the RSS fits file you want to compare
              dfilter = 'filters/' --> Directory where the filters are located
              norm = 'sqrt' --> Normalize function for the color map
              default_filter = 'Halpha_KPNO-NOAO' --> Root name of the default filter
              fo = 1.0 --> Multiplicative factor for the original file
              fc = 1.0 --> Multiplicative factor for "fitscom" file
              extension = False --> Position Table is given as an extension HDU of
                      "name_fits" and "ptable" indicates the extension (string or integer)
              palpha = 0.95 --> Alpha value of the RSS spaxels
              plw = 0.1 --> Linewidth of the RSS spaxels
              plc = 'k' --> Color of the border of the RSS spaxels
              clw = 1 --> Linewidth of the RSS spaxel selector
              cf = False --> Fill value of the RSS selected spaxels
              ca = 0.8 --> Alpha value of the RSS selected spaxels
              slw = 2 --> Linewidth of the spectra-->spaxel-spaxel
              sf = False --> Fill value of the spectra-->spaxel-spaxel
              sa = 0.9 --> Alpha value of the spectra-->spaxel-spaxel
              colorbar = False --> Set colorbar in the RSS spaxel viewer
              fits = False --> Save files in fits type
              txt = True --> Save files in ASCII type
              integrated = True --> Save the integrated spectrum of the selected spaxels
              individual = False --> Save the individual spectra of the selected spaxels
              wlim = None --> 2D-Tuple with the wavelength limits of the spectra inspector
              flim = None --> 2D-Tuple with the flux limits of the spectra inspector
              iclm = True --> DS9-like option for changing the dynamical range of the
                      RSS spaxel Viewer
              fig_spaxel_size = (7,6) --> 2D-Tuple with the size of the RSS spaxel viewer
              fig_spectra_size = (8,5) --> 2D-Tuple with the size of the Spectral inspector
              fig_window_manager = (5,5) --> 2D-Tuple with the size of the Window Manager
              fp = 1.2 --> The filter is multipy by the maximum value of the spectra and
                      by this multiplicative constant
              ft = 'C' --> Fiber type (C = circle; H = hexagon)
              hex_scale = None --> Scale of the hexagon if ft = 'H'
                      (ex: MEGARA = 0.443 mm from center to center, upwards)
              extent = None --> (xlim, xmax, ylim, ymax) in the spatial (spaxel) viewer
              cfilter = False --> Center filter in wavelength range
              remove_cont = False --> Remove continuum from adjacent positions (right and left)
              angle = None -> Rotation angle of the position table (only for RSS)
              skycoord = True --> Use skycoord to compute relative distance in fibers;
                  use X, Y position otherwise
              masked = True --> Use masked arrays for flux (flag = mask)
              vflag = 0 --> Flags with values larger than "vflag" are considered flagged
        # -----------------------------------------------------------------------------
        """

        # Disable default matplotlib shorcuts
        rcParams["keymap.save"] = ""
        rcParams["keymap.yscale"] = ""
        rcParams["keymap.xscale"] = ""

        # Version
        self.version = "RSSViewer version: 1.0.0"

        # Set variables
        self.name_fits = name_fits
        self.bname_fits = os.path.basename(name_fits)
        self.ptable = ptable
        self.extension = extension
        self.fitscom = fitscom
        self.errcom = False
        self.dfilter = dfilter
        self.norm = norm
        self.default_filter = default_filter
        self.list_filters = None
        self.palpha = palpha
        self.plw = plw
        self.plc = plc
        self.colorbar = colorbar
        self.extent = extent
        self.fc = fc
        self.fo = fo
        # Constant
        self.c = c
        # Position table
        self.radius = None
        self.xs = None
        self.ys = None
        self.x = None
        self.y = None
        # Select Patch Properties
        self.hex_scale = hex_scale
        self.ft = ft
        self.clw = clw
        self.cf = cf
        self.ca = ca
        self.cc = cc
        # Select Spaxel Properties
        self.slw = slw
        self.sf = sf
        self.sa = sa
        self.ids = None
        # Color Map dynamic range
        self.iclm = iclm
        # Spaxel viewer props
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        # Size Windows
        self.fig1_size = fig_spaxel_size
        self.fig2_size = fig_spectra_size
        self.winman_size = fig_window_manager
        # Labels
        self.fig1_label = "Spaxel Viewer"
        self.fig2_label = "Spectral Viewer"
        self.winman_label = "Window Manager"
        self.sint = "Integrated"
        # Flux
        self.flim = flim
        # Lambda
        self.wlim = wlim
        self.wrest = False
        # Save Spectra
        self.fits = fits
        self.txt = txt
        # Save type
        self.integrated = integrated
        self.individual = individual
        # Color Buttons
        self.axcolor1 = "lightgoldenrodyellow"  # "Off mode"
        self.axcolor2 = "#D0A9F5"  # Clicked color "On" mode
        self.hovercolor = "0.975"
        self.ccc = self.axcolor1 if not self.cf else self.axcolor2
        self.scc = self.axcolor1 if not self.sf else self.axcolor2
        self.tcc = self.axcolor1 if not self.integrated else self.axcolor2
        self.icc = self.axcolor1 if not self.individual else self.axcolor2
        self.xcc = self.axcolor1 if not self.txt else self.axcolor2
        self.fcc = self.axcolor1 if not self.fits else self.axcolor2
        self.clab = "Fill " + ("Off" if not self.cf else "On")
        self.slab = "Fill " + ("Off" if not self.sf else "On")
        self.intlab = "Integrated " + ("Off" if not self.integrated else "On")
        self.indlab = "Individual " + ("Off" if not self.individual else "On")
        self.txtlab = "txt " + ("Off" if not self.txt else "On")
        self.fitlab = "fits " + ("Off" if not self.fits else "On")
        # Pick of the Passband respect to the maximum value
        self.fp = fp
        # Modes
        self.mode = True  # Spectra mode
        self.list = []
        self.pat = []
        self.cir = []
        self.spec_mode = 2
        self.pressevent = None
        self.cfilter = cfilter
        self.remove_cont = remove_cont
        self.angle = angle
        # Fit Spec
        self.fitspec = True
        self.pyraf = PYRAF
        self.pyspec = PYSPEC
        self.fitspec_mode = 0  # 0 = PYRAF | 1 = PYSPECKIT
        if not self.pyraf and not self.pyspec:
            self.fitspec = False

        # Basic Info
        self.fobj = LoadFits(self.name_fits, flip2disp=False, **kwargs)
        self.dat = self.fobj.data * self.fo
        self.syn = self.fobj.syn
        self.err = self.fobj.error
        if self.err is not None:
            self.err *= self.fo
        self.flag = self.fobj.flag
        self.perr = None
        if self.flag is not None:
            if self.flag.dtype.kind != "b":
                self.flag = np.where(self.flag > 0, True, False)
            self.perr = self.dat.copy().astype(float)
            self.perr[~self.flag] = np.nan
            if self.err is not None:
                self.err[self.flag] = np.nan
        if mval is not None:
            self.dat[~np.isfinite(self.dat)] = mval
        if masked and self.flag is not None:
            self.dat = np.ma.array(
                self.dat, mask=np.logical_or(self.flag > 0, ~np.isfinite(self.dat))
            )

        self.hd = self.fobj.hdr
        self.crval = self.fobj.crval
        self.cdelt = self.fobj.cdelt
        self.redshift = self.fobj.redshift
        self.velocity = self.fobj.velocity
        shp = np.shape(self.dat)
        naxes = len(shp)
        self.root = ".".join(self.name_fits.split(".")[:-1])
        # min, max = self.dat.min(), self.dat.max()
        self.wl2 = None
        if self.fitscom is not None:
            # self.dat2 = ofits(self.fitscom,hd=False)
            self.fobj2 = LoadFits(self.fitscom, flip2disp=False, **kwargs)
            self.wl2 = self.fobj2.wave
            self.dat2 = self.fc * self.fobj2.data
            self.err2 = self.fc * self.fobj2.error if self.fobj2.error is not None else None

        # Lambda
        if naxes == 2:
            self.ns = shp[0]  # Number of spectra
            self.nl = shp[1]  # Number of lambdas
            # crval = self.hd['CRVAL1']
            # cdelt = self.hd['CDELT1']
            # crpix = self.hd['CRPIX1']

        self.orig_wl = self.fobj.wave
        self.orig_wl_rest = self.fobj.wave_rest
        self.wl = self.fobj.wave
        self.wl_rest = self.fobj.wave_rest
        if self.syn is not None and self.velocity is not None:
            self.orig_wl_rest = self.fobj.wave_rest
            self.orig_wl = self.vredshift(-self.velocity)
        if self.wl is None:
            print("*** WAVELENGTH DATA NOT FOUND!!! ***")
            self.wl = np.arange(self.nl)

        self.wlmin = self.wl.min()
        self.wlmax = self.wl.max()
        self.pbline = None
        self.pband = None

        # Filtros
        if self.dfilter is not None:
            self.list_filters = lsfiles("*txt", self.dfilter, path=False)
        self.ifil = None
        if self.list_filters is not None:
            self.ifil = GetIdFilter(self.list_filters, self.default_filter, self.dfilter)
            self.nfil = len(self.list_filters)
        self.set_fff(verb=False)
        self.gdl = 0  # Global delta position filter
        self.dl = None
        # color = dat[:,0:1900].sum(axis=1)

        # Read Position Table
        self.readPositionTable(extension, skycoord=skycoord)

        # Create Figure 2 (Spectral Viewer) ---------------------------------
        self.fig2 = plt.figure(2, self.fig2_size)
        self.fig2.set_label(self.fig2_label)
        self.setWindowTitle(self.fig2, self.fig2_label)
        self.ax2 = self.fig2.add_subplot(111)
        self.awlmin, self.awlmax = GetLambdaLimits((self.wl, self.wl2), 0.05, wlim=self.wlim)
        self.fmin, self.fmax = GetFluxLimits(self.flim)
        self.ax2.set_xlim((self.awlmin, self.awlmax))
        self.ax2.set_ylim((self.fmin, self.fmax))
        # self.ax2.xaxis.get_major_formatter()
        # self.fig2.canvas.set_window_title("Spectral Viewer by RGB")

        # Create Figure 1 (Spaxel Viewer) ---------------------------------
        self.fig = plt.figure(1, self.fig1_size)
        self.fig.set_label(self.fig1_label)
        self.setWindowTitle(self.fig, self.fig1_label)
        self.ax = self.fig.add_subplot(111)
        if self.xmin is None or self.xmax is None or self.ymin is None or self.ymax is None:
            self.xmin, self.xmax, self.ymin, self.ymax = GetSpaxelLimits(
                self.x, self.y, self.radius
            )
        self.ax.axis([self.xmin, self.xmax, self.ymin, self.ymax])
        self.ax.set_title(self.bname_fits)
        self.p = self.getPatches(
            cmap=plt.cm.gray_r,
            norm=rnorm(self.norm),
            alpha=self.palpha,
            lw=self.plw,
            color=self.plc,
        )
        self.p.set_array(self.color)
        self.ax.add_collection(self.p)
        if self.colorbar is True:
            self.fig.colorbar(self.p)
        self.cmin, self.cmax = self.color.min(), self.color.max()
        self.p.set_clim([self.cmin, self.cmax])

        # Bind Events
        self.fig.canvas.mpl_connect("key_press_event", self.PressKey)
        self.fig2.canvas.mpl_connect("key_press_event", self.PressKey)
        self.fig.canvas.mpl_connect("key_press_event", self.ChangeFilter)
        self.fig2.canvas.mpl_connect("key_press_event", self.ChangeFilter)
        self.fig.canvas.mpl_connect("motion_notify_event", self.SpectraViewer)
        self.fig.canvas.mpl_connect("button_press_event", self.SpectraViewer)
        self.fig2.canvas.mpl_connect("button_press_event", self.SpectraViewer)
        self.fig2.canvas.mpl_connect("pick_event", self.GetSpectraInfo)
        self.fig2.canvas.mpl_connect("button_press_event", self.PassBandPress)
        self.fig2.canvas.mpl_connect("motion_notify_event", self.PassBandMove)
        self.fig2.canvas.mpl_connect("button_release_event", self.PassBandRelease)
        self.fig2.canvas.mpl_connect("key_press_event", self.ErrorSpec)
        self.fig.canvas.mpl_connect("key_press_event", self.ErrorSpec)
        self.fig2.canvas.mpl_connect("key_press_event", self.Redshift)
        self.fig2.canvas.mpl_connect("key_press_event", self.RestWave)
        self.fig2.canvas.mpl_connect("key_press_event", self.FitSpec)

        # DS9-like dynamic range of colormap
        if self.iclm:
            self.icmp = IntColorMap(self.p)
        plt.show()

    def readCalifaPosTable(self):
        try:
            ptarray, pthd = pyfits.getdata(self.name_fits, extname=self.ptable, header=True)
        except:
            ptarray, pthd = pyfits.getdata(self.name_fits, ext=int(self.ptable), header=True)
        self.ft, self.xs, self.ys = pthd["FIBSHAPE"].strip(), pthd["FIBSIZEX"], pthd["FIBSIZEY"]
        self.x, self.y = ptarray[pthd["TTYPE1"]], ptarray[pthd["TTYPE2"]]

    def readMegaraPosTable(self, sky=False):
        t = readMegaraPosTable(self.name_fits, table=True)
        self.x = t["X"]
        self.y = t["Y"]
        self.ft = "H"
        if not sky:
            self.xmin = -6.5
            self.xmax = 6.5
            self.ymin = -6.0
            self.ymax = 6.0

    def readLIFUPosTable(self, sky=False, skycoord=True):
        t = readLIFUPosTable(self.name_fits, skycoord=skycoord)
        hdr = pyfits.getheader(self.name_fits, ext=0)
        mode = hdr.get("OBSMODE")
        if isinstance(mode, str):
            mode = mode.strip().upper()
            if "LIFU" in mode:
                self.xs = 2.6 / 2.0
            else:
                self.xs = 1.3 / 2.0
        self.x = t["X"]
        self.y = t["Y"]
        self.ft = "C"
        if not skycoord:
            # Hard coded, by eye, check!
            self.xs = 0.06
        if not sky and "LIFU" in mode:
            ids = t["TARGCLASS"] != "SKY"
            d = abs(self.x[ids].max() - self.x[ids].min()) / 10.0
            self.xmin = self.x[ids].min() - d
            self.xmax = self.x[ids].max() + d
            self.ymin = self.y[ids].min() - d
            self.ymax = self.y[ids].max() + d

    def readPositionTable(self, extension=False, skycoord=True):
        if not extension:
            ckfiles(self.ptable)
            with open(self.ptable, "r") as f:  # Close file
                self.ft, self.xs, self.ys, eid = f.readline().split()
            self.id, self.x, self.y = np.loadtxt(
                self.ptable, unpack=True, usecols=(0, 1, 2), skiprows=1
            )
        else:
            try:
                self.readCalifaPosTable()
            except:
                pass
            try:
                self.readMegaraPosTable()
            except:
                pass
            try:
                self.readLIFUPosTable(skycoord=skycoord)
            except:
                pass
        if self.angle is not None:
            self.x, self.y = rotatePosTable(self.x, self.y, angle=self.angle)

        if self.ft == "C":
            self.radius = float(self.xs)
            self.fiber_size = 2.0 * self.radius

        if self.ft == "H" and self.hex_scale is None:
            self.radius = get_radius(self.x, self.y)
            self.hex_scale = self.radius

        self.nsfm = len(str(self.ns))  # Number of digits of number of spaxels

    def PressKey(self, event):
        if event.key == "*":
            self.ax.cla()
            self.ax2.cla()
            self.ax.add_collection(self.p)
            self.ax.set_title(self.bname_fits)
            self.ax.axis([self.xmin, self.xmax, self.ymin, self.ymax])
            self.list = []
            self.fig.canvas.draw()
            self.fig2.canvas.draw()
        if event.key == "s":
            self.mode = not self.mode
        # Save selected Spectra
        if len(self.list) > 0 and event.key == "S":
            self.SaveFile()
        if event.key == "w":
            self.WindowManager()
        if event.key == "l":
            self.LambdaLimits()
        if event.key == "Y":
            self.FluxLimits()
        if event.key == "E":
            self.xyLimits()
        if event.key == "q":
            sys.exit()

    def setWindowTitle(self, fig, label):
        # self.backend = matplotlib.get_backend()
        try:
            fig.canvas.setWindowTitle(label)
        except:
            try:
                fig.canvas.set_window_title(label)
            except:
                pass

    def SaveFile(self):
        print('***** Actual Save File Options for "' + self.name_fits + '" *****')
        print("Input example:  ngc6211, fits=True, individual=True")
        print("%10s = %-5s | %10s = %-5s" % ("fits", str(self.fits), "txt", str(self.txt)))
        print(
            "%10s = %-5s | %10s = %-5s"
            % ("integrated", str(self.integrated), "individual", str(self.individual))
        )
        iname = input("Enter Root Name (Enter to abort): ")
        if len(iname) == 0:
            print("Nothing to save")
        else:
            # fname = self.root if len(name) == 0 else name
            if iname.find(",") != -1:
                fname = iname.split(",")[0]
                for item in iname.split(",")[1:]:
                    if item.lower().find("fits") != -1:
                        self.fits = eval(item.split("=")[1])
                    if item.lower().find("txt") != -1:
                        self.txt = eval(item.split("=")[1])
                    if item.lower().find("integrated") != -1:
                        self.integrated = eval(item.split("=")[1])
                    if item.lower().find("individual") != -1:
                        self.individual = eval(item.split("=")[1])
            else:
                fname = iname
            if (self.txt is not True and self.fits is not True) or (
                self.integrated is not True and self.individual is not True
            ):
                print("Nothing to save")
            else:
                if self.root == fname:
                    stint = "_int"
                else:
                    stint = ""
                if self.integrated:
                    infotxt = [
                        'Integrated Spectra extracted from "' + self.name_fits + '"',
                        "Sum of Spaxels (ID): " + " | ".join([str(id) for id in self.list]),
                    ]
                    infohd = [["RSSVWR_1", infotxt[0]], ["RSSVWR_2", infotxt[1]]]
                    self.intspec = np.take(self.dat, self.list, axis=0).sum(0)
                    save_spec(
                        self.wl,
                        self.intspec,
                        fname + stint,
                        fits=self.fits,
                        txt=self.txt,
                        hd=self.hd,
                        infohd=infohd,
                        infotxt=infotxt,
                    )
                if self.individual:
                    for item in self.list:
                        infotxt = [
                            'Spectra extracted from "' + self.name_fits + '"',
                            "Spaxel (ID): " + str(item),
                        ]
                        infohd = [["RSSVWR_1", infotxt[0]], ["RSSVWR_2", infotxt[1]]]
                        save_spec(
                            self.wl,
                            self.dat[item],
                            fname + "_" + str(item).zfill(self.nsfm),
                            fits=self.fits,
                            hd=self.hd,
                            txt=self.txt,
                            infohd=infohd,
                            infotxt=infotxt,
                        )
                print("Files Saved")

    def vredshift(self, cz):
        return self.wl / (1.0 + (cz / self.c))

    def LambdaLimits(self):
        print("***** Lambda Limits *****")
        print("Input example:  4000, None")
        lm = input("Enter lambda limits (Enter to abort): ")
        if len(lm) == 0:
            self.awlmin, self.awlmax = GetLambdaLimits((self.wl, self.wl2), 0.05, wlim=self.wlim)
        else:
            self.awlmin, self.awlmax = GetLambdaLimits((self.wl, self.wl2), 0.05, wlim=eval(lm))
        print("*** Lambda limits set to: (%6.1f, %6.1f) ***" % (self.awlmin, self.awlmax))
        self.ax2.set_xlim((self.awlmin, self.awlmax))
        self.fig2.canvas.draw()

    def FluxLimits(self):
        print("***** Flux Limits *****")
        print("Input example:  1e-16, None")
        lm = input("Enter flux limits (Enter to abort): ")
        if len(lm) == 0:
            pass
        else:
            self.fmin, self.fmax = GetFluxLimits(eval(lm))
            print("*** Flux limits set to: (%s, %s) ***" % (self.fmin, self.fmax))
            self.ax2.set_ylim((self.fmin, self.fmax))
            self.fig2.canvas.draw()

    def xyLimits(self):
        print("***** XY Limits [spaxel viewer (extent)] *****")
        print("Input example:  -6.5, None, -6, 6  [xmin, xmax, ymin, ymax]")
        lm = input("Enter extent limits (Enter to abort): ")
        if len(lm) == 0:
            pass
        else:
            lm = eval(lm)
            if not isinstance(lm, (tuple, list)) or len(lm) != 4:
                print(
                    "Flux limits should be a tuple or list of 4 items: ex. --> -6.5, None, -6, 6  [xmin, xmax, ymin, ymax]"
                )
                return
            xmin, xmax, ymin, ymax = lm
            if xmin is not None:
                self.xmin = xmin
            if xmax is not None:
                self.xmax = xmax
            if ymin is not None:
                self.ymin = ymin
            if ymax is not None:
                self.ymax = ymax
            print(
                "*** XY limits set to: (%s, %s, %s, %s) ***"
                % (self.xmin, self.xmax, self.ymin, self.ymax)
            )
            self.ax.set_xlim((self.xmin, self.xmax))
            self.ax.set_ylim((self.ymin, self.ymax))
            self.fig.canvas.draw()

    def Redshift(self, event):
        if event.key == "k":
            cz = input("Enter Redshift in km/s (Enter to abort): ")
            if len(cz) == 0:
                print("No Redshift provided")
            else:
                velocity = float(cz)
                self.wl = self.vredshift(float(cz))
                if self.velocity is None:
                    self.velocity = velocity
                if self.redshift is None:
                    self.redshift = velocity / self.c
                if self.orig_wl_rest is None:
                    self.orig_wl_rest = self.wl
                self.set_fff(verb=False, dl=self.gdl)
                self.updatePatch()
                self.PlotSpec()
                self.fig.canvas.draw()
                self.fig2.canvas.draw()

    def RestWave(self, event):
        if event.key == "z" and self.orig_wl_rest is not None:
            if self.wrest and self.velocity is not None:
                print(
                    "Rest Wavelength (redshift = %8.5f | velocity = %7.1f km/s)"
                    % (self.redshift, self.velocity)
                )
                self.wl = self.orig_wl_rest
            if not self.wrest and self.velocity is not None:
                print(
                    "Observed Wavelength (redshift = %8.5f | velocity = %7.1f km/s)"
                    % (self.redshift, self.velocity)
                )
                self.wl = self.orig_wl
            self.set_fff(verb=False, dl=self.gdl)
            self.updatePatch()
            self.PlotSpec()
            self.fig.canvas.draw()
            self.fig2.canvas.draw()
            self.wrest = ~self.wrest

    def updatePatch(self):
        color = self.color
        if self.fitscom is not None:
            if (self.wl[0] <= self.wmax) and (self.wmax <= self.wl[-1]):
                color = self.color
            if (self.wl2[0] <= self.wmax) and (self.wmax <= self.wl2[-1]):
                color = self.color2
        self.cmin, self.cmax = get_min_max(color)
        self.p.set_array(color)
        self.p.set_clim([self.cmin, self.cmax])

    def updatePassBand(self, remove=False):
        wl = self.wl if self.wl.ndim == 1 else self.wl[self.ids]
        fff = self.fff if self.fff.ndim == 1 else self.fff[self.ids]
        if self.pbline is not None and remove:
            self.pbline.pop(0).remove()
            self.pbline = None
            self.pband.remove()
            self.pband = None
        if self.fitscom is None:
            self.pbline = self.ax2.plot(wl, fff * self.dat[self.ids].max() * self.fp, "g")
            self.pband = self.ax2.fill_between(
                wl, fff * self.dat[self.ids].max() * self.fp, color="g", alpha=0.25
            )
        else:
            wl2 = self.wl2 if self.wl2.ndim == 1 else self.wl2[self.ids]
            if (wl[0] <= self.wmax) and (self.wmax <= wl[-1]):
                self.pbline = self.ax2.plot(wl, fff * self.dat[self.ids].max() * self.fp, "g")
                self.pband = self.ax2.fill_between(
                    wl, fff * self.dat[self.ids].max() * self.fp, color="g", alpha=0.25
                )
            elif (wl2[0] <= self.wmax) and (self.wmax <= wl2[-1]):
                self.pbline = self.ax2.plot(
                    wl2, self.fff2 * self.dat2[self.ids].max() * self.fp, "g"
                )
                self.pband = self.ax2.fill_between(
                    wl2, self.fff2 * self.dat2[self.ids].max() * self.fp, color="g", alpha=0.25
                )
            else:
                pass
        self.ax2.set_xlim((self.awlmin, self.awlmax))
        self.ax2.set_ylim((self.fmin, self.fmax))

    def ChangeFilter(self, event):
        if event.key in ["t", "T", "a", "c"]:
            if event.key == "t":
                self.ifil += 1
                if self.ifil >= self.nfil:
                    self.ifil = 0
            if event.key == "T":
                self.ifil -= 1
                if self.ifil < 0:
                    self.ifil = self.nfil - 1
            if event.key == "a":
                self.cfilter = not self.cfilter
            if event.key == "c":
                self.remove_cont = not self.remove_cont
            if event.key != "c":
                self.gdl = 0
            self.set_fff(verb=True, dl=self.gdl)
            self.updatePatch()
            self.updatePassBand(remove=True)
            self.fig.canvas.draw()
            self.fig2.canvas.draw()

    def getPatches(self, **kwargs):
        if self.ft == "C":
            patches = [Circle((self.x[i], self.y[i]), self.radius) for i in range(self.ns)]
        if self.ft == "H":
            patches = hexagon(self.x, self.y, scale=self.hex_scale, angle=self.angle)
        cpatch = PatchCollection(patches, **kwargs)
        return cpatch

    def getPatch(self, idf, **kwargs):
        if self.ft == "C":
            patch = Circle((self.x[idf], self.y[idf]), self.radius, **kwargs)
        if self.ft == "H":
            patch = hexagon(
                self.x[idf], self.y[idf], scale=self.hex_scale, angle=self.angle, **kwargs
            )[0]
        return patch

    def SpectraViewer(self, event):
        tb = plt.get_current_fig_manager().toolbar
        if event.inaxes and tb.mode == "":
            dist = np.hypot(self.x - event.xdata, self.y - event.ydata)
            idf = np.argmin(dist)
            if dist[idf] < self.radius:
                # Selection of Spectra --------------------------------
                if event.button == 1:
                    if len(self.cir) > 0:
                        [item.set_visible(False) for item in self.cir]  # item.remove() ?
                    if event.key == "d":
                        if idf in self.list:
                            idp = self.list.index(idf)
                            self.pat[idp].remove()  # self.pat[idp].set_visible(False) ?
                            self.list.remove(idf)
                            self.pat.pop(idp)
                    else:
                        if idf not in self.list:
                            self.pat.append(
                                self.getPatch(
                                    idf, color=self.cc, alpha=self.ca, lw=self.clw, fill=self.cf
                                )
                            )
                            self.list.append(idf)
                            self.ax.add_patch(self.pat[-1])
                        # xr,yr = PRectangle(x[idf],y[idf],radius); ax.fill(xr,yr,fill=False,ec='r')
                    self.fig.canvas.draw()
                # Mode Spectra Viewer ('s') ---------------------------------
                if self.mode:
                    self.ids = idf
                    self.PlotSpec()
                    self.fig2.canvas.draw()
            # View selected Spectra  -------------------------------
            # Right click
            if len(self.list) > 0 and event.button == 3:
                self.pband = None
                # Change mode of spectra view if we are in the SpectraViewer Figure
                if event.canvas.figure.get_label() == self.fig2_label:
                    self.spec_mode += 1
                self.ax2.cla()
                # Set selection spectra mode so other modes like error visualization do not activate
                self.ids = None
                # self.ax.add_collection(self.p)
                if (self.spec_mode % 3) == 0 or (self.spec_mode % 3) == 1:
                    for i in self.list:
                        wl = self.wl if self.wl.ndim == 1 else self.wl[i]
                        p = self.ax2.plot(wl, self.dat[i], label=str(i + 1), picker=True)
                        if self.fitscom is not None:
                            wl2 = self.wl2 if self.wl2.ndim == 1 else self.wl2[i]
                            self.ax2.plot(
                                wl2,
                                self.dat2[i],
                                label=str(i + 1),
                                picker=True,
                                c=p[0].get_color(),
                                alpha=0.7,
                            )
                if (self.spec_mode % 3) == 1 or (self.spec_mode % 3) == 2:
                    if self.wl.ndim == 1:
                        self.intspec = np.take(self.dat, self.list, axis=0).sum(0)
                        p = self.ax2.plot(self.wl, self.intspec, label=self.sint, picker=True)
                        if self.fitscom is not None:
                            intspec = np.ma.array([self.dat2[i] for i in self.list]).sum(0)
                            self.ax2.plot(
                                self.wl2,
                                intspec,
                                label=self.sint,
                                picker=True,
                                c=p[0].get_color(),
                                alpha=0.7,
                            )
                    else:
                        print(">>> Integrated spectrum NOT available! Wavelength array ndim = 2")
                self.fig2.canvas.draw()

    def PlotSpec(self):
        if self.ids is None:
            return
        self.ax2.cla()
        wl = self.wl if self.wl.ndim == 1 else self.wl[self.ids]
        self.ax2.plot(wl, self.dat[self.ids, :])
        if self.fitscom is not None:
            wl2 = self.wl2 if self.wl2.ndim == 1 else self.wl2[self.ids]
            if self.dat2.ndim == 2:
                self.ax2.plot(wl2, self.dat2[self.ids, :], label=self.fitscom)
            if self.dat2.ndim == 1:
                self.ax2.plot(wl2, self.dat2, label=self.fitscom)
        if self.errcom:
            self.ax2.errorbar(
                wl, self.dat[self.ids, :], yerr=self.err[self.ids, :], fmt="none", ecolor="grey"
            )
            if self.fitscom is not None and self.err2 is not None:
                self.ax2.errorbar(
                    wl2,
                    self.dat2[self.ids, :],
                    yerr=self.err2[self.ids, :],
                    fmt="none",
                    ecolor="grey",
                )
        if self.perr is not None:
            self.ax2.plot(wl, self.perr[self.ids, :], "r")
        # self.pbline, = self.ax2.plot(self.wlpbline,self.fff*self.dat[idf].max()*self.fp,'g')
        self.updatePassBand(remove=False)
        self.ax2.set_title(
            (("ID spaxel = %" + str(self.nsfm) + "i") % (self.ids + 1)),
            ha="left",
            position=(0.4, 1),
        )

    def ErrorSpec(self, event):
        if event.key == "e" and self.ids is not None:
            self.errcom = ~self.errcom
            if self.errcom and self.err is None:
                self.errcom = False
                print("*** No error spectra data found ***")
            if self.errcom and self.err is not None:
                wl = self.wl if self.wl.ndim == 1 else self.wl[self.ids]
                if self.ids is not None:
                    self.ax2.errorbar(
                        wl,
                        self.dat[self.ids, :],
                        yerr=self.err[self.ids, :],
                        fmt="none",
                        ecolor="grey",
                    )
                    if self.fitscom is not None and self.err2 is not None:
                        wl2 = self.wl2 if self.wl2.ndim == 1 else self.wl2[self.ids]
                        self.ax2.errorbar(
                            wl2,
                            self.dat2[self.ids, :],
                            yerr=self.err2[self.ids, :],
                            fmt="none",
                            ecolor="grey",
                        )
                self.fig2.canvas.draw()
            if not self.errcom and self.err is not None:
                self.errcom = False
                if self.ids is not None:
                    self.ax2.cla()
                    self.PlotSpec()
                    self.updatePassBand()
                    self.fig2.canvas.draw()

    def GetSpectraInfo(self, event):
        if len(self.cir) > 0:
            [item.set_visible(False) for item in self.cir]  # [item.remove() for item in self.cir]
        idl = event.artist.get_label()
        col = event.artist.get_color()
        self.ax2.set_title("Spaxel ID = %s" % str(idl), color=col)
        if str(idl) != self.sint:
            sp = self.getPatch(
                int(idl) - 1, color=col, ec=col, lw=self.slw, fill=self.sf, alpha=self.sa
            )
            self.cir.append(sp)
            self.ax.add_patch(sp)
            self.fig.canvas.draw()
        self.fig2.canvas.draw()

    def BoxFilter(self, wini=4500.0, wend=5000.0, dl=1.0):
        lf = np.arange(wini - dl, wend + dl, dl)
        ff = np.ones(lf.shape)
        ff[0] = 0.0
        ff[-1] = 0.0
        return lf, ff

    def IntFilter(
        self,
        ifi,
        fil,
        lm,
        dat,
        pasb=False,
        verb=True,
        dfil="",
        dl=None,
        center=False,
        remove_cont=True,
    ):
        def ftrapz(lm, fff, dat, ax=None):
            if dat.ndim == 3:
                # tile(lm,80*75).reshape((1970,75,80))
                lmd = lm[:, np.newaxis, np.newaxis]
                fffd = fff[:, np.newaxis, np.newaxis]
                ax = 0 if ax is None else ax
            else:
                ax = 1 if ax is None else ax
                lmd = lm
                fffd = fff
            return np.trapz(lmd * fffd * dat, lm, axis=ax) / np.trapz(fff * lm, lm)

        # Turn Off: "Warning: overflow encountered in multiply"
        np.seterr(all="ignore")
        if fil is not None:
            lf, ff = np.loadtxt(os.path.join(dfil, fil[ifi]), unpack=True, usecols=(0, 1))
        else:
            lf, ff = self.BoxFilter()
        if verb is True and fil is not None:
            print("Selected Filter: " + fil[ifi])
        iffmax = np.argmax(ff)
        if center:
            wmax = lf[iffmax]
            wcen = (lm.min() + lm.max()) / 2.0
            lf += wcen - wmax
        if dl is not None:
            lf = lf + dl
        self.wmax = lf[iffmax]
        fff = np.interp(lm, lf, ff)
        ival = ftrapz(lm, fff, dat)
        if remove_cont:
            if verb:
                print(">>> Continuum removed")
            dlw = lf.max() - lf.min()
            lfff = np.interp(lm, lf - dlw, ff)
            rfff = np.interp(lm, lf + dlw, ff)
            cval = (ftrapz(lm, lfff, dat) + ftrapz(lm, rfff, dat)) / 2.0
            ival -= cval
        if pasb is True:
            return ival, fff
        else:
            return ival

    def set_fff(self, **kwargs):
        kwargs.pop("pasb", None)
        dpars = dict(
            ifi=self.ifil,
            fil=self.list_filters,
            pasb=True,
            dfil=self.dfilter,
            center=self.cfilter,
            remove_cont=self.remove_cont,
        )
        dpars.update(kwargs)
        self.color, self.fff = self.IntFilter(lm=self.wl, dat=self.dat, **dpars)
        if self.fitscom is not None:
            dpars["verb"] = False
            self.color2, self.fff2 = self.IntFilter(lm=self.wl2, dat=self.dat2, **dpars)

    def PassBandPress(self, event):
        if self.pband is None:
            return
        tb = plt.get_current_fig_manager().toolbar
        self.dl = 0.0
        if tb.mode == "":
            if not self.pband.contains(event)[0]:
                return
            self.pressevent = event
            self.wlpbline = self.pbline[0].get_xdata()

    def PassBandMove(self, event):
        if self.pressevent is None or event.inaxes != self.pressevent.inaxes:
            return
        self.dl = event.xdata - self.pressevent.xdata
        self.pbline[0].set_xdata(self.wlpbline + self.dl)
        self.fig2.canvas.draw()

    def PassBandRelease(self, event):
        if self.pband is None:
            return
        tb = plt.get_current_fig_manager().toolbar
        if tb.mode == "" and self.dl is not None and self.pressevent is not None:
            self.pressevent = None
            self.gdl = self.gdl + self.dl
            self.set_fff(verb=False, dl=self.gdl)
            self.updatePatch()
            self.updatePassBand(remove=True)
            if self.iclm:
                if self.icmp is not None:
                    self.zmode = self.icmp.zmode
                self.icmp = IntColorMap(self.p, zmode=self.zmode)
            # self.ax2.set_xlim(self.ax2.get_xlim())
            self.fig.canvas.draw()
            self.fig2.canvas.draw()

    def FitSpec(self, event):
        if not self.fitspec and ((event.key == "i") or (event.key == "x")):
            print('*** You need module "pyspeckit" or "pyraf" to use fitting features!!! ***')
            return
        if event.key == "i":
            # 0 = PYRAF | 1 = PYSPECKIT
            if self.pyraf and not self.pyspec:
                self.fitspec_mode = 0
                print(">>> Only PyRAF module available")
            elif not self.pyraf and self.pyspec:
                self.fitspec_mode = 1
                print(">>> Only PySPECKIT module available")
            else:
                self.fitspec_mode += 1
                if self.fitspec_mode % 2 == 0:
                    print(">>> PyRAF fitting mode selected")
                else:
                    print(">>> PySPECKIT fitting mode selected")
        if event.key == "x" and len(self.list) > 0:
            ids = self.list[0]
            if self.fitspec_mode % 2 == 0:
                if plt.get_backend() != "TkAgg":
                    print(
                        '*** You need to set backend to "TkAgg" if you want to use PyRAF interactive fitting ***'
                    )
                else:
                    # sname = '%s_%s' % ('.'.join(self.name_fits.split('.')[0:-1]),self.ids.replace(',','_').replace(' ',''))
                    sname = "%s_%s" % (".".join(self.name_fits.split(".")[0:-1]), ids + 1)
                    tmpfits = tmpName(prefix="tmp_%s" % sname)
                    convert2iraf_spec(tmpfits, self.wl, self.dat[ids], title=sname)
                    print(">>> Spectrum (%s) of %s" % (ids, self.name_fits))
                    splot(tmpfits)
                    if os.path.exists(tmpfits):
                        os.remove(tmpfits)
            else:
                espec = self.err[ids] if self.err is not None else None
                sp = pyspeckit.Spectrum(xarr=self.wl, data=self.dat[ids], err=espec, header=self.hd)
                sp.plotter()
                sp.plotter.axis.set_xlabel(r"Wavelength $(\AA)$")
                sp.plotter.axis.set_ylabel(r"Flux $(\mathrm{erg/s/cm^2/\AA})$")
                sp.plotter.axis.set_title("%s (%s)" % (sp.plotter.title, ids + 1))
                plt.show()

    def WindowManager(self):
        from matplotlib.widgets import CheckButtons, Slider, Button, RadioButtons
        from matplotlib import gridspec

        if not plt.fignum_exists(3):
            self.fig3 = plt.figure(3, self.winman_size)
            self.fig3.set_label(self.winman_label)
            self.setWindowTitle(self.fig3, self.winman_label)
            # Patch Spaxel Properties
            # No podemos anyadir 'axes' antes porque no se ve, y si lo ponemos despues, los widgets No funcionan
            ps = gridspec.GridSpec(
                3, 2, left=0.2, bottom=0.8, top=0.95, wspace=0.5, hspace=0.7, right=0.8
            )
            ax = plt.subplot(ps[0, :], frameon=False, xticks=[], yticks=[])
            ax.text(0.5, 0.5, "Spaxel Properties", ha="center")
            plw = Slider(plt.subplot(ps[1, :]), "Linewidth", 0.0, 5.0, valinit=self.plw)
            pap = Slider(plt.subplot(ps[2, :]), "Alpha", 0.0, 1.0, valinit=self.palpha)
            # Patch Properties ------------------------
            cs = gridspec.GridSpec(
                3, 3, left=0.2, bottom=0.55, top=0.7, wspace=0.7, hspace=0.5, right=0.9
            )
            ax = plt.subplot(cs[0, :], frameon=False, xticks=[], yticks=[])
            ax.text(0.5, 0.5, "Spaxel Selector Properties", ha="center")
            clw = Slider(plt.subplot(cs[1, :-1]), "Linewidth", 0.0, 5.0, valinit=self.clw)
            cap = Slider(plt.subplot(cs[2, :-1]), "Alpha", 0.0, 1.0, valinit=self.ca)
            ctf = Button(plt.subplot(cs[1:, -1]), self.clab, hovercolor=self.ccc, color=self.ccc)
            # Spectra-->Spaxel Identify Properties -------------------
            ss = gridspec.GridSpec(
                3, 3, left=0.2, bottom=0.3, top=0.45, wspace=0.7, hspace=0.5, right=0.9
            )
            ax = plt.subplot(ss[0, :], frameon=False, xticks=[], yticks=[])
            ax.text(0.5, 0.5, "Spectra --> Spaxel Identifier Properties", ha="center")
            slw = Slider(plt.subplot(ss[1, :-1]), "Linewidth", 0.0, 5.0, valinit=self.slw)
            sap = Slider(plt.subplot(ss[2, :-1]), "Alpha", 0.0, 1.0, valinit=self.sa)
            stf = Button(plt.subplot(ss[1:, -1]), self.slab, hovercolor=self.scc, color=self.scc)
            # Save File Options + Colorbar ---------------------------
            svp = gridspec.GridSpec(
                3, 2, left=0.2, bottom=0.05, top=0.2, wspace=0.7, hspace=0.5, right=0.9
            )
            ax = plt.subplot(svp[0, 0], frameon=False, xticks=[], yticks=[])
            ax.text(0.5, 0.5, "Spectra Type", ha="center")
            ax = plt.subplot(svp[0, 1], frameon=False, xticks=[], yticks=[])
            ax.text(0.5, 0.5, "File Type", ha="center")
            sint = Button(plt.subplot(svp[1, 0]), self.intlab, hovercolor=self.tcc, color=self.tcc)
            sind = Button(plt.subplot(svp[2, 0]), self.indlab, hovercolor=self.icc, color=self.icc)
            stxt = Button(plt.subplot(svp[1, 1]), self.txtlab, hovercolor=self.xcc, color=self.xcc)
            sfit = Button(plt.subplot(svp[2, 1]), self.fitlab, hovercolor=self.fcc, color=self.fcc)

            def update(val):
                self.plw = plw.val
                self.p.set_linewidth(self.plw)
                self.palpha = pap.val
                self.p.set_alpha(self.palpha)
                self.clw = clw.val
                self.ca = cap.val
                self.slw = slw.val
                self.sa = sap.val
                [(item.set_alpha(self.ca), item.set_lw(self.clw)) for item in self.pat]
                [(item.set_alpha(self.sa), item.set_lw(self.slw)) for item in self.cir]
                self.fig.canvas.draw()

            plw.on_changed(update)
            pap.on_changed(update)
            clw.on_changed(update)
            cap.on_changed(update)
            slw.on_changed(update)
            sap.on_changed(update)

            def press(event):
                if event.inaxes == ctf.ax.axes:
                    self.cf = not self.cf
                    self.ccc = self.axcolor2 if self.cf else self.axcolor1
                    ctf.color, ctf.hovercolor = self.ccc, self.ccc
                    ctf.ax.set_axis_bgcolor(self.ccc)
                    self.clab = "Fill On" if self.cf else "Fill Off"
                    ctf.label.set_text(self.clab)
                    [item.set_fill(self.cf) for item in self.pat]
                    self.fig3.canvas.draw()
                    self.fig.canvas.draw()
                if event.inaxes == stf.ax.axes:
                    self.sf = not self.sf
                    self.scc = self.axcolor2 if self.sf else self.axcolor1
                    stf.color, stf.hovercolor = self.scc, self.scc
                    stf.ax.set_axis_bgcolor(self.scc)
                    self.slab = "Fill On" if self.sf else "Fill Off"
                    stf.label.set_text(self.slab)
                    [item.set_fill(self.sf) for item in self.cir]
                    self.fig3.canvas.draw()
                    self.fig.canvas.draw()
                if event.inaxes == sint.ax.axes:
                    self.integrated = not self.integrated
                    self.tcc = self.axcolor1 if not self.integrated else self.axcolor2
                    sint.color, sint.hovercolor = self.tcc, self.tcc
                    sint.ax.set_axis_bgcolor(self.tcc)
                    self.intlab = "Integrated On" if self.integrated else "Integrated Off"
                    sint.label.set_text(self.intlab)
                    self.fig3.canvas.draw()
                if event.inaxes == sind.ax.axes:
                    self.individual = not self.individual
                    self.icc = self.axcolor1 if not self.individual else self.axcolor2
                    sind.color, sind.hovercolor = self.icc, self.icc
                    sind.ax.set_axis_bgcolor(self.icc)
                    self.indlab = "Individual On" if self.individual else "Individual Off"
                    sind.label.set_text(self.indlab)
                    self.fig3.canvas.draw()
                if event.inaxes == stxt.ax.axes:
                    self.txt = not self.txt
                    self.xcc = self.axcolor1 if not self.txt else self.axcolor2
                    stxt.color, stxt.hovercolor = self.xcc, self.xcc
                    stxt.ax.set_axis_bgcolor(self.xcc)
                    self.txtlab = "Txt On" if self.txt else "Txt Off"
                    stxt.label.set_text(self.txtlab)
                    self.fig3.canvas.draw()
                if event.inaxes == sfit.ax.axes:
                    self.fits = not self.fits
                    self.fcc = self.axcolor1 if not self.fits else self.axcolor2
                    sfit.color, sfit.hovercolor = self.fcc, self.fcc
                    sfit.ax.set_axis_bgcolor(self.fcc)
                    self.fitlab = "Fits On" if self.fits else "Fits Off"
                    sfit.label.set_text(self.fitlab)
                    self.fig3.canvas.draw()

            self.fig3.canvas.mpl_connect("button_press_event", press)
            plt.show()


# ----------------------------------------------------------------------------
