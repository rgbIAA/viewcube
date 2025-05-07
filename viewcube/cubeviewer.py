############################################################################
#                              VIEWCUBE                                    #
#                              PYTHON 3                                    #
#                                                                          #
# RGB@IAA ---> Last Change: 2024/10/10                                     #
############################################################################
#
#
#
################################ VERSION ###################################
VERSION = "0.3.6"                                                          #
############################################################################
#
from matplotlib.collections import PatchCollection, PolyCollection
from .utils import lsfiles, ckfiles, LoadFits, image_max_pixel
from .utils import save_spec, convert2iraf_spec, get_min_max
from matplotlib.widgets import RectangleSelector
from matplotlib.patches import Circle, Rectangle
from setuptools._distutils.version import LooseVersion
import astropy.io.fits as pyfits
from matplotlib import rcParams
from astropy import units as u
from astropy.wcs import WCS
import matplotlib, sys, os
import numpy as np
import itertools
import argparse
import string
import random
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

# Last pylab modules to import (after pyraf)
from .rgbmpl import rnorm, IntColorMap
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def get_wavelength_coordinates(w, Nwave):
    w = w.sub([3])
    pix_coords = np.arange(Nwave)
    wave_coords = w.wcs_pix2world(pix_coords[:, np.newaxis], 0)
    if w.wcs.cunit[0] == "m":
        wave_coords *= 1e10
    return np.squeeze(wave_coords)


# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def GetIdFilter(list, filter, dfil="."):
    # if any(filter in item for item in list):
    if filter is not None:
        lfil = lsfiles("*" + filter + "*", dfil)
    else:
        lfil = []
    if len(lfil) == 0:
        if filter is not None:
            print('"' + filter + '" NOT found. Set to "' + list[0] + '"')
        else:
            print('Set filer to "' + list[0] + '"')
        return 0
    else:
        print("Filter: " + ".".join(lfil[0].split(".")[:-1]))
        return list.index(lfil[0])


# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def GetSpaxelLimits(x, y, radius):
    spax_fac = radius * 7
    xbar = abs(max(x) - min(x)) * 0.5
    ybar = abs(max(y) - min(y)) * 0.5
    xmed = xbar + min(x)
    ymed = ybar + min(y)

    xfbar = 1.2 if xbar > spax_fac else 4.0
    yfbar = 1.2 if ybar > spax_fac else 4.0

    # PPAK special cases
    if len(x) == 331 or len(x) == 993:
        yfbar = 1.3
    if len(x) == 382:
        xfbar = 1.3

    xmax_mosaic = round(xmed + xbar * xfbar)
    xmin_mosaic = round(xmed - xbar * xfbar)
    ymax_mosaic = round(ymed + ybar * yfbar)
    ymin_mosaic = round(ymed - ybar * yfbar)

    return [xmin_mosaic, xmax_mosaic, ymin_mosaic, ymax_mosaic]


# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def GetLambdaLimits(wl, pt=0.05, wlim=None):
    if isinstance(wl, (tuple, list)):
        wl = [(np.min(item), np.max(item)) for item in wl if item is not None]
    wmin = np.min(wl)
    wmax = np.max(wl)
    if wlim is not None:
        if type(wlim) not in [list, tuple] or len(wlim) != 2:
            print("Wavelength limits should be a tuple or list of two items: ex. --> (None, 6200)")
        else:
            wlimmin, wlimmax = wlim
            if wlimmin is not None:
                wmin = wlimmin
            if wlimmax is not None:
                wmax = wlimmax
    range = abs(wmax - wmin)
    return wmin - range * pt, wmax + range * pt


# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def GetFluxLimits(flim):
    fmin = None
    fmax = None
    if flim is not None:
        if type(flim) not in [list, tuple] or len(flim) != 2:
            print("Flux limits should be a tuple or list of two items: ex. --> (None, 1e-18)")
        else:
            fmin, fmax = flim
    return fmin, fmax


# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def PRectangle(x, y, r):
    # xv = x + (-r/2.,-r/2,r/2,r/2)
    # yv = y + (-r/2.,r/2,r/2,-r/2)
    if isinstance(x, (list, tuple)):
        x = np.array(x)
        y = np.array(y)
    if isinstance(x, (int, float)):
        xv = x + np.array([0.0, 0.0, r, r])
        yv = y + np.array([0.0, r, r, 0.0])
    else:
        xv = x[:, np.newaxis] + np.array([0.0, 0.0, r, r])
        yv = y[:, np.newaxis] + np.array([0.0, r, r, 0.0])
    return xv, yv


# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
def tmpName(prefix="tmp", char=8, suffix="fits"):
    schar = "".join(random.choice(string.letters + string.digits) for i in range(char))
    return "%s_%s.%s" % (prefix, schar, suffix)


# ------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
class CubeViewer:
    def __init__(
        self,
        name_fits,
        ptable=None,
        fitscom=None,
        syncom=False,
        default_filter="Halpha_KPNO-NOAO",
        exdata=None,
        exhdr=0,
        exwave=None,
        exflag=None,
        exerror=None,
        specaxis=None,
        dfilter="filters/",
        norm="sqrt",
        fo=1.0,
        fc=1.0,
        mval=0.0,
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
        cspec="#1f77b4",
        lspec=1,
        ccom="#ff7f0e",
        lcom=1,
        cflag="r",
        lflag=1,
        colorbar=True,
        fits=False,
        txt=True,
        integrated=True,
        individual=False,
        wlim=None,
        flim=None,
        iclm=True,
        fp=1.2,
        fig_spaxel_size=(7.1, 6),
        fig_spectra_size=(8, 5),
        fig_window_manager=(5, 5),
        c=299792.458,
        cfilter=False,
        remove_cont=False,
        masked=True,
        vflag=0,
        dsoni=None,
        ref_mode="crpix",
        soni_start=False,
        **kwargs
    ):
        """
        # -----------------------------------------------------------------------------
              USE:
                      rssv = CubeViewer('ngc.fits','ngc.pt.txt')
        # -----------------------------------------------------------------------------
              name_fits --> Name of the fits file
              ptbale --> Name of the position table
              fitscom = None --> Name of the RSS fits file you want to compare
              syncom = False --> If file is a Pycasso FITS file, plot also synthetic
                      spectra (as "fitscom" object)
              dfilter = 'filters/' --> Directory where the filters are located
              norm = 'sqrt' --> Normalize function for the color map
              default_filter = 'Halpha_KPNO-NOAO' --> Root name of the default filter
              fo = 1.0 --> Multiplicative factor for original file
              fc = 1.0 --> Multiplicative factor for "fitscom" file
              mval = 0.0 --> Not Finite masked values
              palpha = 0.95 --> Alpha value of the RSS spaxels
              plw = 0.1 --> Linewidth of the RSS spaxels
              plc = 'k' --> Color of the border of the RSS spaxels
              clw = 1 --> Linewidth of the RSS spaxel selector
              cf = False --> Fill value of the RSS selected spaxels
              ca = 0.8 --> Alpha value of the RSS selected spaxels
              slw = 2 --> Linewidth of the spectra-->spaxel-spaxel
              sf = False --> Fill value of the spectra-->spaxel-spaxel
              sa = 0.9 --> Alpha value of the spectra-->spaxel-spaxel
              cspec = '#1f77b4' --> Color for plotting the spectra
              lspec = 1 --> Linewidth for plotting the spectra
              ccom = '#ff7f0e' --> Color for plotting the comparison spectra
              lcom = 1 --> Linewidth for plotting the comparison spectra
              cflag = 'r' --> Color for plotting the flags
              lflag = 1 --> Linewidth for plotting the flags
              colorbar = False --> Set colorbar in the RSS spaxel viewer
              fits = False --> Save files in fits type
              txt = True --> Save files in ASCII type
              integrated = True --> Save the integrated spectrum of the selected spaxels
              individual = False --> Save the individual spectra of the selected spaxels
              wlim = None --> 2D-Tuple with the wavelength limits of the spectra inspector
              flim = None --> 2D-Tuple with the flux limits of the spectra inspector
              iclm = True --> DS9-like option for changing the dynamical range of the
                      RSS spaxel Viewer
              fp = 1.2 --> The filter is multipy by the maximum value of the spectra and
                      by this multiplicative constant
              fig_spaxel_size = (7,6) --> 2D-Tuple with the size of the RSS spaxel viewer
              fig_spectra_size = (8,5) --> 2D-Tuple with the size of the Spectral inspector
              fig_window_manager = (5,5) --> 2D-Tuple with the size of the Window Manager
              cfilter = False --> Center filter in wavelength range
              remove_cont = False --> Remove continuum from adjacent positions (right and left)
              masked = True --> Use masked arrays for flux (flag = mask)
              vflag = 0 --> Flags with values larger than "vflag" are considered flagged
              dsoni -> Directory for sonification module and data
              ref_mode -> Mode for choosing reference pixel: 'crpix' or 'max'
              soni_start -> Activate sonification mode (import libraries and check database)
        # -----------------------------------------------------------------------------
        """

        # Disable default matplotlib shorcuts
        rcParams["keymap.save"] = ""
        rcParams["keymap.yscale"] = ""
        rcParams["keymap.xscale"] = ""
        rcParams["keymap.pan"] = ""

        # Version
        self.version = "CubeViewer version: %s" % VERSION

        # Set variables
        self.name_fits = name_fits
        self.bname_fits = os.path.basename(name_fits)
        self.ptable = ptable
        self.fitscom = fitscom
        self.syncom = syncom
        self.errcom = False
        self.specres = False
        self.specfres = True
        self.dfilter = dfilter
        self.norm = norm
        self.default_filter = default_filter
        self.list_filters = None
        self.palpha = palpha
        self.plw = plw
        self.plc = plc
        self.colorbar = colorbar
        self.ix, self.iy = None, None
        self.cspec = cspec
        self.lspec = lspec
        self.ccom = ccom
        self.lcom = lcom
        self.cflag = cflag
        self.lflag = lflag
        self.rshow = False
        # Select Circle Properties
        self.clw = clw
        self.cf = cf
        self.ca = ca
        self.cc = cc
        # Select Spaxel Properties
        self.slw = slw
        self.sf = sf
        self.sa = sa
        # Color Map dynamic range
        self.iclm = iclm
        self.icmp = None
        self.zmode = False
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
        self.fo = fo
        self.fc = fc
        # Lambda
        self.wlim = wlim
        self.wrest = False
        # Save Spectra
        self.fits = fits
        self.txt = txt
        # Save type
        self.integrated = integrated
        self.individual = individual
        # Selected spectra
        self.spec = None
        self.espec = None
        self.fspec = None
        self.idl = None
        # Integrated spectra
        self.intspec = None
        self.eintspec = None
        self.view_pintspec = False
        # Sonification
        self.dsoni = dsoni
        self.soni_start = soni_start
        self.soni_mode = False
        self.sc = None
        # Constant
        self.c = c
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
        # Fit Spec
        self.fitspec = True
        self.pyraf = PYRAF
        self.pyspec = PYSPEC
        self.fitspec_mode = 0  # 0 = PYRAF | 1 = PYSPECKIT
        if not self.pyraf and not self.pyspec:
            self.fitspec = False

        # Basic Info
        # self.dat, self.hd = ofits(self.name_fits)
        self.fobj = LoadFits(
            self.name_fits,
            exdata=exdata,
            exhdr=exhdr,
            exwave=exwave,
            exflag=exflag,
            exerror=exerror,
            specaxis=specaxis,
            **kwargs
        )
        self.K = self.fobj.K
        self.dat = self.fo * self.fobj.data
        self.syn = self.fobj.syn
        self.err = self.fobj.error
        self.flag = self.fobj.flag
        self.crval = self.fobj.crval
        self.cdelt = self.fobj.cdelt
        self.redshift = self.fobj.redshift
        self.velocity = self.fobj.velocity
        self.pversion = self.fobj.pversion
        self.zones = None
        self.perr = None
        self.pres = None
        if self.err is not None:
            self.err *= self.fo
        if self.syn is not None:
            self.syn *= self.fo
        if self.flag is not None:
            if self.flag.dtype.kind != "b":
                self.flag = np.where(self.flag > vflag, True, False)
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
        shp = np.shape(self.dat)
        naxes = len(shp)
        self.root = ".".join(self.name_fits.split(".")[:-1])
        # min, max = self.dat.min(), self.dat.max()
        self.fobj2 = None
        self.wl2 = None
        self.err2 = None
        if self.fitscom is not None:
            self.fobj2 = LoadFits(
                self.fitscom,
                exdata=exdata,
                exhdr=exhdr,
                exwave=exwave,
                exflag=exflag,
                exerror=exerror,
                specaxis=specaxis,
                **kwargs
            )
            self.wl2 = self.fobj2.wave
            self.dat2 = self.fc * self.fobj2.data
            if self.fobj2.error is not None:
                self.err2 = self.fc * self.fobj2.error
        if self.syncom and self.syn is not None:
            self.wl2 = self.fobj.wave
            self.dat2 = self.fobj.syn * self.fc
        if self.syn is not None:
            if LooseVersion(self.pversion) < LooseVersion("2.0.0"):
                self.K.zres = (self.K.f_obs - self.K.f_syn) / self.K.fobs_norm
                self.K.tres = (self.K.f_obs - self.K.f_syn).sum(axis=1) / self.K.fobs_norm.sum()
                self.zones = self.fobj.K.qZones.copy()
                self.K.res = self.K.zoneToYX(self.K.zres, extensive=False)
            else:
                if self.K.hasSegmentationMask:
                    self.K.zres = (
                        (self.K.f_obs - self.K.f_syn) * self.K.flux_unit / self.K.fobs_norm
                    )
                    self.K.tres = (self.K.flux_unit * (self.K.f_obs - self.K.f_syn)).sum(
                        axis=1
                    ) / self.K.fobs_norm.sum()
                    self.zones = self.K.spatialize(
                        np.arange(1, self.K.segmentationMask.shape[0] + 1), extensive=False
                    )
                    self.K.res = self.K.spatialize(self.K.zres, extensive=False)
                    self.K.fres = self.K.spatialize(self.K.zres * self.K.fobs_norm, extensive=False)
                else:
                    self.K.res = (
                        (self.K.f_obs - self.K.f_syn).filled(0.0)
                        * self.K.flux_unit
                        / self.K.fobs_norm
                    )
                    self.K.fres = (self.K.f_obs - self.K.f_syn).filled(0.0) * self.K.flux_unit
                    self.K.zres = (
                        (self.K.f_obs - self.K.f_syn).filled(0.0)[:, ~self.K.synthImageMask]
                        * self.K.flux_unit
                        / self.K.fobs_norm[~self.K.synthImageMask]
                    )
                    self.K.tres = (
                        (self.K.flux_unit * (self.K.f_obs - self.K.f_syn))[
                            :, ~self.K.synthImageMask
                        ]
                    ).sum(axis=1) / (self.K.fobs_norm[~self.K.synthImageMask]).sum()
                    self.zones = np.zeros(self.K.synthImageMask.shape)
                    self.zones[~self.K.synthImageMask] = np.arange(
                        1, (~self.K.synthImageMask).sum() + 1
                    )
        # Lambda
        if naxes == 2:
            self.ns = shp[0]  # Number of spectra
            self.nl = shp[1]  # Number of lambdas
        if naxes == 3:
            self.xx = shp[2]  # X
            self.yy = shp[1]  # Y
            self.nl = shp[0]  # Number of lambdas
            try:
                self.cdelt1 = self.hd["CDELT1"]
            except:
                self.cdelt1 = 1.0
            try:
                self.cdelt2 = self.hd["CDELT2"]
            except:
                self.cdelt2 = 1.0
            try:
                self.crpix1 = self.hd["CRPIX1"]
                self.crpix2 = self.hd["CRPIX2"]
            except:
                self.crpix1 = int(self.xx / 2.0)
                self.crpix2 = int(self.yy / 2.0)
            # crval = self.hd['CRVAL3']
            # cdelt = self.hd['CDELT3']
            self.fx = len(str(self.xx))  # Number of digits of the number of pixels in X axis
            self.fy = len(str(self.yy))  # Number of digits of the number of pixels in Y axis
            #  Spatial resolution
            self.sr = 1.0 if (self.cdelt1 == 0 or self.cdelt2 == 0) else float(self.cdelt1)
            self.wcs = WCS(self.hd).celestial.wcs
            if "deg" in self.wcs.cunit:
                self.sr = (self.cdelt1 * u.deg).to("arcsec").value
            # Choose reference pixel and ext
            self.get_ref_pix(ref_mode)
        # Lambda
        # self.wl = crval + arange(self.nl)*cdelt
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
        if self.fitscom is not None:
            self.wlmin = min(self.wl.min(), self.wl2.min())
            self.wlmax = max(self.wl.max(), self.wl2.max())
        else:
            self.wlmin = self.wl.min()
            self.wlmax = self.wl.max()
        self.pbline = None
        self.pband = None
        # Data for color imshow map
        self.dcolor = self.dat

        # Filtros
        if self.dfilter is not None:
            self.list_filters = lsfiles("*txt", self.dfilter, path=False)
        self.ifil = None
        if self.list_filters is not None:
            self.ifil = GetIdFilter(self.list_filters, self.default_filter, self.dfilter)
            self.nfil = len(self.list_filters)
        self.set_fff(verb=False)
        self.gdl = 0  # Global delta position filter
        # color = dat[:,0:1900].sum(axis=1)

        # Read Position Table
        if self.ptable is not None:
            ckfiles(self.ptable)
            with open(self.ptable, "r") as f:  # Cierra automaticamente el fichero
                mt, xs, ys, eid = f.readline().split()

            if mt == "C":
                self.radius = float(xs)
                self.fiber_size = 2.0 * self.radius

            self.id, self.x, self.y, self.flag = np.loadtxt(self.ptable, unpack=True, skiprows=1)
            self.nsfm = len(str(self.ns))  # Number of digits of number of spaxels

        else:
            pass

        # Create Figure 2 (Spectral Viewer) ---------------------------------
        self.fig2 = plt.figure(2, self.fig2_size)
        self.fig2.set_label(self.fig2_label)
        #self.setWindowTitle(self.fig2, self.fig2_label)
        self.fig2.canvas.manager.set_window_title(self.fig2_label)
        self.ax2 = self.fig2.add_subplot(111)
        self.awlmin, self.awlmax = GetLambdaLimits((self.wl, self.wl2), 0.05, wlim=self.wlim)
        self.fmin, self.fmax = GetFluxLimits(self.flim)
        self.ax2.set_xlim((self.awlmin, self.awlmax))
        self.ax2.set_ylim((self.fmin, self.fmax))
        # self.ax2.xaxis.get_major_formatter()
        # self.fig2.canvas.set_window_title("Spectral Viewer by RGB")

        # Create Figure 1 (Spaxel Viewer) ---------------------------------
        # self.exts = [self.ext[0]-0.5,self.ext[1]+0.5,self.ext[2]-0.5,self.ext[3]+0.5]
        self.fig = plt.figure(1, self.fig1_size)
        self.fig.set_label(self.fig1_label)
        #self.setWindowTitle(self.fig, self.fig1_label)
        self.fig.canvas.manager.set_window_title(self.fig1_label)
        self.ax = self.fig.add_subplot(111)
        self.p = self.ax.imshow(
            self.color,
            alpha=self.palpha,
            extent=self.ext,
            norm=rnorm(self.norm),
            interpolation="nearest",
            aspect="auto",
            origin="lower",
        )
        self.ax.axis(self.ext)
        self.ax.set_title(self.bname_fits)
        # self.ax.set_xlim(self.ax.get_xlim()[::-1])
        # self.ax = self.fig.add_subplot(111,frameon=False,xticks=[],yticks=[])
        # self.ax.axis([0,self.xx,0,self.yy])
        # plt.setp(self.par.get_yticklabels(),visible=False)
        # for item in self.ax.get_xticklabels():
        # self.ax.set_xticks((1,35))
        # self.ax.set_xticklabels(('pepe','cosa'))
        #  self.xmin,self.xmax,self.ymin,self.ymax = GetSpaxelLimits(self.x,self.y,self.radius)
        #  self.ax.axis([self.xmax,self.xmin,self.ymin,self.ymax])

        #  self.patches = [Circle((self.x[i],self.y[i]),self.radius) for i in range(self.ns)]
        #  self.p = PatchCollection(self.patches,cmap=plt.cm.gray_r,norm=rnorm(self.norm),
        # 	alpha=self.palpha,lw=self.plw,color=self.plc)
        #  self.p.set_array(self.color)
        #  self.ax.add_collection(self.p)
        if self.colorbar:
            self.fig.colorbar(self.p)
        self.cmin, self.cmax = get_min_max(self.color)
        self.p.set_clim([self.cmin, self.cmax])

        # Bind Events
        # Se pone el primero SpectraViewer con Motion para que el mensaje
        # informativo de la posicion de la toolbar se cambie
        self.fig.canvas.mpl_connect("motion_notify_event", self.SpectraViewer)
        self.fig.canvas.mpl_connect("key_press_event", self.PressKey)
        self.fig2.canvas.mpl_connect("key_press_event", self.PressKey)
        self.fig.canvas.mpl_connect("key_press_event", self.ChangeFilter)
        self.fig2.canvas.mpl_connect("key_press_event", self.ChangeFilter)
        self.fig.canvas.mpl_connect("button_press_event", self.SpectraViewer)
        self.fig2.canvas.mpl_connect("button_press_event", self.SpectraViewer)
        self.fig2.canvas.mpl_connect("pick_event", self.GetSpectraInfo)
        self.fig2.canvas.mpl_connect("button_press_event", self.PassBandPress)
        self.fig2.canvas.mpl_connect("motion_notify_event", self.PassBandMove)
        self.fig2.canvas.mpl_connect("button_release_event", self.PassBandRelease)
        # Argument drawtype='box in RectableSelector is deprecated in 3.5.0 as the only behaviour is box
        self.RS = RectangleSelector(self.ax, self.onselect, button=1)
        self.fig.canvas.mpl_connect("key_press_event", self.toggle_selector)
        self.fig.canvas.mpl_connect("key_press_event", self.ChangeSpaxelViewer)
        self.fig2.canvas.mpl_connect("key_press_event", self.Redshift)
        self.fig2.canvas.mpl_connect("key_press_event", self.SynthSpec)
        self.fig.canvas.mpl_connect("key_press_event", self.SynthSpec)
        self.fig2.canvas.mpl_connect("key_press_event", self.ErrorSpec)
        self.fig.canvas.mpl_connect("key_press_event", self.ErrorSpec)
        self.fig2.canvas.mpl_connect("key_press_event", self.ResSpec)
        self.fig.canvas.mpl_connect("key_press_event", self.ResSpec)
        self.fig2.canvas.mpl_connect("key_press_event", self.RestWave)
        self.fig2.canvas.mpl_connect("key_press_event", self.FitSpec)
        self.fig.canvas.mpl_connect("key_press_event", self.GetPycassoRes)
        self.fig2.canvas.mpl_connect("key_press_event", self.GetPycassoRes)
        self.fig.canvas.mpl_connect("key_press_event", self.selectZone)
        self.RS.set_active(False)

        # DS9-like dynamic range of colormap
        if self.iclm:
            self.icmp = IntColorMap(self.p)
        plt.show()

    def res2pix(self, x, y):
        return int((x / self.sr) + self.x_ref), int((y / self.sr) + self.y_ref)

    def pix2res(self, x, y):
        return (x - self.x_ref) * self.sr, (y - self.y_ref) * self.sr

    def vredshift(self, cz):
        return self.wl / (1.0 + (cz / self.c))

    def get_ref_pix(self, mode="max", **kwargs):
        if mode.lower() == "max":
            image = np.nanmedian(self.dat, axis=0)
            self.y_ref, self.x_ref = image_max_pixel(image, **kwargs)
        else:
            self.y_ref = self.crpix2
            self.x_ref = self.crpix1
        self.ext = [
            -self.x_ref * self.sr,
            (self.xx - self.x_ref) * self.sr,
            -self.y_ref * self.sr,
            (self.yy - self.y_ref) * self.sr,
        ]
        # self.ext = [-self.ext[0],-self.ext[1],self.ext[2],self.ext[3]]
        # self.ext = [0,self.xx,0,self.yy]

    def onselect(self, eclick, erelease):
        # 'eclick and erelease are matplotlib events at press and release'
        ix, iy = self.res2pix(eclick.xdata, eclick.ydata)
        fx, fy = self.res2pix(erelease.xdata, erelease.ydata)
        # print ' used button   : ', eclick.button

        ix, fx = sorted([ix, fx])
        iy, fy = sorted([iy, fy])
        si = itertools.product(np.arange(ix, fx), np.arange(iy, fy))

        # TODO
        # self.ax.add_collection(PatchCollection(self.r,alpha=0.1))
        # x,y=self.pix2res(np.tile(arange(ix,fx),fy-iy),np.repeat(np.arange(iy,fy),fx-ix))
        # xr,yr=PRectangle(x,y,self.sr)
        # pol=PolyCollection(zip(xr,yr))
        # self.ax.add_collection([po])
        # rec = self.ax.fill(xl,yl,fill=self.cf,color=self.cc,ec=self.cc,alpha=self.ca,lw=self.clw)
        # si=itertools.product(arange(ix,fx),arange(iy,fy))
        # cx=[];cy=[]
        # for x,y in si: cx.append(x); cy.append(y)
        # cx = array(cx); cy = array(cy)
        # cx,cy=self.pix2res(cx,cy)
        # xr,yr = PRectangle(cx,cy,self.sr)
        # yl = []; xl=[]
        # for mx,my in zip(xr,yr):
        # xl.extend(mx)
        # xl.append(None)
        # yl.extend(my)
        # yl.append(None)
        # rec = self.ax.fill(xl,yl,fill=self.cf,color=self.cc,ec=self.cc,alpha=self.ca,lw=self.clw)
        # rec[0].remove()
        # for x in range(ix,fx):
        # for y in range(iy,fy):
        #  cx,cy = self.pix2res(x,y)
        #  xr,yr = PRectangle(cx,cy,self.sr)
        #  rec = self.ax.fill(xr,yr,fill=self.cf,color=self.cc,ec=self.cc,alpha=self.ca,lw=self.clw)
        #  self.list.append((x,y))
        #  self.pat.append(rec)
        self.fig.canvas.draw()

    def toggle_selector(self, event):
        if event.key in ["B"] and self.RS.active:
            print(" RectangleSelector deactivated.")
            self.RS.set_active(False)
        if event.key in ["b"] and not self.RS.active:
            print(" RectangleSelector activated.")
            self.RS.set_active(True)

    def PressKey(self, event):
        if event.key == "*":
            self.ax2.cla()
            self.ax.cla()
            norm = rnorm(self.norm) if self.icmp.scale is None else rnorm(self.icmp.scale)
            ne = self.p.get_clim()
            idc = self.icmp.cmindx
            self.p = self.ax.imshow(
                self.color,
                alpha=self.palpha,
                extent=self.ext,
                cmap=self.p.get_cmap(),
                norm=norm,
                interpolation="nearest",
                origin="lower",
                aspect="auto",
            )
            self.p.set_clim(ne)
            # self.ax = self.fig.add_subplot(111,frameon=False,xticks=[],yticks=[])
            # plt.setp(self.axp.get_yticklabels(),visible=False)
            # plt.setp(self.axp.get_xticklabels(),visible=False)
            # self.ax.axis([0,self.xx,0,self.yy])
            self.ax.axis(self.ext)
            self.ax.set_title(self.bname_fits)
            if self.iclm:
                self.icmp = IntColorMap(self.p)
                self.icmp.cmindx = idc
            self.list = []
            self.pbline = None
            self.fig.canvas.draw()
            self.fig2.canvas.draw()
        if event.key == "s":
            self.mode = not self.mode
        # Sonification
        if event.key == "h":
            if not self.soni_start:
                self.soni_start = True
                self.Sonification()
            self.soni_mode = not self.soni_mode
            if self.sc is None:
                self.soni_mode = False
            if not self.soni_mode and self.sc is not None and self.sc.cs is not None:
                self.sc.stop_sound()
        # Save selected Spectra
        if len(self.list) > 0 and event.key == "S":
            self.SaveFile()
        if event.key == "w":
            self.WindowManager()
        if event.key == "l":
            self.LambdaLimits()
        if event.key == "Y":
            self.FluxLimits()
        if event.key == "I" and self.fobj.K is not None:
            self.view_pintspec = not self.view_pintspec
        if event.key == "q":
            if self.sc is not None and self.sc.cs is not None:
                self.sc.close_sound()
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
                        'Integrated Spectra extracted from: "' + self.name_fits + '"',
                        "Sum of Spaxels (ID): " + " | ".join([str(id) for id in self.list]),
                    ]
                    infohd = [
                        ["3DVWR_1", infotxt[0]],
                        ["3DVWR_2", infotxt[1]],
                        ["CRVAL1", self.crval],
                        ["CDELT1", self.cdelt],
                    ]
                    self.intspec = np.array([self.dat[:, y, x] for x, y in self.list]).sum(0)
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
                            'Spectra extracted from: "' + self.name_fits + '"',
                            "Spaxel (ID): " + str(item),
                        ]
                        infohd = [
                            ["3DVWR_1", infotxt[0]],
                            ["3DVWR_2", infotxt[1]],
                            ["CRVAL1", self.crval],
                            ["CDELT1", self.cdelt],
                        ]
                        x, y = item
                        strl = ("%" + str(self.fx) + "i_%" + str(self.fy) + "i") % (x, y)
                        save_spec(
                            self.wl,
                            self.dat[:, y, x],
                            fname + "_" + strl,
                            fits=self.fits,
                            hd=self.hd,
                            txt=self.txt,
                            infohd=infohd,
                            infotxt=infotxt,
                        )
                print("Files Saved")

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

    def Redshift(self, event):
        if event.key == "k":
            cz = input("Enter Redshift in km/s (Enter to abort): ")
            if len(cz) == 0:
                print("No Redshift provided")
            else:
                velocity = float(cz)
                self.wl = self.vredshift(velocity)
                if self.velocity is None:
                    self.velocity = velocity
                if self.redshift is None:
                    self.redshift = velocity / self.c
                if self.orig_wl_rest is None:
                    self.orig_wl_rest = self.wl
                if self.ix is not None and self.iy is not None:
                    self.set_fff(verb=False, dl=self.gdl)
                    self.PlotSpec()
                    self.updatePatch()
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
            if self.ix is not None and self.iy is not None:
                self.set_fff(verb=False, dl=self.gdl)
                self.PlotSpec()
                self.updatePatch()
                self.fig.canvas.draw()
                self.fig2.canvas.draw()
                self.wrest = ~self.wrest

    def SynthSpec(self, event):
        if event.key == "p":
            self.syncom = not self.syncom
            if self.specres:
                self.syncom = False
                return
            if self.syncom and self.syn is None:
                print("*** No synthetic spectra data found ***")
            if self.syncom and self.syn is not None:
                self.fitscom = "Synthetic Spectra"
                self.wl2 = self.wl
                self.dat2 = self.syn
                if self.ix is not None and self.iy is not None:
                    self.ax2.plot(
                        self.wl2, self.dat2[:, self.iy, self.ix], c=self.ccom, lw=self.lcom
                    )
                self.fig2.canvas.draw()
            if not self.syncom and self.syn is not None:
                self.fitscom = None
                if self.ix is not None and self.iy is not None:
                    self.ax2.cla()
                    self.ax2.plot(
                        self.wl, self.dat[:, self.iy, self.ix], c=self.cspec, lw=self.lspec
                    )
                    self.updatePassBand()
                    self.fig2.canvas.draw()

    def ErrorSpec(self, event):
        if event.key == "e" and self.ix is not None and self.iy is not None:
            self.errcom = ~self.errcom
            if self.errcom and self.err is None:
                self.errcom = False
                print("*** No error spectra data found ***")
            if self.errcom and self.err is not None:
                if self.ix is not None and self.iy is not None:
                    self.ax2.errorbar(
                        self.wl,
                        self.dat[:, self.iy, self.ix],
                        yerr=self.err[:, self.iy, self.ix],
                        fmt="none",
                        ecolor="grey",
                    )
                    if self.fitscom is not None and self.err2 is not None:
                        self.ax2.errorbar(
                            self.wl2,
                            self.dat2[:, self.iy, self.ix],
                            yerr=self.err2[:, self.iy, self.ix],
                            fmt="none",
                            ecolor="grey",
                        )
                self.fig2.canvas.draw()
            if not self.errcom and self.err is not None:
                self.errcom = False
                if self.ix is not None and self.iy is not None:
                    self.ax2.cla()
                    self.PlotSpec()
                    self.updatePassBand()
                    self.fig2.canvas.draw()

    def ResSpec(self, event):
        if event.key == "r" or event.key == "u":
            if event.key == "r":
                self.specres = not self.specres
            if event.key == "u":
                self.specfres = not self.specfres
            if self.specres and self.K is None:
                self.specres = False
                print("*** No residual spectra data found ***")
            if self.specres:
                if self.specfres:
                    self.res = self.K.fres
                else:
                    self.res = self.K.res
                self.pres = self.res.copy()
                self.pres[~self.flag] = np.nan
                self.dcolor = self.res
            else:
                self.dcolor = self.dat
            self.ax2.cla()
            self.PlotSpec()
            self.updatePassBand()
            self.fig2.canvas.draw()
            self.updateAx1()

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
        remove_cont=False,
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
        self.color, self.fff = self.IntFilter(lm=self.wl, dat=self.dcolor, **dpars)
        if self.fitscom is not None and "Synthetic" not in self.fitscom:
            dpars["verb"] = False
            self.color2, self.fff2 = self.IntFilter(lm=self.wl2, dat=self.dat2, **dpars)

    def updatePatch(self):
        color = self.color
        if not (self.fitscom is None or (self.fitscom is not None and "Synthetic" in self.fitscom)):
            if (self.wl[0] <= self.wmax) and (self.wmax <= self.wl[-1]):
                color = self.color
            if (self.wl2[0] <= self.wmax) and (self.wmax <= self.wl2[-1]):
                color = self.color2
        self.cmin, self.cmax = get_min_max(color)
        self.p.set_array(color)
        self.p.set_clim([self.cmin, self.cmax])

    def updatePassBand(self, remove=False):
        if self.pbline is not None and remove:
            self.pbline.pop(0).remove()
            self.pbline = None
            self.pband.remove()
            self.pband = None
        if self.specres:
            self.pbline = self.ax2.plot(
                self.wl, self.fff * self.res[:, self.iy, self.ix].max() * self.fp, "g"
            )
            self.pband = self.ax2.fill_between(
                self.wl,
                self.fff * self.res[:, self.iy, self.ix].max() * self.fp,
                color="g",
                alpha=0.25,
            )
        else:
            if self.fitscom is None or (self.fitscom is not None and "Synthetic" in self.fitscom):
                self.pbline = self.ax2.plot(
                    self.wl, self.fff * self.dat[:, self.iy, self.ix].max() * self.fp, "g"
                )
                self.pband = self.ax2.fill_between(
                    self.wl,
                    self.fff * self.dat[:, self.iy, self.ix].max() * self.fp,
                    color="g",
                    alpha=0.25,
                )
            else:
                if (self.wl[0] <= self.wmax) and (self.wmax <= self.wl[-1]):
                    self.pbline = self.ax2.plot(
                        self.wl, self.fff * self.dat[:, self.iy, self.ix].max() * self.fp, "g"
                    )
                    self.pband = self.ax2.fill_between(
                        self.wl,
                        self.fff * self.dat[:, self.iy, self.ix].max() * self.fp,
                        color="g",
                        alpha=0.25,
                    )
                elif (self.wl2[0] <= self.wmax) and (self.wmax <= self.wl2[-1]):
                    self.pbline = self.ax2.plot(
                        self.wl2, self.fff2 * self.dat2[:, self.iy, self.ix].max() * self.fp, "g"
                    )
                    self.pband = self.ax2.fill_between(
                        self.wl2,
                        self.fff2 * self.dat2[:, self.iy, self.ix].max() * self.fp,
                        color="g",
                        alpha=0.25,
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

    def SpectraViewer(self, event):
        tb = plt.get_current_fig_manager().toolbar
        # plt.get_current_fig_manager().toolbar.set_message('Pepep')
        if event.inaxes and tb.mode == "" and not self.RS.active:
            # s = event.inaxes.format_coord(event.xdata, event.ydata)
            x, y = event.xdata, event.ydata
            ix, iy = self.res2pix(x, y)
            # ix,iy = int(round(event.xdata)), int(round(event.ydata)) #self.res2pix(x,y)
            # self.ax.transLimits.transform((x,y))
            # inv = self.ax.transData.inverted(); inv.transform((event.xdata,event.ydata))
            if (
                (0 <= ix < self.xx)
                and (0 <= iy < self.yy)
                and not isinstance(self.p._A[iy, ix], np.ma.core.MaskedConstant)
            ):
                self.fig.canvas.toolbar.set_message(
                    "(RA,DEC) = (%4.2f,%4.2f) | (x,y) = (%2i,%2i) | z = %.2f"
                    % (x, y, ix, iy, self.p._A[iy, ix])
                )
                # Selection of Spectra --------------------------------
                if event.key == "d":
                    if (ix, iy) in self.list:
                        idp = self.list.index((ix, iy))
                        self.pat[idp][0].remove()
                        self.list.remove((ix, iy))
                        self.pat.pop(idp)
                        self.fig.canvas.draw()
                if event.button == 1:
                    if len(self.cir) > 0:
                        [item[0].set_visible(False) for item in self.cir]  # item.remove() ?
                    if event.key == "d":
                        if (ix, iy) in self.list:
                            idp = self.list.index((ix, iy))
                            self.pat[idp][0].remove()
                            self.list.remove((ix, iy))
                            self.pat.pop(idp)
                    else:
                        if (ix, iy) not in self.list:
                            ir, id = self.pix2res(ix, iy)
                            xr, yr = PRectangle(ir, id, self.sr)
                            rec = self.ax.fill(
                                xr,
                                yr,
                                fill=self.cf,
                                color=self.cc,
                                ec=self.cc,
                                alpha=self.ca,
                                lw=self.clw,
                            )
                            self.list.append((ix, iy))
                            self.pat.append(rec)
                    self.fig.canvas.draw()
                # Mode Spectra Viewer ('s') ---------------------------------
                if self.mode:
                    self.ix, self.iy = ix, iy
                    self.PlotSpec()
                    self.fig2.canvas.draw()
                if (
                    self.zones is not None
                    and self.view_pintspec
                    and (
                        self.zones[iy, ix] <= 0
                        or isinstance(self.zones[iy, ix], np.ma.core.MaskedConstant)
                    )
                ):
                    self.PlotPycassoIntSpec()
                    self.fig2.canvas.draw()
            # Sonification ('h') ----------------------------------
            if self.soni_mode and self.sc is not None and self.sc.cs is not None:
                if (
                    (0 <= ix < self.xx)
                    and (0 <= iy < self.yy)
                    and not isinstance(self.p._A[iy, ix], np.ma.core.MaskedConstant)
                ):
                    self.sc.sonify(self.iy, self.ix, self.spec)
                else:
                    self.sc.stop_sound()
            # View selected Spectra  -------------------------------
            # Right click = 3
            if len(self.list) > 0 and event.button == 3:  # and event.key == 'shift':
                self.pband = None
                # Change mode of spectra view if we are in the SpectraViewer Figure
                if event.canvas.figure.get_label() == self.fig2_label:
                    self.spec_mode += 1
                self.ax2.cla()
                # Set selection spectra mode so other modes like error visualization do not activate
                self.iy = None
                self.ix = None
                if (self.spec_mode % 3) == 0 or (self.spec_mode % 3) == 1:
                    for x, y in self.list:
                        p = self.ax2.plot(
                            self.wl, self.dat[:, y, x], label=str(x) + "," + str(y), picker=True
                        )
                        if self.fitscom is not None:
                            self.ax2.plot(
                                self.wl2,
                                self.dat2[:, y, x],
                                label=str(x) + "," + str(y),
                                picker=True,
                                c=p[0].get_color(),
                                alpha=0.7,
                            )
                if (self.spec_mode % 3) == 1 or (self.spec_mode % 3) == 2:
                    self.intspec = np.ma.array([self.dat[:, y, x] for x, y in self.list]).sum(0)
                    self.eintspec = np.sqrt(
                        np.ma.array([self.err[:, y, x] ** 2 for x, y in self.list]).sum(0)
                    )
                    p = self.ax2.plot(self.wl, self.intspec, label=self.sint, picker=True)
                    if self.fitscom is not None:
                        intspec = np.ma.array([self.dat2[:, y, x] for x, y in self.list]).sum(0)
                        self.ax2.plot(
                            self.wl2,
                            intspec,
                            label=self.sint,
                            picker=True,
                            c=p[0].get_color(),
                            alpha=0.7,
                        )
                self.fig2.canvas.draw()
        if not event.inaxes and self.soni_mode and self.sc is not None and self.sc.cs is not None:
            self.sc.stop_sound()

    def getSpec(self, ix, iy):
        if ix is None or iy is None:
            return
        if self.specres:
            self.spec = self.res[:, iy, ix]
        else:
            self.spec = self.dat[:, iy, ix]
        self.espec = self.err[:, iy, ix] if self.err is not None else None
        self.fspec = self.flag[:, iy, ix] if self.flag is not None else None
        self.slabel = ("Spaxel ID = %" + str(self.fx) + "i , %" + str(self.fy) + "i") % (ix, iy)

    def PlotSpec(self):
        if self.ix is None or self.iy is None:
            return
        self.ax2.cla()
        self.getSpec(self.ix, self.iy)
        if self.zones is not None and not isinstance(
            self.zones[self.iy, self.ix], np.ma.core.MaskedConstant
        ):
            self.slabel = "%s (#%i)" % (self.slabel, self.zones[self.iy, self.ix])
        # self.ilab = ('%'+str(self.fx)+'i , %'+str(self.fy)+'i') % (self.ix,self.iy)
        self.ax2.plot(
            self.wl, self.spec, c=self.cspec, lw=self.lspec, label=self.slabel, picker=True
        )
        self.ax2.set_title(self.slabel, ha="left", position=(0.35, 1))
        if self.fitscom is not None and not self.specres:
            if self.dat2.ndim == 3:
                self.ax2.plot(
                    self.wl2,
                    self.dat2[:, self.iy, self.ix],
                    c=self.ccom,
                    lw=self.lcom,
                    label=self.fitscom,
                )
                if self.errcom and self.err2 is not None:
                    self.ax2.errorbar(
                        self.wl2,
                        self.dat2[:, self.iy, self.ix],
                        yerr=self.err2[:, self.iy, self.ix],
                        fmt="none",
                        ecolor="grey",
                    )
            if self.dat2.ndim == 1:
                self.ax2.plot(self.wl2, self.dat2, c=self.ccom, lw=self.lcom, label=self.fitscom)
                if self.errcom and self.err2 is not None:
                    self.ax2.errorbar(
                        self.wl2, self.dat2, yerr=self.err2, fmt="none", ecolor="grey"
                    )
        if self.errcom and self.espec is not None:
            self.ax2.errorbar(self.wl, self.spec, yerr=self.espec, fmt="none", ecolor="grey")
        if self.perr is not None and not self.specres:
            self.ax2.plot(self.wl, self.perr[:, self.iy, self.ix], c=self.cflag, lw=self.lflag)
        if self.pres is not None and self.specres:
            self.ax2.plot(self.wl, self.pres[:, self.iy, self.ix], c=self.cflag, lw=self.lflag)
        self.updatePassBand(remove=False)

    def PlotPycassoIntSpec(self):
        self.ax2.cla()
        strid = "Pycasso Integrated Spectra (Total: %i zones)" % np.max(self.zones)
        if self.fobj.pycasso == 1:
            self.spec = self.fobj.K.integrated_f_obs.copy() * self.fo
            self.espec = self.fobj.K.integrated_f_err.copy() * self.fo
            self.fspec = self.fobj.K.integrated_f_flag.copy().astype(bool)
        else:
            self.spec = self.fobj.K.integ_f_obs.copy() * self.fobj.K.flux_unit * self.fo
            self.espec = self.fobj.K.integ_f_err.copy() * self.fobj.K.flux_unit * self.fo
            self.fspec = self.fobj.K.integ_f_flag.copy().astype(bool)
        self.espec[self.fspec] = np.nan
        self.ax2.plot(self.wl, self.spec, c=self.cspec, lw=self.lspec, label=strid, picker=True)
        self.ax2.set_title(strid, ha="center")
        if self.fitscom is not None and self.fobj2 is not None and self.fobj2.K is not None:
            if self.fobj.pycasso == 1:
                cspec = self.fobj2.K.integrated_f_obs * self.fc * self.fobj.K.flux_unit
                cspec[self.fobj2.K.integrated_f_flag.astype(bool)] = np.nan
            else:
                cspec = self.fobj2.K.integ_f_obs * self.fc
                cspec[self.fobj2.K.integ_f_flag.astype(bool)] = np.nan
            self.ax2.plot(self.wl2, cspec, c=self.ccom, lw=self.lcom, label=self.fitscom)
        if self.syncom:
            if self.fobj.pycasso == 1:
                cspec = self.fobj.K.integrated_f_syn.copy() * self.fc
            else:
                cspec = self.fobj.K.integ_f_syn.copy() * self.fc * self.fobj.K.flux_unit
            self.ax2.plot(self.wl, cspec, c=self.ccom, lw=self.lcom, label=self.fitscom)
        if self.espec is not None and self.errcom:
            self.ax2.errorbar(self.wl, self.spec, yerr=self.espec, fmt="none", ecolor="grey")
        if self.fobj.pycasso == 1:
            fintspec = self.fobj.K.integrated_f_obs.copy()
        else:
            fintspec = self.fobj.K.integ_f_obs.copy() * self.fobj.K.flux_unit
        fintspec[~self.fspec] = np.nan
        self.ax2.plot(self.wl, fintspec, c=self.cflag, lw=self.lflag)
        self.pbline = self.ax2.plot(self.wl, self.fff * self.spec.max() * self.fp, "g")
        self.pband = self.ax2.fill_between(
            self.wl, self.fff * self.spec.max() * self.fp, color="g", alpha=0.25
        )
        self.ax2.set_xlim((self.awlmin, self.awlmax))
        self.ax2.set_ylim((self.fmin, self.fmax))

    def GetSpectraInfo(self, event):
        if self.view_pintspec:
            return
        if len(self.cir) > 0:
            [
                item[0].set_visible(False) for item in self.cir
            ]  # [item.remove() for item in self.cir]
        self.idl = event.artist.get_label()
        col = event.artist.get_color()
        if self.idl.find("=") != -1:
            self.idl = self.idl.split("=")[1].split("(")[0]
        stit = self.idl
        if str(self.idl) != self.sint and "," in self.idl:
            ilx, ily = list(map(int, self.idl.split(",")))
            lx, ly = self.pix2res(ilx, ily)
            xr, yr = PRectangle(lx, ly, self.sr)
            rec = self.ax.fill(xr, yr, fill=self.sf, color=col, ec=col, lw=self.slw, alpha=self.sa)
            self.cir.append(rec)
            self.spec = self.dat[:, ily, ilx]
            self.espec = self.err[:, ily, ilx] if self.err is not None else None
            self.fspec = self.flag[:, ily, ilx] if self.flag is not None else None
            self.fig.canvas.draw()
            if self.zones is not None:
                stit = "%s (#%i)" % (stit, self.zones[ily, ilx])
        else:
            self.spec = self.intspec
            self.espec = self.eintspec
        stit = "Spaxel ID = %s" % stit
        self.ax2.set_title(stit, ha="left", position=(0.35, 1), color=col)
        self.fig2.canvas.draw()

    def PassBandPress(self, event):
        if self.pband is None:
            return
        if not self.pband.contains(event)[0]:
            return
        tb = plt.get_current_fig_manager().toolbar
        # Inicializar self.dl para que cuando se clicka en el mismo sitio no haga nada
        self.dl = 0.0
        if tb.mode == "":
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
        if tb.mode == "" and self.pressevent is not None:
            self.pressevent = None
            self.gdl = self.gdl + self.dl
            self.set_fff(verb=False, dl=self.gdl)
            self.updatePassBand(remove=True)
            self.updatePatch()
            # self.ax2.set_xlim(self.ax2.get_xlim())
            if self.iclm:
                if self.icmp is not None:
                    self.zmode = self.icmp.zmode
                self.icmp = IntColorMap(self.p, zmode=self.zmode)
            self.fig.canvas.draw()
            self.fig2.canvas.draw()

    def ChangeSpaxelViewer(self, event):
        if event.key == "y":
            if self.K is None:
                print("*** Non Pycasso Object found! ***")
                return
            print("***** Pycasso Property *****")
            print("Input example: McorSD")
            lm = input("Enter property (Enter to abort): ")
            if len(lm) == 0:
                pass
            else:
                newprop = self.GetPycassoProp(lm.strip())
                if newprop is not None:
                    self.color = newprop
                    self.updateAx1(False)
                else:
                    print('Nothing to plot! Property "%s" not found' % lm.strip())

    def updateAx1(self, color=True):
        self.fig.clf()
        self.ax = self.fig.add_subplot(111)
        norm = rnorm(self.norm) if self.icmp.scale is None else rnorm(self.icmp.scale)
        if color:
            self.set_fff(verb=False, dl=self.gdl)
        self.p = self.ax.imshow(
            self.color,
            alpha=self.palpha,
            extent=self.ext,
            cmap=self.p.get_cmap(),
            norm=norm,
            interpolation="nearest",
            origin="lower",
            aspect="auto",
        )
        self.cmin, self.cmax = get_min_max(self.color)
        self.p.set_array(self.color)
        self.p.set_clim([self.cmin, self.cmax])
        if self.iclm:
            self.icmp = IntColorMap(self.p)
        if self.colorbar:
            self.fig.colorbar(self.p)
        self.ax.axis(self.ext)
        self.fig.canvas.draw()

    def GetPycassoProp(self, prop):
        lprop = [item for item in dir(self.K) if not item.startswith("_")]
        sprop = None
        if any([prop == item for item in lprop]):
            sprop = prop
        if LooseVersion(self.pversion) < LooseVersion("2.0.0"):
            if any([prop + "__yx" == item for item in lprop]):
                sprop = prop + "__yx"
        if sprop is None:
            return
        else:
            fprop = eval("self.K." + sprop)
            if LooseVersion(self.pversion) >= LooseVersion("2.0.0"):
                if self.K.hasSegmentationMask:
                    fprop = self.K.spatialize(fprop, extensive=False)
                if fprop.ndim == 3:
                    fprop = fprop.sum(axis=0)
            nz, ny, nx = self.dat.shape
            if (ny, nx) == fprop.shape:
                return fprop
            else:
                return None

    def GetPycassoRes(self, event):
        if event.key == "R":
            self.rshow = ~self.rshow
            if self.K is None:
                print("*** Non Pycasso Object found! ***")
                return
            else:
                if self.rshow:
                    print(">>> Residual Map")
                    self.dcolor = self.K.res
                    self.plotResidualMap()
                else:
                    print(">>> Data Map")
                    self.dcolor = self.dat
                    self.updateAx1()

    def plotResidualMap(self):
        self.fig4 = plt.figure(4, (9, 7))
        self.fig4.set_label("2D Residual Map")
        self.setWindowTitle(self.fig4, "2D Residual Map")
        self.ax4 = self.fig4.add_axes([0.09, 0.30, 0.81, 0.68])
        awlmin, awlmax = self.wl[0], self.wl[-1]
        pr = self.ax4.imshow(
            self.K.zres.T,
            interpolation="nearest",
            extent=[awlmin, awlmax, 0.5, self.K.zres.shape[1] + 0.5],
            cmap=plt.get_cmap("bwr"),
            aspect="auto",
            origin="lower",
        )
        ipr = IntColorMap(pr)
        self.ax4.set_ylabel("# Zone")
        self.ax5 = self.fig4.add_axes([0.09, 0.1, 0.81, 0.15])
        self.ax5.plot(self.wl, self.K.tres)
        self.ax5.set_xlabel("Wavelength ($\AA$)")
        self.ax5.set_ylim(-0.1, 0.1)
        self.ax5.set_xlim(awlmin, awlmax)
        self.ax6 = self.fig4.add_axes([0.92, 0.1, 0.03, 0.88])
        self.fig4.colorbar(pr, cax=self.ax6)
        self.fig4.canvas.draw()
        self.hline = None
        self.updateAx1()
        self.fig4.canvas.mpl_connect("button_press_event", self.ResidualViewer)
        plt.show()

    def ResidualViewer(self, event):
        tb = plt.get_current_fig_manager().toolbar
        if event.inaxes and tb.mode == "" and event.button == 1:
            x, y = event.xdata, int(round(event.ydata))
            if y > 0 and y <= self.zones.max() + 0.5:
                awlmin, awlmax = self.wl[0], self.wl[-1]
                self.ax5.cla()
                self.ax5.set_xlabel("Wavelength ($\AA$)")
                if self.hline is not None:
                    self.hline.remove()
                self.ax5.plot(self.wl, self.K.zres[:, y - 1])
                self.hline = self.ax4.axhline(y, lw=2, c="g", label="#%i" % y)
                self.ax5.set_ylim(-0.1, 0.1)
                self.ax5.set_xlim(awlmin, awlmax)
                self.ax5.text(
                    1.0,
                    -0.36,
                    "Zone #%i" % y,
                    weight="bold",
                    transform=self.ax5.transAxes,
                    ha="right",
                    multialignment="right",
                )
                self.fig4.canvas.draw()
        if not event.inaxes and tb.mode == "" and event.button == 1:
            if self.hline is not None:
                self.hline.remove()
                self.hline = None
            self.ax5.cla()
            self.ax5.set_xlabel("Wavelength ($\AA$)")
            self.ax5.plot(self.wl, self.K.tres)
            self.ax5.set_ylim(-0.1, 0.1)
            self.ax5.text(
                1.0,
                -0.36,
                "Integrated",
                weight="bold",
                transform=self.ax5.transAxes,
                ha="right",
                multialignment="right",
            )
            self.fig4.canvas.draw()

    def selectZone(self, event):
        if event.key == "n":
            if self.K is None:
                print("*** Non Pycasso Object found! ***")
                return
            else:
                print("***** Select Pycasso Zone *****")
                print("Input example: 3")
                zn = input("Enter property (Enter to abort): ")
                if len(zn) == 0:
                    pass
                else:
                    lzy, lzx = np.where(self.zones == int(zn))
                    if lzy.size < 1 or lzx.size < 1:
                        print("*** No available ZONE with ID: %i ***" % int(zn))
                        return
                    else:
                        for ix, iy in self.list:
                            idp = self.list.index((ix, iy))
                            self.pat[idp][0].remove()
                            self.list.remove((ix, iy))
                            self.pat.pop(idp)
                        for zy, zx in zip(lzy, lzx):
                            if (zx, zy) not in self.list:
                                ir, id = self.pix2res(zx, zy)
                                xr, yr = PRectangle(ir, id, self.sr)
                                rec = self.ax.fill(
                                    xr,
                                    yr,
                                    fill=self.cf,
                                    color=self.cc,
                                    ec=self.cc,
                                    alpha=self.ca,
                                    lw=self.clw,
                                )
                                self.list.append((zx, zy))
                                self.pat.append(rec)
                        self.fig.canvas.draw()

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
            ix, iy = self.list[0]
            self.getSpec(ix, iy)
            if self.fitspec_mode % 2 == 0:
                if plt.get_backend() != "TkAgg":
                    print(
                        '*** You need to set backend to "TkAgg" if you want to use PyRAF interactive fitting ***'
                    )
                else:
                    # sname = '%s_%s' % ('.'.join(self.name_fits.split('.')[0:-1]),self.idl.replace(',','_').replace(' ',''))
                    sname = "%s_%s" % (
                        ".".join(self.name_fits.split(".")[0:-1]),
                        "%s_%s" % (ix, iy),
                    )
                    tmpfits = tmpName(prefix="tmp_%s" % sname)
                    convert2iraf_spec(tmpfits, self.wl, self.spec, title=sname)
                    print(">>> Spectrum (%s) of %s" % (self.idl, self.name_fits))
                    splot(tmpfits)
                    if os.path.exists(tmpfits):
                        os.remove(tmpfits)
            else:
                sp = pyspeckit.Spectrum(
                    xarr=self.wl, data=self.spec, error=self.espec, header=self.hd
                )
                sp.plotter()
                sp.plotter.axis.set_xlabel(r"Wavelength $(\AA)$")
                sp.plotter.axis.set_ylabel(r"Flux $(\mathrm{erg/s/cm^2/\AA})$")
                sp.plotter.axis.set_title("%s (%s, %s)" % (sp.plotter.title, ix, iy))
                plt.show()

    def WindowManager(self):
        from matplotlib.widgets import CheckButtons, Slider, Button, RadioButtons
        from matplotlib import gridspec

        if not plt.fignum_exists(3):
            self.fig3 = plt.figure(3, self.winman_size)
            self.fig3.set_label(self.winman_label)
            self.setWindowTitle(self.fig3, self.winman_label)
            plt.get_current_fig_manager().toolbar.set_message = lambda x: None
            # Patch Spaxel Properties
            # No podemos anyadir 'axes' antes porque no se ve, y si lo ponemos despues, los widgets No funcionan
            ps = gridspec.GridSpec(
                3, 2, left=0.2, bottom=0.8, top=0.95, wspace=0.5, hspace=0.7, right=0.8
            )
            ax = plt.subplot(ps[:2, :], frameon=False, xticks=[], yticks=[])
            ax.text(0.5, 0.5, "Spaxel Properties", ha="center")
            # plw = Slider(plt.subplot(ps[1,:]), 'Linewidth', 0.0, 5.0, valinit=self.plw)
            pap = Slider(plt.subplot(ps[2, :]), "Alpha", 0.0, 1.0, valinit=self.palpha)
            # Circle Properties ------------------------
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
                # self.plw = plw.val; self.p.set_lw(self.plw)
                self.palpha = pap.val
                self.p.set_alpha(self.palpha)
                self.clw = clw.val
                self.ca = cap.val
                self.slw = slw.val
                self.sa = sap.val
                [(item[0].set_alpha(self.ca), item[0].set_lw(self.clw)) for item in self.pat]
                [(item[0].set_alpha(self.sa), item[0].set_lw(self.slw)) for item in self.cir]
                self.p.set_clim(self.p.get_clim())  # Hay que actualizar clim para que cambie alpha
                self.fig.canvas.draw()

            # plw.on_changed(update);
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
                    ctf.ax.set_facecolor(self.ccc)
                    self.clab = "Fill On" if self.cf else "Fill Off"
                    ctf.label.set_text(self.clab)
                    [item[0].set_fill(self.cf) for item in self.pat]
                    self.fig3.canvas.draw()
                    self.fig.canvas.draw()
                if event.inaxes == stf.ax.axes:
                    self.sf = not self.sf
                    self.scc = self.axcolor2 if self.sf else self.axcolor1
                    stf.color, stf.hovercolor = self.scc, self.scc
                    stf.ax.set_facecolor(self.scc)
                    self.slab = "Fill On" if self.sf else "Fill Off"
                    stf.label.set_text(self.slab)
                    [item[0].set_fill(self.sf) for item in self.cir]
                    self.fig3.canvas.draw()
                    self.fig.canvas.draw()
                if event.inaxes == sint.ax.axes:
                    self.integrated = not self.integrated
                    self.tcc = self.axcolor1 if not self.integrated else self.axcolor2
                    sint.color, sint.hovercolor = self.tcc, self.tcc
                    sint.ax.set_facecolor(self.tcc)
                    self.intlab = "Integrated On" if self.integrated else "Integrated Off"
                    sint.label.set_text(self.intlab)
                    self.fig3.canvas.draw()
                if event.inaxes == sind.ax.axes:
                    self.individual = not self.individual
                    self.icc = self.axcolor1 if not self.individual else self.axcolor2
                    sind.color, sind.hovercolor = self.icc, self.icc
                    sind.ax.set_facecolor(self.icc)
                    self.indlab = "Individual On" if self.individual else "Individual Off"
                    sind.label.set_text(self.indlab)
                    self.fig3.canvas.draw()
                if event.inaxes == stxt.ax.axes:
                    self.txt = not self.txt
                    self.xcc = self.axcolor1 if not self.txt else self.axcolor2
                    stxt.color, stxt.hovercolor = self.xcc, self.xcc
                    stxt.ax.set_facecolor(self.xcc)
                    self.txtlab = "Txt On" if self.txt else "Txt Off"
                    stxt.label.set_text(self.txtlab)
                    self.fig3.canvas.draw()
                if event.inaxes == sfit.ax.axes:
                    self.fits = not self.fits
                    self.fcc = self.axcolor1 if not self.fits else self.axcolor2
                    sfit.color, sfit.hovercolor = self.fcc, self.fcc
                    sfit.ax.set_facecolor(self.fcc)
                    self.fitlab = "Fits On" if self.fits else "Fits Off"
                    sfit.label.set_text(self.fitlab)
                    self.fig3.canvas.draw()

            self.fig3.canvas.mpl_connect("button_press_event", press)
            plt.show()

    def Sonification(self):
        if self.dsoni is None:
            return
        sys.path.append(self.dsoni)
        from .sonicube import SoniCube

        self.sc = SoniCube(
            self.fig,
            file=self.name_fits,
            data=self.dat,
            base_dir=self.dsoni,
            ref=(self.y_ref, self.x_ref),
        )
