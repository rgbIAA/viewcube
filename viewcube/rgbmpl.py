# Image module:
# 					ilinear(), isqrt(), ilog(), iasinh(),
# 					isfun(), rnorm(), zscale(),
# 					zsc_sample(), zsc_fit_line(),
# 					zsc_compute_sigma(), ColorMapList(),
# 					IntColorMap()
#
# System utilities:			lsfiles()
#
# ---------------------------------------------------------------------------
#
# Version 1.2.4
# Ruben Garcia Benito (RGB) - UAM/KIAA-PKU 2011
# PYTHON 3: 2024/10/10
#
# Scaling functions written by Min-Su Shin (ILOG IRAF opcion by RGB)
# Department of Astrophysical Sciences, Princeton University
# RNORM basada en APLpy
#
# Zscale function adapted from Numdisplay
# ---------------------------------------------------------------------------
#
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
import astropy.io.fits as pyfits
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy as np
import matplotlib
import fnmatch
import math
import sys
import os


MAX_REJECT = 0.5
MIN_NPIXELS = 5
GOOD_PIXEL = 0
BAD_PIXEL = 1
KREJ = 2.5
MAX_ITERATIONS = 5


##############################################################################
# ---------------------------- Subrutina ILINEAR -----------------------------
# ----------------------------------------------------------------------------
#
def ilinear(inputArray, scale_min=None, scale_max=None, verb=False):
    """Performs linear scaling of the input numpy array.

    @type inputArray: numpy array
    @param inputArray: image data array
    @type scale_min: float
    @param scale_min: minimum data value
    @type scale_max: float
    @param scale_max: maximum data value
    @verb = False --> Print type of scaling
    @rtype: numpy array
    @return: image data array

    """
    if verb is True:
        print("Scale : linear")

    imageData = np.array(inputArray, copy=True)

    if scale_min == None:
        scale_min = imageData.min()
    if scale_max == None:
        scale_max = imageData.max()

    imageData = imageData.clip(min=scale_min, max=scale_max)
    imageData = (imageData - scale_min) / (scale_max - scale_min)
    indices = np.where(imageData < 0)
    imageData[indices] = 0.0
    indices = np.where(imageData > 1)
    imageData[indices] = 1.0

    return imageData


# ----------------------------------------------------------------------------


##############################################################################
# ---------------------------- Subrutina ISQRT -------------------------------
# ----------------------------------------------------------------------------
#
def isqrt(inputArray, scale_min=None, scale_max=None, verb=False):
    """Performs sqrt scaling of the input numpy array.

    @type inputArray: numpy array
    @param inputArray: image data array
    @type scale_min: float
    @param scale_min: minimum data value
    @type scale_max: float
    @param scale_max: maximum data value
    @verb = False --> Print type of scaling
    @rtype: numpy array
    @return: image data array

    """

    if verb is True:
        print("Scale : sqrt")

    imageData = np.array(inputArray, copy=True)

    if scale_min == None:
        scale_min = imageData.min()
    if scale_max == None:
        scale_max = imageData.max()

    imageData = imageData.clip(min=scale_min, max=scale_max)
    imageData = imageData - scale_min
    indices = np.where(imageData < 0)
    imageData[indices] = 0.0
    imageData = np.sqrt(imageData)
    imageData = imageData / math.sqrt(scale_max - scale_min)

    return imageData


# ----------------------------------------------------------------------------


##############################################################################
# ---------------------------- Subrutina ILOG --------------------------------
# ----------------------------------------------------------------------------
#
def ilog(inputArray, scale_min=None, scale_max=None, oiraf=True, z1=1.0, z2=1000.0, verb=False):
    """Performs log10 scaling of the input numpy array.

    @type inputArray: numpy array
    @param inputArray: image data array
    @type scale_min: float
    @param scale_min: minimum data value
    @type scale_max: float
    @param scale_max: maximum data value
    @rtype: numpy array
    @return: image data array

    oiraf = True --> Opcion de usar el algoritmo de IRAF para
         escalar el array. Se mapean los valores del array
         desde 'z1' a 'z2' y se aplica el logaritmo
    """

    if verb is True:
        print("Scale : log")

    imageData = np.array(inputArray, copy=True)

    if scale_min == None:
        scale_min = imageData.min()
    if scale_max == None:
        scale_max = imageData.max()

    indices0 = np.where(imageData < scale_min)
    indices1 = np.where((imageData >= scale_min) & (imageData <= scale_max))
    indices2 = np.where(imageData > scale_max)

    if oiraf is False:
        factor = math.log10(scale_max - scale_min)
        imageData[indices0] = 0.0
        imageData[indices2] = 1.0
        try:
            imageData[indices1] = np.log10(imageData[indices1]) / factor
        except:
            print("Error on math.log10 for ", (imageData[i][j] - scale_min))

    else:
        a = (z2 - z1) / (scale_max - scale_min)
        imageData[imageData == 0.0] = np.nan
        imageData = (imageData - scale_min) * a + z1
        imageData[indices0] = z1
        imageData[indices2] = z2
        factor = np.log10(z2 - z1)
        imageData = np.log10(imageData) / factor

    imageData[np.isnan(imageData)] = 0.0
    imageData[np.isinf(imageData)] = 0.0

    return imageData


# ----------------------------------------------------------------------------


##############################################################################
# --------------------------- Subrutina IASINH -------------------------------
# ----------------------------------------------------------------------------
#
def iasinh(inputArray, scale_min=None, scale_max=None, non_linear=2.0, verb=False):
    """Performs asinh scaling of the input numpy array.

    @type inputArray: numpy array
    @param inputArray: image data array
    @type scale_min: float
    @param scale_min: minimum data value
    @type scale_max: float
    @param scale_max: maximum data value
    @type non_linear: float
    @param non_linear: non-linearity factor
    @rtype: numpy array
    @return: image data array

    """

    if verb is True:
        print("Scale : asinh")

    imageData = np.array(inputArray, copy=True)

    if scale_min == None:
        scale_min = imageData.min()
    if scale_max == None:
        scale_max = imageData.max()
    factor = np.arcsinh((scale_max - scale_min) / non_linear)
    indices0 = np.where(imageData < scale_min)
    indices1 = np.where((imageData >= scale_min) & (imageData <= scale_max))
    indices2 = np.where(imageData > scale_max)
    imageData[indices0] = 0.0
    imageData[indices2] = 1.0
    imageData[indices1] = np.arcsinh((imageData[indices1] - scale_min) / non_linear) / factor

    return imageData


# ----------------------------------------------------------------------------


##############################################################################
# --------------------------- Subrutina ISFUN --------------------------------
# ----------------------------------------------------------------------------
#
def isfun(
    input,
    scale="linear",
    vmin=None,
    vmax=None,
    zsc=True,
    oiraf=True,
    z1=1.0,
    z2=1000.0,
    nlf=2.0,
    omap=False,
    verb=False,
):
    """
    ##############################################################################
    # --------------------------- Subrutina ISFUN --------------------------------
    # ----------------------------------------------------------------------------
    #     Subrutina para aplicar un escalado a un array
    # ----------------------------------------------------------------------------
    #     USO:
    #		rsqrt = isfun(r,'sqrt',1.,100.)
    #
    #	input --> String con nombre del fichero FITS o array numpy
    #	scale = 'linear' --> A elegir entre:
    #		'linear' | 'sqrt' | 'log' | 'asinh'
    #	vmin = None --> Valor por debajo del cual se enmascara
    #	vmax = None --> Valor por encima del cual se enmascara
    #	zsc = True --> Si 'vmin' y 'vmax' son 'None' y 'zsc = True',
    #		se aplica el algoritmo 'zscale' de 'numdisplay' (~ IRAF)
    #	oiraf = True --> Opcion de usar el algoritmo de IRAF para
    #		escalar el array. Se mapean los valores del array
    #		desde 'z1' a 'z2' y se aplica el logaritmo
    #	z1 = 1. --> Valor minimo para mapear el array en 'log' segun
    #		el algoritmo de IRAF
    #	z2 = 1000. --> Valor maximo para mapear el array en 'log' segun
    #		el algoritmo de IRAF
    #	nlf = 2.0 --> Nonlinear factor para la opcion scale = 'asinh'
    #	omap = False --> Opcion para devolver un objeto 'maputils.FITSimage'
    #	verb = False --> Imprime por pantalla la opcion de 'scale'
    # ----------------------------------------------------------------------------
    """
    from kapteyn import maputils

    lops = ["linear", "sqrt", "log", "asinh"]

    if scale not in lops:
        print('Opcion "scale = ' + scale + '" NO encontrada')
        print(lops)
        sys.exit()

    if type(input) is str:
        mstr = maputils.FITSimage(input)
        iarray, ihdr = mstr.dat, mstr.hdr
    else:
        iarray = input

    if vmin is None and vmax is None and zsc is True:
        vmin, vmax = zscale(iarray)

    if scale == "linear":
        dat = ilinear(iarray, scale_min=vmin, scale_max=vmax, verb=verb)

    elif scale == "sqrt":
        dat = isqrt(iarray, scale_min=vmin, scale_max=vmax, verb=verb)

    elif scale == "log":
        dat = ilog(iarray, scale_min=vmin, scale_max=vmax, oiraf=oiraf, z1=z1, z2=z2, verb=verb)

    else:
        dat = iasinh(iarray, scale_min=vmin, scale_max=vmax, non_linear=nlf, verb=verb)

    if omap is True:
        return maputils.FITSimage(externaldata=dat, externalheader=ihdr)
    else:
        return dat


# ----------------------------------------------------------------------------


##############################################################################
# --------------------------- Subrutina RNORM --------------------------------
#    A Normalize class for imshow that allows different scaling functions
#    for astronomical images.
# ----------------------------------------------------------------------------
#
class rnorm(Normalize):
    """
    ##############################################################################
    # --------------------------- Subrutina RNORM --------------------------------
    #    A Normalize class for imshow that allows different scaling functions
    #    for astronomical images.
    # ----------------------------------------------------------------------------
    #     USO:
    #		norm = rnorm('ilog'); imshow(z,norm=norm)
    #
    #	scale = 'linear' --> A elegir entre:
    #
    #		'linear' | 'sqrt' | 'log' | 'ilog' | 'asinh' | 'power'
    #
    #		'ilog' es el algoritmo usado en IRAF
    #	vmin = None --> Valor por debajo del cual se enmascara
    #	vmax = None --> Valor por encima del cual se enmascara
    #	zsc = None --> Array Numpy 2-D de la imagen para calcular max y min
    #		aplicando el algoritmo 'zscale' de 'numdisplay' (~ IRAF)
    #	z1 = 1. --> Valor minimo para mapear el array en 'ilog' segun
    #		el algoritmo de IRAF
    #	z2 = 1000. --> Valor maximo para mapear el array en 'ilog' segun
    #		el algoritmo de IRAF
    #	exponent = 2. --> Exponente de la opcion 'power'
    #	vmid = None --> Valor del pixel medio para 'log' y 'asinh'
    #		(0.05 y -0.033 por defecto)
    #	clip = False --> Si 'clip = True' y el valor dado cae fuera del
    #		rango, se convierte a 0 o 1, el mas cercano
    # ----------------------------------------------------------------------------
    """

    def __init__(
        self,
        scale="linear",
        vmin=None,
        vmax=None,
        zsc=None,
        z1=1.0,
        z2=1000.0,
        exponent=2.0,
        vmid=None,
        clip=False,
    ):
        # Call original initalization routine
        Normalize.__init__(self, vmin=vmin, vmax=vmax, clip=clip)

        # Save parameters
        self.z1, self.z2 = z1, z2
        self.exponent = exponent
        self.scale = scale
        self.zsc = zsc

        if np.equal(vmid, None):
            if scale == "log":
                self.midpoint = 0.05
            elif scale == "asinh":
                self.midpoint = -0.033
            else:
                self.midpoint = None
        else:
            self.midpoint = vmid

    def __call__(self, value, clip=None):
        vmin, vmax = self.vmin, self.vmax
        exponent = self.exponent
        z1, z2 = self.z1, self.z2
        midpoint = self.midpoint
        zsc = self.zsc

        if zsc is not None:
            vmin, vmax = zscale(zsc)

        # ORIGINAL MATPLOTLIB CODE

        if clip is None:
            clip = self.clip

        val = ma.asarray(value)

        self.autoscale_None(val)

        if vmin > vmax:
            raise ValueError("minvalue must be less than or equal to maxvalue")
        elif vmin == vmax:
            return 0.0 * val
        else:
            if clip:
                mask = ma.getmask(val)
                val = ma.array(np.clip(val.filled(vmax), vmin, vmax), mask=mask)
            result = (val - vmin) * (1.0 / (vmax - vmin))

            # RUMPL

            if self.scale == "linear":
                pass

            elif self.scale == "log":
                result = ma.log10((result / self.midpoint) + 1.0) / ma.log10(
                    (1.0 / self.midpoint) + 1.0
                )

            elif self.scale == "ilog":
                a = (z2 - z1) / (vmax - vmin)
                val = (val - vmin) * a + z1
                factor = ma.log10(z2 - z1)
                result = ma.log10(val) / factor
                result[np.isinf(result)] = 0.0

            elif self.scale == "sqrt":
                result = ma.sqrt(result)

            elif self.scale == "asinh":
                result = ma.arcsinh(result / self.midpoint) / ma.arcsinh(1.0 / self.midpoint)

            elif self.scale == "power":
                result = ma.power(result, exponent)

            else:
                raise Exception("Unknown scale in 'rnorm': %s" % self.scale)

        return result

    def inverse(self, value):
        # ORIGINAL MATPLOTLIB CODE
        if not self.scaled():
            raise ValueError("Not invertible until scaled")

        vmin, vmax = self.vmin, self.vmax
        z1, z2 = self.z1, self.z2
        zsc = self.zsc

        if zsc is not None:
            vmin, vmax = zscale(zsc)

        # CUSTOM APLPY CODE
        val = ma.asarray(value)

        if self.scale == "linear":
            pass

        elif self.scale == "log":
            val = self.midpoint * (
                ma.power(10.0, (val * ma.log10(1.0 / self.midpoint + 1.0))) - 1.0
            )

        elif self.scale == "ilog":
            a = (z2 - z1) / (vmax - vmin)
            factor = ma.log10(z2 - z1)
            val = vmin + (ma.power(10.0, (val * factor)) - z1) / a

        elif self.scale == "sqrt":
            val = val * val

        elif self.scale == "asinh":
            val = self.midpoint * ma.sinh(val * ma.arcsinh(1.0 / self.midpoint))

        elif self.scale == "power":
            val = ma.power(val, (1.0 / self.exponent))

        else:
            raise Exception("Unknown scale in 'rnorm': %s" % self.scale)

        if self.scale == "ilog":
            return val
        else:
            return vmin + val * (vmax - vmin)


# ----------------------------------------------------------------------------


##############################################################################
# --------------------------- Subrutina ZSCALE -------------------------------
# ----------------------------------------------------------------------------
#
def zscale(image, nsamples=1000, contrast=0.25, bpmask=None, zmask=None):
    """Implement IRAF zscale algorithm nsamples=1000 and contrast=0.25 are the
    IRAF display task defaults bpmask and zmask not implemented yet image is a
    2-d numpy array (from Numdisplay)
    returns (z1, z2)
    """

    # Sample the image
    samples = zsc_sample(image, nsamples, bpmask, zmask)
    npix = len(samples)
    samples.sort()
    # zmin = samples[0]; zmax = samples[-1]
    # Avoid getting Nan values as min/max
    # zmin = nanmin(samples)
    # zmax = nanmax(samples)
    zmin = ma.min(ma.masked_invalid(samples))
    zmax = ma.max(ma.masked_invalid(samples))
    # For a zero-indexed array
    center_pixel = (npix - 1) // 2
    if npix % 2 == 1:
        median = samples[center_pixel]
    else:
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1])
    #
    # Fit a line to the sorted array of samples
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    ngrow = max(1, int(npix * 0.01))
    ngoodpix, zstart, zslope = zsc_fit_line(samples, npix, KREJ, ngrow, MAX_ITERATIONS)

    if ngoodpix < minpix:
        z1 = zmin
        z2 = zmax
    else:
        if contrast > 0:
            zslope = zslope / contrast
        z1 = max(zmin, median - (center_pixel - 1) * zslope)
        z2 = min(zmax, median + (npix - center_pixel) * zslope)

    return sorted((z1, z2))


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------- Subrutina ZSC_SAMPLE -----------------------------
# ----------------------------------------------------------------------------
#
def zsc_sample(image, maxpix, bpmask=None, zmask=None):
    """Adapted from Numdisplay"""
    # Figure out which pixels to use for the zscale algorithm
    # Returns the 1-d array samples
    # Don't worry about the bad pixel mask or zmask for the moment
    # Sample in a square grid, and return the first maxpix in the sample
    lis = len(image.shape)
    nc = image.shape[0]
    if lis < 2:
        nl = 2
    else:
        nl = image.shape[1]
    stride = max(1.0, math.sqrt((nc - 1) * (nl - 1) / float(maxpix)))
    stride = int(stride)
    if lis < 2:
        samples = image[::stride].flatten()
    else:
        samples = image[::stride, ::stride].flatten()
    return samples[:maxpix]


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------- Subrutina ZSC_FIT_LINE ---------------------------
# ----------------------------------------------------------------------------
#
def zsc_fit_line(samples, npix, krej, ngrow, maxiter):
    """(from Numdisplay)"""
    #
    # First re-map indices from -1.0 to 1.0
    xscale = 2.0 / (npix - 1)
    xnorm = np.arange(npix)
    xnorm = xnorm * xscale - 1.0

    ngoodpix = npix
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    last_ngoodpix = npix + 1

    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad
    badpix = np.zeros(npix, dtype="int32")

    #
    #  Iterate

    for niter in range(maxiter):
        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix):
            break

        # Accumulate sums to calculate straight line fit
        goodpixels = np.where(badpix == GOOD_PIXEL)
        sumx = xnorm[goodpixels].sum()
        sumxx = (xnorm[goodpixels] * xnorm[goodpixels]).sum()
        sumxy = (xnorm[goodpixels] * samples[goodpixels]).sum()
        sumy = samples[goodpixels].sum()
        sum = len(goodpixels[0])

        delta = sum * sumxx - sumx * sumx
        # Slope and intercept
        intercept = (sumxx * sumy - sumx * sumxy) / delta
        slope = (sum * sumxy - sumx * sumy) / delta

        # Subtract fitted line from the data array
        fitted = xnorm * slope + intercept
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        ngoodpix, mean, sigma = zsc_compute_sigma(flat, badpix, npix)

        threshold = sigma * krej

        # Detect and reject pixels further than k*sigma from the fitted line
        lcut = -threshold
        hcut = threshold
        below = np.where(flat < lcut)
        above = np.where(flat > hcut)

        badpix[below] = BAD_PIXEL
        badpix[above] = BAD_PIXEL

        # Convolve with a kernel of length ngrow
        kernel = np.ones(ngrow, dtype="int32")
        badpix = np.convolve(badpix, kernel, mode="same")

        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0])

        niter += 1

    # Transform the line coefficients back to the X range [0:npix-1]
    zstart = intercept - slope
    zslope = slope * xscale

    return ngoodpix, zstart, zslope


# ----------------------------------------------------------------------------


##############################################################################
# ----------------------- Subrutina ZSC_COMPUTE_SIGMA ------------------------
# ----------------------------------------------------------------------------
#
def zsc_compute_sigma(flat, badpix, npix):
    """(from Numdisplay)"""
    # Compute the rms deviation from the mean of a flattened array.
    # Ignore rejected pixels

    # Accumulate sum and sum of squares
    goodpixels = np.where(badpix == GOOD_PIXEL)
    sumz = flat[goodpixels].sum()
    sumsq = (flat[goodpixels] * flat[goodpixels]).sum()
    ngoodpix = len(goodpixels[0])
    if ngoodpix == 0:
        mean = None
        sigma = None
    elif ngoodpix == 1:
        mean = sumz
        sigma = None
    else:
        mean = sumz / ngoodpix
        temp = sumsq / (ngoodpix - 1) - sumz * sumz / (ngoodpix * (ngoodpix - 1))
        if temp < 0:
            sigma = 0.0
        else:
            sigma = math.sqrt(temp)

    return ngoodpix, mean, sigma


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------ Subrutina COLORMAPLIST ----------------------------
# ----------------------------------------------------------------------------
#
class ColorMapList(object):
    """
    List with names of colormaps as used in combination
    with keyword parameter *cmap* in the constructor

    Adapted from Kapteyn Package
    """

    def __init__(self):
        from pylab import cm

        # A list with available Matplotlib color maps
        # The '_r' entries are reversed versions
        self.colormaps = sorted([m for m in cm.datad.keys() if not m.endswith("_r")])

    def add(self, clist):
        if not issequence(clist):
            clist = [clist]
        for c in clist[::-1]:
            self.colormaps.insert(0, c)


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------- Subrutina INTCOLORMAP ----------------------------
# ----------------------------------------------------------------------------
#
class IntColorMap:
    def __init__(
        self,
        im=None,
        data=None,
        z1=1.0,
        z2=1000.0,
        exponent=2.0,
        ccms=True,
        cdrcm=True,
        opf=False,
        zmode=False,
    ):
        """Class to get Interactive Colormap and Scaling settings in Matplotlib.

             im = None --> The class gets the actual image [gci()]. "im" could be a
                also a Patch Collection
             data = None --> Data to get the maximum and minimum values to set the
                colormap settings. The class gets de values from "im"
             z1 = 1. --> Valor minimo para mapear el array en 'ilog' segun
                el algoritmo de IRAF [rnorm()]
             z2 = 1000. --> Valor maximo para mapear el array en 'ilog' segun
                el algoritmo de IRAF [rnorm()]
             exponent = 2. --> Exponente de la opcion 'power' de rnorm()
             ccms = True --> Opcion para habilitar "ChangeColorMapScale". Se cambia el mapa
                de colores.
             cdrcm = True --> Opcion para habilitar "ChangeDynamicRangeColorMap". Cambia el
                rango dinamico y la escala.
             opf = False --> Forzar el rango dinamico al valor minimo y maximo de los
                datos. Por defecto, se pueden superar estos valores con "cdrcm"
             zmode = False -> Zscale mode

        # ----------------------------------------------------------------------------
             Active keys:
                * Right click hover mouse: Move the mouse over the axes of the figure while
                  holding down the right mouse button
                * 'z' --> Press 'z' to use the "Zscale" scaling
                * '+/-' --> Press '+' or '-' to change the ColorMap
                * 'm' --> Press 'm' to print the available ColorMaps and select one
                * 'i' --> Press 'i' to invert the ColorMap
                * '0' --> Press '0' to print the available scalings
                * 'v' --> Press 'v' to choose the Min/Max values of the Scaling Map
        # ----------------------------------------------------------------------------
             USE:
                        int = IntColorMap()

                DO NOT call it simply as "IntColorMap()", since the interactive
                options will not work.
        # ----------------------------------------------------------------------------
        """
        from pylab import gci  # connect, gcf

        self.data = data
        self.im = im
        self.ccms = ccms
        self.cdrcm = cdrcm
        self.opf = opf
        self.zmode = zmode

        if self.im is None:
            self.im = gci()

        if self.data is not None:
            self.cmin, self.cmax = self.data.min(), self.data.max()
        else:
            self.data = self.im.get_array()
            self.cmin, self.cmax = self.data.min(), self.data.max()
        self.crange = abs(self.cmax - self.cmin)
        self.rmin, self.rmax = self.cmin, self.cmax

        self.scales = {"1": "linear", "2": "ilog", "3": "sqrt", "4": "power", "5": "asinh"}
        self.scale = None
        self.z1 = z1
        self.z2 = z2
        self.exponent = exponent
        self.cmindx = 0

        if self.zmode:
            self.ColorMapZscale()

        if self.ccms is True:
            self.im.axes.figure.canvas.mpl_connect("key_press_event", self.ChangeColorMapScale)
        if self.cdrcm is True:
            self.im.axes.figure.canvas.mpl_connect(
                "motion_notify_event", self.ChangeDynamicRangeColorMap
            )
            self.im.axes.figure.canvas.mpl_connect("key_press_event", self.setColorMapZscale)
            self.im.axes.figure.canvas.mpl_connect("button_release_event", self.DynamicRangeInfo)

    def ChangeDynamicRangeColorMap(self, event):
        tb = plt.get_current_fig_manager().toolbar
        if event.button == 3 and event.inaxes and tb.mode == "":
            # event.canvas.figure.clear()
            # ax = event.canvas.figure.gca()
            x, y = event.xdata, event.ydata
            xy = self.im.axes.transData.transform((x, y))
            x, y = self.im.axes.transAxes.inverted().transform(xy)
            # x, y = fig.transFigure.inverted().transform((event.x, event.y))
            scale = self.crange * (y + 0.0001)  # /2.0
            offset = self.rmin + x * self.crange

            self.vmin = offset - scale
            self.vmax = offset + scale

            if self.opf is True:
                if self.vmin < self.rmin:
                    self.vmin = self.rmin
                if self.vmax > self.rmax:
                    self.vmax = self.rmax

            self.im.set_clim(self.vmin, self.vmax)
            self.im.axes.figure.canvas.draw()
            # gcf().canvas.draw()

    def DynamicRangeInfo(self, event):
        tb = plt.get_current_fig_manager().toolbar
        if event.button == 3 and event.inaxes and tb.mode == "":
            print(
                "Min/Max = %8.4e / %8.4e  ||  Set to %8.4e / %8.4e"
                % (self.cmin, self.cmax, self.vmin, self.vmax)
            )

    def setColorMapZscale(self, event):
        omode = self.zmode
        if event.key == "z":
            self.zmode = not self.zmode
        if self.zmode and self.zmode != omode:
            self.ColorMapZscale()
        if not self.zmode and self.zmode != omode:
            self.ColorMapMinMax()

    def ColorMapZscale(self):
        self.rmin, self.rmax = zscale(self.data)
        self.crange = abs(self.rmax - self.rmin)
        self.im.set_clim(self.rmin, self.rmax)
        print(
            "Min/Max = %8.4e / %8.4e  ||  ZSCALE set to %8.4e / %8.4e"
            % (self.cmin, self.cmax, self.rmin, self.rmax)
        )
        self.im.axes.figure.canvas.draw()

    def ColorMapMinMax(self):
        self.rmin, self.rmax = self.cmin, self.cmax
        self.crange = abs(self.cmax - self.cmin)
        self.im.set_clim(self.rmin, self.rmax)
        print(
            "Min/Max = %8.4e / %8.4e  ||  Dynamical range set to MIN/MAX" % (self.cmin, self.cmax)
        )
        self.im.axes.figure.canvas.draw()

    def ChangeColorMapScale(self, event):
        colormaps = ColorMapList().colormaps
        # This can happen when the caps lock is on.
        if event.key is None:
            return
        # Request for another color map with page up/down keys
        if event.key in ["pageup", "pagedown", "+", "-"]:
            lm = len(colormaps)
            if event.key in ["pageup", "+"]:
                self.cmindx += 1
                if self.cmindx >= lm:
                    self.cmindx = 0
            if event.key in ["pagedown", "-"]:
                self.cmindx -= 1
                if self.cmindx < 0:
                    self.cmindx = lm - 1
            newmap = colormaps[self.cmindx]
            self.im.set_cmap(newmap)
            self.im.axes.figure.canvas.draw()
            mes = "Color map set to '%s'" % newmap
            print(mes)

        elif event.key == "m":
            print("Available Color Maps:\n")
            print(colormaps)
            newmap = input("\nWrite Colormap Name (Enter to exit): ")
            if len(newmap) == 0:
                pass
            elif newmap not in colormaps:
                print('Colormap "' + newmap + '" does not exists')
            else:
                self.im.set_cmap(newmap)
                self.im.axes.figure.canvas.draw()
                mes = "Color map set to '%s'" % newmap
                print(mes)

        elif event.key == "i":
            map = self.im.cmap.name
            if map.endswith("_r"):
                imap = map.replace("_r", "")
            else:
                imap = map + "_r"
            self.im.set_cmap(imap)
            self.im.axes.figure.canvas.draw()

        # Request for another scale, linear, logarithmic etc.
        elif event.key in self.scales:
            self.scale = self.scales[event.key]
            nclim = self.im.get_clim()
            self.im.set_norm(rnorm(self.scale, z1=self.z1, z2=self.z2, exponent=self.exponent))
            self.im.set_clim(nclim)
            self.im.axes.figure.canvas.draw()
            mes = "Color Norm/Scale set to '%s'" % self.scale
            print(mes)

        elif event.key == "0":
            print("Press the associated number to change the Scale. Available Scales: ")
            for key in sorted(self.scales.keys()):
                print("%s: %s") % (key, self.scales[key])

        elif event.key == "v":
            print('Input new Min/Max [tuple: (min,max)] values ("None" also possible)')
            nval = input("[Enter for nothing to do]: ")
            if len(nval) == 0:
                print("Actual values: %8.4e / %8.4e" % self.im.get_clim())
                # self.im.set_clim([self.cmin, self.cmax])
                # self.im.axes.figure.canvas.draw()
            else:
                try:
                    if type(eval(nval)) in (float, int):
                        print("*** You should write 2-item tuple (min,max) ***")
                    elif len(eval(nval)) != 2:
                        print("*** You should write 2-item tuple (min,max) ***")
                    else:
                        self.im.set_clim(eval(nval))
                        self.im.axes.figure.canvas.draw()
                        smin, smax = self.im.get_clim()
                        print(
                            "Min/Max = %8.4e / %8.4e  ||  Set to %8.4e / %8.4e"
                            % (self.cmin, self.cmax, smin, smax)
                        )
                except (NameError, TypeError, SyntaxError):
                    print("*** You should write 2-item tuple (min,max) ***")


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------------ Subrutina LSFILES ---------------------------
# ----------------------------------------------------------------------------
def lsfiles(strfind, dir=".", path=False):
    if (path is False) or (dir == "."):
        return sorted(fnmatch.filter(os.listdir(dir), strfind))
    if (path is not False) and (dir != "."):
        return [dir + item for item in sorted(fnmatch.filter(os.listdir(dir), strfind))]


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------------ Subrutina CKFILES ---------------------------
# ----------------------------------------------------------------------------
def ckfiles(file):
    if type(file) is str:
        delf = [file]

    if type(file) is list:
        delf = file

    for item in delf:
        if os.path.exists(item) == False:
            print("**********************************************************************")
            print('********** EL ARCHIVO "' + item + '" NO EXISTE **********')
            print("**********************************************************************")
            sys.exit()


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------------ Subrutina OFITS -----------------------------
# ----------------------------------------------------------------------------
def ofits(name, ax=0, hd=True):
    hdu = pyfits.open(name)
    if hd is True:
        data, header = hdu[ax].data, hdu[ax].header
        hdu.close()
        return data, header
    else:
        data = hdu[ax].data
        hdu.close()
        return data


# ----------------------------------------------------------------------------


##############################################################################
# ----------------------------- Subrutina SAVESPEC ---------------------------
# ----------------------------------------------------------------------------
def SaveSpec(
    wl,
    spec,
    name,
    fits=False,
    txt=True,
    infotxt=None,
    hd=None,
    infohd=None,
    endtxt=".txt",
    fmt="%-10.2f %14.6e",
):
    if txt:
        f = open(name + endtxt, "w")
        if infotxt is not None:
            if type(infotxt) is not list:
                infotxt = [infotxt]
            for item in infotxt:
                f.write("# " + item + "\n")
        for iwl, ispec in zip(wl, spec):
            f.write((fmt + "\n") % (iwl, ispec))
        f.close()
    if fits:
        nhd = hd.copy()
        if infohd is not None:
            for key, value in infohd:
                nhd.update(key, value)
        pyfits.writeto(name + ".fits", spec, nhd, clobber=True)


# ----------------------------------------------------------------------------


##############################################################################
# -------------------------- Subrutina GETIDFILTER ---------------------------
# ----------------------------------------------------------------------------
def GetIdFilter(list, filter, dfil="."):
    # if any(filter in item for item in list):
    lfil = lsfiles("*" + filter + "*", dfil)
    if len(lfil) == 0:
        print('"' + filter + '" NOT found. Set to "' + list[0] + '"')
        return 0
    else:
        print("Filter: " + ".".join(lfil[0].split(".")[:-1]))
        return list.index(lfil[0])


# ----------------------------------------------------------------------------


##############################################################################
# -------------------------- Subrutina GETIDFILTER ---------------------------
# ----------------------------------------------------------------------------
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


# ----------------------------------------------------------------------------


##############################################################################
# ------------------------ Subrutina GETLAMBDALIMITS -------------------------
# ----------------------------------------------------------------------------
def GetLambdaLimits(wl, pt=0.05, wlim=None):
    min = wl.min()
    max = wl.max()
    if wlim is not None:
        if type(wlim) not in [list, tuple] or len(wlim) != 2:
            print("Wavelength limits should be a tuple or list of two items: ex. --> (None, 6200)")
        else:
            wlimmin, wlimmax = wlim
            if wlimmin is not None:
                min = wlimmin
            if wlimmax is not None:
                max = wlimmax
    range = abs(max - min)
    return min - range * pt, max + range * pt


# ----------------------------------------------------------------------------


##############################################################################
# ---------------------------- Subrutina INTFILTER ---------------------------
# ----------------------------------------------------------------------------
def IntFilter(ifi, fil, lm, dat, pasb=False, verb=True, dfil="", dl=None):
    lf, ff = np.loadtxt(dfil + fil[ifi], unpack=True)
    if dl is not None:
        lf = lf + dl
    if verb is True:
        print("Selected Filter: " + fil[ifi])
    fff = np.interp(lm, lf, ff)
    if pasb is True:
        return np.trapz(lm * fff * dat, lm, axis=1) / np.trapz(fff * lm, lm), fff
    else:
        return np.trapz(lm * fff * dat, lm, axis=1) / np.trapz(fff * lm, lm)


# ----------------------------------------------------------------------------
