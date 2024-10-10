############################################################################
#                              VIEWCUBE UTILS                              #
#                                 PYTHON 3                                 #
#                                                                          #
# RGB@IAA ---> Last Change: 2024/10/10                                     #
############################################################################
#
#
#
################################ VERSION ###################################
VERSION = "0.0.3"                                                          #
############################################################################
#
from astropy.io import fits as pyfits
import numpy as np
import fnmatch
import sys
import os


# -------------------------------------------------------------------------------------
def check_list(var, return_None=True, include_array=False):
    if not (return_None and var is None):
        instances = (list, tuple, np.ndarray) if include_array else (list, tuple)
        if not isinstance(var, instances):
            var = [var]
    return var


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def lsfiles(strfind, dir=".", path=False):
    if os.path.exists(dir):
        if not path or (dir == "."):
            return sorted(fnmatch.filter(os.listdir(dir), strfind))
        if path and (dir != "."):
            return [dir + item for item in sorted(fnmatch.filter(os.listdir(dir), strfind))]
    else:
        return None


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
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


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def ofits(name, ax=0, hd=True, clip=True, valclip=1e100):
    hdu = pyfits.open(name)
    if hd is True:
        data, header = hdu[ax].data, hdu[ax].header
        hdu.close()
        if clip is True:
            data[data < -valclip] = numpy.inf
            data[data > valclip] = numpy.inf
        return data, header
    else:
        data = hdu[ax].data
        hdu.close()
        if clip is True:
            data[data < -valclip] = numpy.inf
            data[data > valclip] = numpy.inf
        return data


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def convert2iraf_spec(name, wave, flux, title=None):
    crval = wave[0]
    cdelt = (wave[-1] - crval) / wave.size
    nwave = crval + cdelt * np.arange(wave.size)
    nflux = np.interp(nwave, wave, flux)
    hdu = pyfits.PrimaryHDU(nflux)
    hdu.header["CRVAL1"] = crval
    hdu.header["CDELT1"] = cdelt
    if title is not None:
        hdu.header["OBJECT"] = title
    hdu.writeto(name, clobber=True)


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def save_spec(
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


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def get_min_max(nar, mask=np.nan):
    if nar[np.isfinite(nar)].size > 0:
        return nar[np.isfinite(nar)].min(), nar[np.isfinite(nar)].max()
    else:
        return mask, mask
    # return numpy.nanmin(nar[~isinf(nar)]), numpy.nanmax(nar[~isinf(nar)])


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def image_quadrant(image, xc=None, yc=None, chunks=3.0, indices=False):
    if isinstance(image, np.ndarray):
        ny, nx = image.shape
    else:
        ny, nx = image
    if isinstance(chunks, (list, tuple)):
        fx, fy = chunks
    else:
        fx, fy = chunks, chunks

    if xc is None:
        dx = np.round(nx / fx, 0)
        if dx < 1.0:
            dx = 1.0
        xc = int(np.round(nx / 2.0, 0) - dx / 2.0)
        xc = slice(xc, int(xc + dx))
    if yc is None:
        dy = np.round(ny / fy, 0)
        if dy < 1.0:
            dy = 1.0
        yc = int(np.round(ny / 2.0, 0) - dy / 2.0)
        yc = slice(yc, int(yc + dy))

    mask = np.zeros((ny, nx), dtype=bool)
    mask[yc, xc] = True

    if indices:
        return yc, xc
    else:
        return mask


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
def image_max_pixel(image, integer=True, **kwargs):
    y, x = image_quadrant(image, indices=True, **kwargs)
    yc, xc = np.unravel_index(np.nanargmax(image[y, x], axis=None), np.shape(image[y, x]))
    yc += y.start
    xc += x.start
    return yc, xc


# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
class LoadFits:
    def __init__(
        self,
        name,
        exdata=None,
        exflag=None,
        exerror=None,
        exerw=None,
        exfibc=None,
        exflat=None,
        exwave=None,
        exsensf=None,
        exhdr=0,
        error_file=None,
        flag_file=None,
        wave_file=None,
        sensf_file=None,
        hdisp="DISPAXIS",
        hcrval="CRVAL",
        hcdelt="CDELT",
        crval=None,
        cdelt=None,
        haxis=True,
        specaxis=None,
        hvel=None,
        z=None,
        velocity=None,
        guess=True,
        flip2disp=True,
        rss=False,
        f2e=None,
        bflag=None,
        c=299792.458,
        dist_pc=None,
        filters=None,
        filter_props=None,
        filter_mode=None,
        sensf=True,
        units=None,
        ivar=False,
        verbo=True,
    ):
        """
        Reads a FITS file and returns an object with some attributes

        Parameters
        ----------
        name : string
              String with the name of the FITS file
        exdata : None or int, optional
              Number of the HDU to read the data. By default, 0
        exflag : None or int, optional
              Number of the HDU to read the mask/flag data. If it is CALIFA data, it will search
              for 'BADPIX' automatically if "getcalifa = True"
        exerror : None or int, optional
              Number of the HDU to read the error data. If it is CALIFA data, it will search
              for 'ERROR' automatically if "getcalifa = True"
        exerw : None or int, optional
              Number of the HDU to read the error weighted data. If it is CALIFA data, it will search
              for 'ERRWEIGHT' automatically if "getcalifa = True"
        exfibc : None or int, optional
              Number of the HDU to read the fiber cover data. If it is CALIFA data, it will search
              for 'FIBCOVER' automatically if "getcalifa = True"
        exflat : None or int, optional
              Number of the HDU to read the flat data. If it is CALIFA data, it will search
              for 'FLAT' automatically if "getcalifa = True"
        exwave : string of integer
              Number or name (string) of the HDU where to read the wavelength array
        exsensf : string of integer
              Number or name (string) of the HDU where to read the SENSFUNC
        exhdr : int
              Number of the HDU where to read the Header
        error_file : string
              External FITS file to read the ERRORs (see function ReadError())
        flag_file : string
              External FITS file to read the FLAGs (see function ReadFlag())
        wave_file : string
              External FITS file to read the WAVELENGTH solution (see function ReadWave())
        sensf_file : string
              External FITS file to read the SENSFUNC solution (see function ReadSensf())
        hdisp : string
              String for reading the dispersion axis axis ('DISPAXIS')
        hcrval : string
              String for reading the 'CRVAL' information
        hdelt : string
              String for reading the 'CDELT' information
        crval : float
              Value for 'CRVAL' information. It overwrites "hcrval"
        haxis : Bool
              Uses "specaxis" to add the dimension integer to the header "hcrval" and "hcdelt"
              keywords. If "hcrval" and "hcdelt" are given and are absolute (nothing to add),
              set "haxis = False"
        delt : float
              Value for 'CDELT' information. It overwrites "cdelt"
        specaxis : None or int, optional
              Number to read the wavelenght (CRVAL/CDELT) information. If 'DISPAXIS'
              exists, it will be used. If does not exists, by default it will search
              the 1 axis (i.e.: CRVAL1 --> CRVAL+specaxis). If it fails (because there
              is not wavelength information or the provided numbers are wrong), it
              will give "None" results for the wavelength attributes
        hvel : string
              String for reading the recessional velocity (km/s) from the header
        z : None or float, optional
              If the redshift of the object is provided (or is readed from the header), a
              "wave_rest" attribute will be written with the rest wavelength for that redshift
              (see RestWave())
        velocity : None or float, optinal
              Velocity in km/s to obtain the "wave_rest" attribute (used instead of z if both provided)
        guess : Bool
              For CALIFA data, it will search for the 'BADPIX', 'ERROR' and 'ERRWEIGHT' HDUs.
              Also, it will calculate the "wave_rest" wavelength from the recessional velocity
              of the CALIFA header ('MED_VEL' or 'V500 MED_VEL'), if present. It will be overwritten
              if "z" is used. It can also read MANGA cubes and FITS Pycasso Files
        flip2disp : Bool
              If the data is 2D, it will transpose (flip) the Matrix according to the DISPAXIS
              axis so that a matrix in shape (wavelength, number of spectra) is provided
        rss : Bool
              It gets an RSS object instead of a Cube object. Applies only to Pycasso objects
        f2e : float or None
              Conversion factor (default: 0.1) of Flux to Error in case the ERROR is not provided:
                                      f_error = f2e * np.abs(flux)
        bflag : float or None
              Default NON-flagged value. If is not None, it creates a matrix fill with this value
        c : Float
              Speed of light in km/s. Used to obtain the redshift from the recessional velocity of
              CALIFA cubes
        dist_pc : Float
              Distance of the object in pc (to estimate luminosities)
        filters : List of list or dictionary
              filters = {'r': 'sloarn_r.txt'}   |   ['sloan_r.txt']   |   [['r', 'sloan_r.txt']]
        filter_mode : String or list of strings
              'Obs' to estimate the filter in the observed frame or 'ObsR' in rest frame
        filter_props : String or list of strings
              List of properties to calculate the magnitudes. By default:
                      filter_props = ['Flux','Mag','AB','Lum','L']
        sensf : Bool
              In case there is a SENSFUNC HDU (WEAVE), apply it to the DATA HDU
        units : Float
              Units of the flux (i.e. 1e-16)
        isivar : Bool
              Convert error HDU from IVAR to error (1/sqrt/ivar)
        verbo : Bool
              Verbosity for some Warnings

        Returns
        -------
        Object: obj.data, obj.hdr, obj.phdr, obj.flag, obj.error, obj.errorw, obj.fibcover, obj.wave, obj.wave_rest,
                      obj.redshift, obj.velocity, obj.cdelt, obj.crval, obj.nwave, self.syn, self.dispaxis
        """
        if not os.path.exists(name):
            sys.exit('*** File "' + os.path.basename(name) + '" does NOT exists!!! ***')
        try:
            hdu = pyfits.open(name)
        except:
            sys.exit('*** File "' + os.path.basename(name) + '" has a PROBLEM!!! ***')
        self.hdr = hdu[exhdr].header
        self.phdr = hdu[0].header
        self.hdu_names = [ihdu.name for ihdu in hdu]

        # Setting and reading variables
        naxis = self.hdr.get("NAXIS")
        self.califaid = self.hdr.get("CALIFAID")
        self.kalifaid = "K" + str(self.califaid).zfill(4) if self.califaid is not None else None
        self.survey = self.phdr.get("SURVEY")
        self.instrument = self.phdr.get("INSTRUME")
        self.isivar = ivar
        self.pycasso = None
        if "QVERSION" in self.phdr:
            self.pycasso = 1
        if "PYCASSO VERSION" in self.phdr:
            self.pycasso = 2
        if self.survey is not None:
            self.survey = self.survey.strip()
        if self.instrument is not None:
            self.instrument = self.instrument.strip()
        self.filter_mode = checkList(filter_mode) if filter_mode is not None else ["Obs", "ObsR"]
        self.filter_props = (
            ["Flux", "Mag", "AB", "Lum", "L"] if filter_props is None else checkList(filter_props)
        )
        self.dist_pc = dist_pc
        self.filters = filters
        self.units = units
        self.dpass = None
        self.dmag = None
        self.pversion = None
        self.dispaxis = None
        self.velocity = None
        self.redshift = None
        self.error = None
        self.ivar = None
        self.flag = None
        self.errorw = None
        self.fibcover = None
        self.sensfunc = None
        self.flat = None
        self.wave = None
        self.syn = None
        self.K = None

        # Velocity read in case of CALIFA
        # QBICK redshift. Sometimes MED_VEL is external, but is stored as CZ
        if "CZ" in list(self.hdr.keys()):
            self.velocity = self.hdr["CZ"]
        if "MED_VEL" in list(self.hdr.keys()):
            self.velocity = self.hdr["MED_VEL"]
        if "V500 MED_VEL" in list(self.hdr.keys()):  # COMBO V500+V1200
            self.velocity = self.hdr["V500 MED_VEL"]
        if hvel is not None and hvel in list(self.hdr.keys()):
            self.velocity = float(self.hdr[hvel])
        if self.velocity is not None:
            self.redshift = float(self.velocity) / c

        # Read HDUs
        self.data = hdu[0].data
        if len(hdu) > 1 and guess:
            for i in range(1, len(hdu)):
                if "EXTNAME" in list(hdu[i].header.keys()):
                    if hdu[i].header["EXTNAME"].split()[0] == "ERROR":
                        self.error = hdu[i].data
                    elif hdu[i].header["EXTNAME"].split()[0] == "BADPIX":
                        self.flag = hdu[i].data.astype("bool")
                    elif hdu[i].header["EXTNAME"].split()[0] == "ERRWEIGHT":
                        self.errorw = hdu[i].data
                    elif hdu[i].header["EXTNAME"].split()[0] == "FIBCOVER":
                        self.fibcover = hdu[i].data
                    elif hdu[i].header["EXTNAME"].split()[0] == "WAVE":
                        self.wave = hdu[i].data
                    elif hdu[i].header["EXTNAME"].split()[0] == "FLAT":
                        self.flat = hdu[i].data
                    # MaNGA
                    if self.instrument == "MaNGA":
                        if hdu[i].header["EXTNAME"].split()[0] == "FLUX":
                            self.data = hdu["FLUX"].data
                        if hdu[i].header["EXTNAME"].split()[0] == "IVAR":
                            self.ivar = hdu["IVAR"].data
                            self.error = 1.0 / np.sqrt(hdu["IVAR"].data)
                        if hdu[i].header["EXTNAME"].split()[0] == "WAVE":
                            self.wave = hdu["WAVE"].data
                        if hdu[i].header["EXTNAME"].split()[0] == "MASK":
                            self.flag = hdu["MASK"].data
                    # MUSE
                    if self.instrument == "MUSE":
                        if hdu[i].header["EXTNAME"].split()[0] == "DATA":
                            self.data = hdu["DATA"].data
                            self.hdr = hdu["DATA"].header
                            naxis = self.hdr.get("NAXIS")
                        if hdu[i].header["EXTNAME"].split()[0] == "STAT":
                            self.error = hdu["STAT"].data
                    # WEAVE
                    if isinstance(self.instrument, str) and "WEAVE" in self.instrument:
                        if hdu[i].header["EXTNAME"].split()[0].endswith("_DATA"):
                            self.data = hdu[i].data
                            self.hdr = hdu[i].header
                            naxis = self.hdr.get("NAXIS")
                        if hdu[i].header["EXTNAME"].split()[0].endswith("_IVAR"):
                            self.ivar = hdu[i].data
                            self.error = 1.0 / np.sqrt(hdu[i].data)
                        if hdu[i].header["EXTNAME"].split()[0].endswith("_SENSFUNC"):
                            self.sensfunc = hdu[i].data

        if exdata is not None:
            self.data = hdu[exdata].data
        if exflag is not None:
            self.flag = hdu[exflag].data
        if exerror is not None:
            self.error = hdu[exerror].data
        if exerw is not None:
            self.errorw = hdu[exerw].data
        if exfibc is not None:
            self.fibcover = hdu[exfibc].data
        if exflat is not None:
            self.flat = hdu[exflat].data
        if exsensf is not None:
            self.sensfunc = hdu[exsensf].data

        # Read external files for ERROR, FLAG and SENSFUNC
        self.ReadError(error_file)
        self.ReadFlag(flag_file)
        self.ReadSensf(sensf_file)

        if self.isivar and self.error is not None:
            self.error = 1.0 / np.sqrt(self.error)

        if (
            isinstance(self.instrument, str)
            and "WEAVE" in self.instrument
            and self.sensfunc is not None
            and sensf
        ):
            sensfunc = self.sensfunc
            if self.data.ndim == 2 and self.sensfunc.ndim == 1:
                sensfunc = self.sensfunc[:, np.newaxis]
            if self.data.ndim == 3 and self.sensfunc.ndim == 1:
                sensfunc = self.sensfunc[:, np.newaxis, np.newaxis]
            self.data *= sensfunc
            if self.error is not None:
                self.error *= sensfunc

        # Close file
        hdu.close()

        # Wavelength
        disp = self.hdr.get(hdisp)
        if disp is not None:
            if specaxis is None:
                specaxis = disp
        else:
            if specaxis is None:
                specaxis = 3 if (naxis == 3) else 1
        self.dispaxis = specaxis

        self.nwave = self.hdr.get("NAXIS" + str(specaxis))
        if haxis:
            scrval = hcrval + str(specaxis)
            scdelt = hcdelt + str(specaxis)
        else:
            scrval = hcrval
            scdelt = hcdelt
        self.crval = self.hdr.get(scrval)
        self.cdelt = self.hdr.get(scdelt)
        if crval is not None:
            self.crval = crval
        if cdelt is not None:
            self.cdelt = cdelt
        if (
            self.cdelt is not None
            and self.crval is not None
            and self.nwave is not None
            and self.wave is None
        ):
            self.wave = self.crval + self.cdelt * np.arange(self.nwave)
        if self.cdelt is None and self.wave is None and self.crval is not None:
            # SDSS/MANGA simulated cubes
            # if ('CTYPE3' in self.hdr.keys()) and (self.hdr['CTYPE3'].find('LOG10') != -1) and ('CD3_3' in self.hdr.keys()):
            # self.cdelt = self.hdr.get('CD3_3') # MANGA
            # self.wave = np.sort(np.power(10.,float(self.crval) - float(self.cdelt)*np.arange(self.nwave)))
            # For MUSE
            if (
                ("CTYPE3" in list(self.hdr.keys()))
                and (self.hdr["CTYPE3"].find("AWAV") != -1)
                and ("CD3_3" in list(self.hdr.keys()))
            ):
                self.cdelt = self.hdr.get("CD3_3")  # MUSE
                self.wave = self.crval + self.cdelt * np.arange(self.nwave)
            # For WEAVE *single* RSS pointings
            if (
                ("CRVAL1" in list(self.hdr.keys()))
                and ("CD1_1" in list(self.hdr.keys()))
                and (self.instrument is not None and "WEAVE" in self.instrument)
                and self.data.ndim == 2
            ):
                self.crval = self.hdr.get("CRVAL1")
                self.cdelt = self.hdr.get("CD1_1")
                self.wave = self.crval + self.cdelt * np.arange(self.nwave)
        # For PYCASSO2 pre-processed cube
        if self.pycasso == 2 and not "F_SYN" in self.hdu_names and self.nwave is not None:
            from astropy.wcs import WCS

            self.wave = get_wavelength_coordinates(WCS(self.hdr), self.nwave)

        # Read external files for WAVE
        self.ReadHDU(name, exwave)
        self.ReadWave(wave_file)

        # External redshift & Rest Wavelength
        self.RestWave(z=z, velocity=velocity)

        # Guess Pycasso Format
        if self.pycasso is not None:
            self.GuessPycasso(name)
        # Read Pycasso Format
        if self.pycasso is not None and guess:
            self.ReadPycasso(name, rss=rss)

        # Units/factor
        if self.units is not None:
            self.data *= self.units
            if self.error is not None:
                self.error *= self.units

        # Create Error and Flag Files if do not exist
        if self.error is None and f2e is not None:
            if verbo:
                print("*** No ERROR found!! We create one as (%2.1f*F_lambda)!! ***" % f2e)
            self.error_file = "No error file provided. Created one as (%2.1f*F_lambda)" % f2e
            self.error = f2e * np.abs(self.data)

        if self.flag is None and bflag is not None:
            if verbo:
                print("*** No FLAG found!! We create one and set all to %s!! ***" % bflag)
            self.flag_file = "No flag file provided. Created one and set all to %s" % bflag
            self.flag = np.zeros((self.data.shape))
            self.flag[:] = bflag

        # Transpose if required
        # if naxis == 2 and disp == 1 and flip2disp:
        if naxis == 2 and specaxis == 1 and flip2disp:
            self.data = self.data.T
            if self.error is not None:
                self.error = self.error.T
            if self.flag is not None:
                self.flag = self.flag.T
            if self.errorw is not None:
                self.errorw = self.errorw.T

        # Magnitudes
        self.setPassbands(filters=self.filters)
        self.getPassbands()

    # Functions -------------
    def OpenSimpleFits(self, name_fits, ax=0):
        im = pyfits.open(name_fits)
        data = im[ax].data
        im.close()
        return data

    def ReadError(self, error_file=None, ax=0):
        if error_file is not None:
            self.error = self.OpenSimpleFits(error_file, ax=ax)
            self.error_file = os.path.basename(error_file)
        else:
            self.error_file = None

    def ReadFlag(self, flag_file=None, ax=0):
        if flag_file is not None:
            self.flag = self.OpenSimpleFits(flag_file, ax=ax)
            self.flag_file = os.path.basename(flag_file)
        else:
            self.flag_file = None

    def ReadSensf(self, sensf_file=None, ax=0):
        if sensf_file is not None:
            self.sensfunc = self.OpenSimpleFits(sensf_file, ax=ax)
            self.sensf_file = os.path.basename(sensf_file)
        else:
            self.sensf_file = None

    def ReadWave(self, wave_file=None, ax=0):
        if wave_file is not None:
            self.wave = self.OpenSimpleFits(wave_file, ax=ax)
            self.wave_file = os.path.basename(wave_file)
        else:
            self.wave_file = None

    def ReadHDU(self, name, exhdu=None):
        if exhdu is not None:
            try:
                self.wave = pyfits.getdata(name, exhdu)
            except:
                print('*** FITS EXTENSION "%s" NOT found!!! ***' % exhdu)

    def RestWave(self, z=None, velocity=None):
        if z is not None:
            self.redshift = z
        if velocity is not None:
            self.velocity = velocity
            self.redshift = float(velocity) / (CST.c / 1e5)
        if self.redshift is not None and self.wave is not None:
            self.wave_rest = self.wave / (1.0 + self.redshift)
        else:
            self.wave_rest = None

    def GuessPycasso(self, name):
        lhpy = ["F_OBS", "F_ERR", "F_FLAG"]
        if not all(ilhpy in self.hdu_names for ilhpy in lhpy):
            self.pycasso = None

    def ReadPycasso(self, name, rss=False):
        if self.pycasso == 1:
            try:
                from pycasso.fitsdatacube import fitsQ3DataCube
                import pycasso

                K = fitsQ3DataCube(name)
                self.pversion = pycasso.__version__
                self.hdr = K.header
                self.wave_rest = K.l_obs
                self.wave = K.l_obs
                if rss:
                    self.data = K.f_obs
                    self.error = K.f_err
                    self.flag = K.f_flag
                    self.errorw = K.f_wei
                    self.syn = K.f_syn
                else:
                    self.data = K.zoneToYX(K.f_obs, extensive=False)
                    self.error = K.zoneToYX(K.f_err, extensive=False)
                    self.flag = K.zoneToYX(K.f_flag, extensive=False)
                    self.errorw = K.zoneToYX(K.f_wei, extensive=False)
                    self.syn = K.zoneToYX(K.f_syn, extensive=False)
                self.K = K
            except:
                sys.exit("*** PyCASSO FITS format. You need PyCasso module! ***")

        elif self.pycasso == 2:
            try:
                from pycasso2.cube import FitsCube
                import pycasso2

                K = FitsCube(name)
                self.pversion = pycasso2.__version__
                self.hdr = K._header
                self.wave_rest = K.l_obs
                self.wave = K.l_obs
                self.K = K
                if K.hasSegmentationMask and not rss:
                    self.data = K.spatialize(K.f_obs * K.flux_unit, extensive=False)
                    self.error = K.spatialize(K.f_err * K.flux_unit, extensive=False)
                    # New spatialize is time consuming if array is integer, so we convert to float before spatializing
                    self.flag = K.spatialize(K.f_flag.astype(float), extensive=False).astype(int)
                    # For preprocessed cubes that do not have synthetic data
                    try:
                        self.syn = K.spatialize(K.f_syn * K.flux_unit, extensive=False)
                        self.errorw = K.spatialize(K.f_wei, extensive=False)
                    except:
                        pass
                else:
                    self.data = K.f_obs * K.flux_unit
                    self.error = K.f_err * K.flux_unit
                    self.flag = K.f_flag
                    # For preprocessed cubes that do not have synthetic data
                    try:
                        self.syn = K.f_syn * K.flux_unit
                        self.errorw = K.f_wei
                    except:
                        pass
            except:
                sys.exit("*** PyCASSO2 FITS format. You need PyCasso v2 module! ***")

            if self.units is not None:
                if self.data is not None:
                    self.data *= self.units
                if self.error is not None:
                    self.error *= self.units
                if self.syn is not None:
                    self.syn *= self.units

    def setPassbands(self, filters=None):
        if filters is None:
            return
        try:
            from pyfu.passband import PassBand
        except:
            print('>>> Need "pyfu" package to convolve filters')
            return
        self.dpass = OrderedDict()
        if isinstance(filters, (str, list, tuple)):
            filters = checkList(filters)
            for filt in filters:
                if isinstance(filt, str):
                    passband = PassBand(file=filt)
                    self.dpass[passband._name_filter] = passband
                elif isinstance(filt, list):
                    passband = PassBand(file=filt[1])
                    self.dpass[filt[0]] = passband
        elif isinstance(filters, dict):
            for filt in list(filters.keys()):
                passband = PassBand(file=filters[filt])
                self.dpass[filt] = passband

    def getPassbands(
        self,
        filter_props=None,
        filter_mode=None,
        suffix=None,
        prefix=None,
        dist_pc=None,
        error_fmt="e_%s",
        **kwargs
    ):
        if self.dpass is None or len(self.dpass) < 1:
            return
        if filter_mode is None:
            filter_mode = self.filter_mode
        if filter_props is None:
            filter_props = self.filter_props
        filter_mode = checkList(filter_mode)
        filter_props = checkList(filter_props)
        dist_pc = dist_pc if dist_pc is not None else self.dist_pc
        if filter_mode is None or filter_props is None:
            return
        if self.dmag is None:
            self.dmag = OrderedDict()
        for filt in self.dpass:
            pb = self.dpass[filt]
            for mode in filter_mode:
                if not "Obs" in mode or "D" in mode:
                    continue
                wave = self.wave_rest if "R" in mode else self.wave
                if wave is None:
                    print('Wave in mode "%s" for filter "%s" NOT available' % (mode, filt))
                    continue
                flux, eflux = pb.getFluxPass(
                    wave, self.data, error=self.error, mask=self.flag, **kwargs
                )
                for prop in filter_props:
                    key = "%s_%s_%s" % (mode, prop, filt)
                    if suffix is not None:
                        key = "%s%s" % (key, suffix)
                    if prefix is not None:
                        key = "%s%s" % (prefix, key)
                    ekey = error_fmt % key
                    if "Flux" == prop:
                        self.dmag[key] = flux
                        self.dmag[ekey] = eflux
                    if "Mag" == prop:
                        mag, emag = pb.fluxToMag(flux, error=eflux, system="Vega")
                        self.dmag[key] = mag
                        self.dmag[ekey] = emag
                    if "AB" == prop:
                        mag, emag = pb.fluxToMag(flux, error=eflux, system="AB")
                        self.dmag[key] = mag
                        self.dmag[ekey] = emag
                    if ("Lum" == prop or "L" == prop) and dist_pc is not None:
                        lum, elum = pb.fluxToLuminosity(d_pc=dist_pc, flux=flux, error=eflux)
                        if "Lum" == prop:
                            self.dmag[key] = lum
                            self.dmag[ekey] = elum
                        if "L" == prop:
                            sun_lum = pb.sunLuminosity()
                            self.dmag[key] = lum / sun_lum
                            self.dmag[ekey] = elum / sun_lum if elum is not None else None


# -------------------------------------------------------------------------------------
