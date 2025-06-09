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
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
                            QLabel, QSpinBox, QComboBox, QFileDialog, QMessageBox, QGroupBox)
from PyQt5.QtCore import Qt, QRectF, QPointF
from PyQt5.QtGui import QColor, QPen, QBrush
import pyqtgraph as pg
from .utils import lsfiles, ckfiles, LoadFits, image_max_pixel
from .utils import save_spec, convert2iraf_spec, get_min_max
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.wcs import WCS
import numpy as np
import itertools
import string
import random
import math
import sys
import os

# Configuración global de pyqtgraph
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

def get_wavelength_coordinates(w, Nwave):
    w = w.sub([3])
    pix_coords = np.arange(Nwave)
    wave_coords = w.wcs_pix2world(pix_coords[:, np.newaxis], 0)
    if w.wcs.cunit[0] == "m":
        wave_coords *= 1e10
    return np.squeeze(wave_coords)

def GetIdFilter(list, filter, dfil="."):
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

def GetFluxLimits(flim):
    fmin = None
    fmax = None
    if flim is not None:
        if type(flim) not in [list, tuple] or len(flim) != 2:
            print("Flux limits should be a tuple or list of two items: ex. --> (None, 1e-18)")
        else:
            fmin, fmax = flim
    return fmin, fmax

def PRectangle(x, y, r):
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

def tmpName(prefix="tmp", char=8, suffix="fits"):
    schar = "".join(random.choice(string.ascii_letters + string.digits) for i in range(char))
    return "%s_%s.%s" % (prefix, schar, suffix)

class CubeViewer(QWidget):
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
        super().__init__()
        
        # Guardar parámetros
        self.name_fits = name_fits
        self.ptable = ptable
        self.fitscom = fitscom
        self.kwargs = kwargs
        
        # Crear layouts principales
        self.main_layout = QVBoxLayout(self)
        
        # Crear widgets para visualización
        self.spaxel_widget = pg.PlotWidget()
        self.spectrum_widget = pg.PlotWidget()
        
        # Configurar widgets
        self.setup_spaxel_widget()
        self.setup_spectrum_widget()
        
        # Cargar datos
        self.load_data()
        
        # Inicializar variables de estado
        self.selected_spaxels = []
        self.current_filter = default_filter
        self.colorbar = colorbar
        self.wlim = wlim
        self.flim = flim
        
    def setup_spaxel_widget(self):
        """Configura el widget de visualización de spaxels"""
        self.spaxel_widget.setBackground('w')
        self.spaxel_widget.showGrid(x=True, y=True)
        self.spaxel_widget.setLabel('left', 'Y')
        self.spaxel_widget.setLabel('bottom', 'X')
        
    def setup_spectrum_widget(self):
        """Configura el widget de visualización del espectro"""
        self.spectrum_widget.setBackground('w')
        self.spectrum_widget.showGrid(x=True, y=True)
        self.spectrum_widget.setLabel('left', 'Flux')
        self.spectrum_widget.setLabel('bottom', 'Wavelength')
        
    def load_data(self):
        """Carga los datos del archivo FITS"""
        try:
            # Cargar archivo FITS usando LoadFits
            self.fits_data = LoadFits(
                self.name_fits,
                exdata=self.kwargs.get('exdata'),
                exhdr=self.kwargs.get('exhdr', 0),
                exerror=self.kwargs.get('exerror'),
                exflag=self.kwargs.get('exflag'),
                ivar=self.kwargs.get('ivar', False)
            )
            
            # Actualizar visualizaciones
            self.update_visualizations()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error cargando archivo FITS: {str(e)}")
            
    def update_visualizations(self):
        """Actualiza las visualizaciones de spaxel y espectro"""
        self.update_spaxel_view()
        self.update_spectrum_view()
        
    def update_spaxel_view(self):
        """Actualiza la vista de spaxels"""
        if hasattr(self, 'fits_data') and self.fits_data.data is not None:
            # Limpiar vista anterior
            self.spaxel_widget.clear()
            
            # Crear imagen
            img = pg.ImageItem(self.fits_data.data[0])
            self.spaxel_widget.addItem(img)
            
            # Ajustar límites
            self.spaxel_widget.autoRange()
            
            # Agregar colorbar si está habilitado
            if self.colorbar:
                self.add_colorbar()
                
    def add_colorbar(self):
        """Agrega una barra de color a la vista de spaxels"""
        colorbar = pg.ColorBarItem(
            values=(np.min(self.fits_data.data[0]), np.max(self.fits_data.data[0])),
            colorMap='viridis'
        )
        colorbar.setImageItem(self.spaxel_widget.getImageItem())
        
    def update_spectrum_view(self):
        """Actualiza la vista del espectro"""
        if hasattr(self, 'fits_data') and self.fits_data.data is not None:
            # Limpiar vista anterior
            self.spectrum_widget.clear()
            
            # Si hay wavelength, usarlo, si no crear array
            if hasattr(self.fits_data, 'wave') and self.fits_data.wave is not None:
                x = self.fits_data.wave
            else:
                x = np.arange(self.fits_data.data.shape[0])
                
            # Obtener espectro promedio
            y = np.mean(self.fits_data.data, axis=(1,2))
            
            # Aplicar límites si están definidos
            if self.wlim is not None:
                wmin, wmax = GetLambdaLimits(x, wlim=self.wlim)
                mask = (x >= wmin) & (x <= wmax)
                x = x[mask]
                y = y[mask]
                
            if self.flim is not None:
                fmin, fmax = GetFluxLimits(self.flim)
                if fmin is not None:
                    y = np.maximum(y, fmin)
                if fmax is not None:
                    y = np.minimum(y, fmax)
            
            # Plotear espectro
            self.spectrum_widget.plot(x, y, pen=pg.mkPen('b', width=2))
            
    def on_spaxel_click(self, event):
        """Maneja el evento de clic en la vista de spaxels"""
        pos = event.scenePos()
        view_pos = self.spaxel_widget.getPlotItem().vb.mapSceneToView(pos)
        x, y = int(view_pos.x()), int(view_pos.y())
        
        if hasattr(self, 'fits_data') and self.fits_data.data is not None:
            if 0 <= x < self.fits_data.data.shape[2] and 0 <= y < self.fits_data.data.shape[1]:
                # Actualizar espectro seleccionado
                self.plot_selected_spectrum(x, y)
                
    def plot_selected_spectrum(self, x, y):
        """Plotea el espectro del spaxel seleccionado"""
        if hasattr(self, 'fits_data') and self.fits_data.data is not None:
            # Obtener espectro del spaxel
            spectrum = self.fits_data.data[:, y, x]
            
            # Obtener wavelength si existe
            if hasattr(self.fits_data, 'wave') and self.fits_data.wave is not None:
                wavelength = self.fits_data.wave
            else:
                wavelength = np.arange(len(spectrum))
                
            # Actualizar vista del espectro
            self.spectrum_widget.clear()
            self.spectrum_widget.plot(wavelength, spectrum, pen=pg.mkPen('r', width=2))
            
    def get_info(self):
        """Retorna información sobre el cubo cargado"""
        if hasattr(self, 'fits_data'):
            info = []
            if hasattr(self.fits_data, 'data'):
                info.append(f"Data shape: {self.fits_data.data.shape}")
            if hasattr(self.fits_data, 'wave'):
                info.append(f"Wavelength range: {self.fits_data.wave[0]:.2f} - {self.fits_data.wave[-1]:.2f}")
            return "\n".join(info)
        return "No data loaded"
        
    def PlotSpec(self):
        """Muestra el espectro actual"""
        if hasattr(self, 'fits_data') and self.fits_data.data is not None:
            # Crear una nueva ventana para el espectro
            spectrum_window = QWidget()
            layout = QVBoxLayout(spectrum_window)
            
            # Crear widget de gráfico
            plot_widget = pg.PlotWidget()
            plot_widget.setBackground('w')
            plot_widget.showGrid(x=True, y=True)
            plot_widget.setLabel('left', 'Flux')
            plot_widget.setLabel('bottom', 'Wavelength')
            
            # Obtener datos
            if hasattr(self.fits_data, 'wave') and self.fits_data.wave is not None:
                x = self.fits_data.wave
            else:
                x = np.arange(self.fits_data.data.shape[0])
                
            y = np.mean(self.fits_data.data, axis=(1,2))
            
            # Plotear datos
            plot_widget.plot(x, y, pen=pg.mkPen('b', width=2))
            
            # Agregar widget al layout
            layout.addWidget(plot_widget)
            
            # Mostrar ventana
            spectrum_window.setWindowTitle("Spectrum Viewer")
            spectrum_window.resize(800, 600)
            spectrum_window.show()
            
    def updateAx1(self, color=True):
        """Actualiza la vista principal de spaxels"""
        self.update_spaxel_view()
        
    def plotResidualMap(self):
        """Muestra el mapa de residuos"""
        if hasattr(self, 'fits_data') and self.fits_data.data is not None and self.fitscom is not None:
            try:
                # Cargar datos de comparación
                comp_data = LoadFits(self.fitscom).data
                
                # Calcular residuos
                residuals = self.fits_data.data - comp_data
                
                # Crear nueva ventana
                residual_window = QWidget()
                layout = QVBoxLayout(residual_window)
                
                # Crear widget de gráfico
                plot_widget = pg.PlotWidget()
                plot_widget.setBackground('w')
                
                # Crear imagen de residuos
                img = pg.ImageItem(residuals[0])
                plot_widget.addItem(img)
                
                # Agregar colorbar
                colorbar = pg.ColorBarItem(
                    values=(np.min(residuals[0]), np.max(residuals[0])),
                    colorMap='viridis'
                )
                colorbar.setImageItem(img)
                
                # Agregar widget al layout
                layout.addWidget(plot_widget)
                
                # Mostrar ventana
                residual_window.setWindowTitle("Residual Map")
                residual_window.resize(800, 600)
                residual_window.show()
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error calculando residuos: {str(e)}")
                
    def ResidualViewer(self, event=None):
        """Muestra el visor de residuos"""
        self.plotResidualMap()
        
    def WindowManager(self):
        """Muestra el administrador de ventanas"""
        # Crear ventana
        manager_window = QWidget()
        layout = QVBoxLayout(manager_window)
        
        # Crear controles
        # Límites de longitud de onda
        wave_group = QGroupBox("Wavelength Limits")
        wave_layout = QHBoxLayout(wave_group)
        
        wave_min = QSpinBox()
        wave_min.setRange(0, 10000)
        wave_min.setValue(4000)
        
        wave_max = QSpinBox()
        wave_max.setRange(0, 10000)
        wave_max.setValue(7000)
        
        wave_layout.addWidget(QLabel("Min:"))
        wave_layout.addWidget(wave_min)
        wave_layout.addWidget(QLabel("Max:"))
        wave_layout.addWidget(wave_max)
        
        # Botón de aplicar
        apply_button = QPushButton("Apply")
        apply_button.clicked.connect(lambda: self.apply_wave_limits(wave_min.value(), wave_max.value()))
        
        # Agregar widgets al layout
        layout.addWidget(wave_group)
        layout.addWidget(apply_button)
        
        # Mostrar ventana
        manager_window.setWindowTitle("Window Manager")
        manager_window.resize(400, 200)
        manager_window.show()
        
    def apply_wave_limits(self, wmin, wmax):
        """Aplica los límites de longitud de onda"""
        self.wlim = [wmin, wmax]
        self.update_spectrum_view()
        
    def SaveFile(self):
        """Guarda el espectro actual"""
        if hasattr(self, 'fits_data') and self.fits_data.data is not None:
            try:
                # Obtener nombre de archivo
                file_path, _ = QFileDialog.getSaveFileName(
                    self,
                    "Guardar espectro",
                    "",
                    "Text Files (*.txt);;FITS Files (*.fits)"
                )
                
                if file_path:
                    # Obtener datos
                    if hasattr(self.fits_data, 'wave') and self.fits_data.wave is not None:
                        x = self.fits_data.wave
                    else:
                        x = np.arange(self.fits_data.data.shape[0])
                        
                    y = np.mean(self.fits_data.data, axis=(1,2))
                    
                    # Guardar archivo
                    if file_path.endswith('.txt'):
                        np.savetxt(file_path, np.column_stack((x, y)))
                    else:
                        hdu = pyfits.PrimaryHDU(np.column_stack((x, y)))
                        hdu.writeto(file_path, overwrite=True)
                        
                    QMessageBox.information(self, "Éxito", "Archivo guardado correctamente")
                    
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error guardando archivo: {str(e)}")
                
    def Sonification(self):
        """Muestra la interfaz de sonificación"""
        QMessageBox.information(self, "Info", "Sonification feature not implemented in PyQt5 version yet")
