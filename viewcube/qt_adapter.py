import pyqtgraph as pg
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QDialog, QVBoxLayout, QPushButton, QLabel, QDoubleSpinBox
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
import numpy as np
from astropy.io import fits

# Configuración global de pyqtgraph
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

class LambdaLimitsDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Límites Lambda")
        layout = QVBoxLayout()
        
        # Límite inferior
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(-1e6, 1e6)
        layout.addWidget(QLabel("Límite inferior:"))
        layout.addWidget(self.min_spin)
        
        # Límite superior
        self.max_spin = QDoubleSpinBox()
        self.max_spin.setRange(-1e6, 1e6)
        layout.addWidget(QLabel("Límite superior:"))
        layout.addWidget(self.max_spin)
        
        # Botón de aplicar
        apply_btn = QPushButton("Aplicar")
        apply_btn.clicked.connect(self.accept)
        layout.addWidget(apply_btn)
        
        self.setLayout(layout)

class CubeViewerAdapter:
    def __init__(self, cube_viewer):
        """Inicializa el adaptador con una instancia de CubeViewer"""
        self.cube_viewer = cube_viewer
        
        # Crear widgets de visualización
        self.spaxel_widget = pg.PlotWidget()
        self.spaxel_widget.getAxis('bottom').setPen('k')
        self.spaxel_widget.getAxis('left').setPen('k')
        self.spaxel_widget.showGrid(x=True, y=True, alpha=0.3)
        
        self.spectrum_widget = pg.PlotWidget()
        self.spectrum_widget.getAxis('bottom').setPen('k')
        self.spectrum_widget.getAxis('left').setPen('k')
        self.spectrum_widget.showGrid(x=True, y=True, alpha=0.3)
        
        # Crear elementos de visualización
        self.image_item = pg.ImageItem()
        self.spaxel_widget.addItem(self.image_item)
        
        # Añadir barra de colores
        self.colorbar = pg.ColorBarItem(
            values=(0, 1),
            colorMap='viridis',
            orientation='vertical'
        )
        self.colorbar.setImageItem(self.image_item)
        
        # Configurar interacción del spaxel
        self.spaxel_widget.scene().sigMouseMoved.connect(self.on_mouse_moved)
        self.spaxel_widget.scene().sigMouseClicked.connect(self.on_mouse_clicked)
        
        # Variable para rastrear si el botón está presionado
        self.mouse_pressed = False
        
        # Límites de longitud de onda
        self.lambda_min = None
        self.lambda_max = None
        
        # Actualizar visualizaciones iniciales
        self.update_visualizations()
    
    def update_visualizations(self):
        """Actualiza todas las visualizaciones"""
        if hasattr(self.cube_viewer, 'data'):
            try:
                # Calcular slice central para visualización inicial
                central_slice = np.nanmean(self.cube_viewer.data, axis=0)
                
                # Actualizar spaxel
                self.plot_spaxel(central_slice)
                
                # Si hay un espectro seleccionado, actualizarlo
                if self.cube_viewer.spec is not None:
                    self.plot_spectrum(
                        self.cube_viewer.wl,
                        self.cube_viewer.spec,
                        self.cube_viewer.speccom
                    )
            except Exception as e:
                print(f"Error actualizando visualizaciones: {e}")
    
    def plot_spaxel(self, data):
        """Dibuja el spaxel"""
        try:
            # Actualizar datos de la imagen
            self.image_item.setImage(data.T)
            
            # Actualizar rango de valores para la barra de colores
            vmin, vmax = np.nanmin(data), np.nanmax(data)
            self.colorbar.setLevels((vmin, vmax))
            
            # Configurar etiquetas
            self.spaxel_widget.setLabel('bottom', 'X')
            self.spaxel_widget.setLabel('left', 'Y')
            
            # Configurar aspecto
            self.spaxel_widget.getViewBox().setAspectLocked(True)
            
            # Ajustar rangos
            self.spaxel_widget.setRange(
                xRange=(0, data.shape[1]),
                yRange=(0, data.shape[0]),
                padding=0
            )
        except Exception as e:
            print(f"Error dibujando spaxel: {e}")
    
    def plot_spectrum(self, wavelength, flux, comparison=None):
        """Dibuja el espectro"""
        try:
            self.spectrum_widget.clear()
            
            # Aplicar límites de longitud de onda si existen
            if self.lambda_min is not None and self.lambda_max is not None:
                mask = (wavelength >= self.lambda_min) & (wavelength <= self.lambda_max)
                wavelength = wavelength[mask]
                flux = flux[mask]
                if comparison is not None:
                    comparison = comparison[mask]
            
            # Aplicar factores de multiplicación
            fo = getattr(self.cube_viewer, 'kwargs', {}).get('fo', 1.0)
            fc = getattr(self.cube_viewer, 'kwargs', {}).get('fc', 1.0)
            
            # Serie principal
            self.spectrum_widget.plot(
                wavelength,
                flux * fo,
                pen=pg.mkPen('b', width=1),
                name='Spectrum'
            )
            
            # Si hay datos de comparación
            if comparison is not None:
                self.spectrum_widget.plot(
                    wavelength,
                    comparison * fc,
                    pen=pg.mkPen('r', width=1),
                    name='Comparison'
                )
            
            # Configurar etiquetas
            self.spectrum_widget.setLabel('bottom', 'Wavelength (Å)')
            self.spectrum_widget.setLabel('left', 'Flux')
            
            # Ajustar rangos
            self.spectrum_widget.enableAutoRange()
        except Exception as e:
            print(f"Error dibujando espectro: {e}")
    
    def on_mouse_moved(self, pos):
        """Maneja el movimiento del ratón sobre el spaxel"""
        try:
            # Convertir posición de la escena a coordenadas de datos
            view_pos = self.spaxel_widget.getViewBox().mapSceneToView(pos)
            x, y = int(view_pos.x()), int(view_pos.y())
            
            # Verificar límites
            if (0 <= x < self.cube_viewer.nx and 
                0 <= y < self.cube_viewer.ny):
                # Obtener y mostrar espectro
                spectrum = self.cube_viewer.getSpec(x, y)
                if spectrum is not None:
                    self.plot_spectrum(
                        self.cube_viewer.wl,
                        spectrum,
                        self.cube_viewer.speccom
                    )
        except Exception as e:
            print(f"Error procesando movimiento del ratón: {e}")
    
    def on_mouse_clicked(self, event):
        """Maneja el clic del ratón sobre el spaxel"""
        if event.button() == Qt.LeftButton:
            # Actualizar estado del botón
            self.mouse_pressed = not self.mouse_pressed
            
            if not self.mouse_pressed:
                # Si se suelta el botón, mantener el último espectro mostrado
                pass
    
    def show_lambda_limits_dialog(self):
        """Muestra el diálogo de límites lambda"""
        dialog = LambdaLimitsDialog(self.spectrum_widget)
        
        # Establecer valores actuales
        if self.lambda_min is not None:
            dialog.min_spin.setValue(self.lambda_min)
        if self.lambda_max is not None:
            dialog.max_spin.setValue(self.lambda_max)
        
        if dialog.exec_() == QDialog.Accepted:
            # Actualizar límites
            self.lambda_min = dialog.min_spin.value()
            self.lambda_max = dialog.max_spin.value()
            
            # Redibujar espectro con nuevos límites
            self.plot_spectrum(
                self.cube_viewer.wl,
                self.cube_viewer.spec,
                self.cube_viewer.speccom
            )
    
    def show_residuals(self):
        """Muestra los residuos entre el espectro original y el de comparación"""
        if self.cube_viewer.speccom is not None:
            try:
                # Calcular residuos
                residuals = (self.cube_viewer.spec * self.cube_viewer.kwargs.get('fo', 1.0) - 
                           self.cube_viewer.speccom * self.cube_viewer.kwargs.get('fc', 1.0))
                
                # Crear nueva ventana para residuos
                residual_widget = pg.PlotWidget()
                residual_widget.setWindowTitle("Residuos")
                residual_widget.plot(
                    self.cube_viewer.wl,
                    residuals,
                    pen=pg.mkPen('g', width=1)
                )
                residual_widget.setLabel('bottom', 'Wavelength (Å)')
                residual_widget.setLabel('left', 'Residuals')
                residual_widget.show()
            except Exception as e:
                print(f"Error mostrando residuos: {e}") 