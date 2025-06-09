import os
import platform
import sys

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QMainWindow, QApplication, QWidget, QTabWidget,
                             QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSpinBox, QComboBox,
                             QFrame, QFileDialog, QAction, QGroupBox, QMessageBox, QDoubleSpinBox, QMenu)
from astropy.io import fits

from config import strings
from viewcube import cubeviewer as cv
from viewcube.qt_adapter import CubeViewerAdapter

# Configuración global de pyqtgraph
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

def fits_file_info(file_path):
    """
    Muestra información del archivo FITS
    """
    info = []
    with fits.open(file_path) as hdul:
        info.append("Extensiones:")
        for i, hdu in enumerate(hdul):
            info.append(f"{i}: {hdu.__class__.__name__}, shape={hdu.data.shape if hasattr(hdu, 'data') else 'No data'}")
    return "\n".join(info)


def file_extensions(file_path):
    """
    :param file_path:  Path to the FITS file
    :return:  Number of extensions in the FITS file
    """
    return len(fits.open(file_path))


class PlotWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        
        # Crear el gráfico
        self.chart = pg.PlotWidget()
        self.chart.setBackgroundVisible(False)
        
        # Añadir la vista al layout
        self.layout.addWidget(self.chart)
        
        # Series de datos
        self.series = None
        
    def clear(self):
        """Limpia el gráfico"""
        self.chart.clear()
        self.series = None
        
    def plot(self, x, y, color=Qt.blue, clear=True):
        """Dibuja una serie de datos"""
        if clear:
            self.clear()
            
        # Crear nueva serie
        self.series = pg.PlotDataItem(x=x, y=y, pen=pg.mkPen(color))
        
        # Añadir serie al gráfico
        self.chart.addItem(self.series)
        
    def set_labels(self, title="", x_label="", y_label=""):
        """Establece las etiquetas del gráfico"""
        self.chart.setTitle(title)
        self.chart.setLabel('left', y_label)
        self.chart.setLabel('bottom', x_label)


class SpaxelWidget(PlotWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.scatter = pg.ScatterPlotItem()
        self.scatter.setData(x=[], y=[], pen=None, symbol='o', size=10, symbolBrush=(255, 0, 0))
        self.chart.addItem(self.scatter)
        
    def plot_spaxel(self, data, clear=True):
        """Dibuja el spaxel"""
        if clear:
            self.clear()
            
        # Crear serie para la imagen
        self.series = pg.PlotDataItem(x=np.arange(data.shape[1]), y=data[0], pen=None)
        
        # Añadir serie al gráfico
        self.chart.addItem(self.series)
        
        # Actualizar rangos
        self.chart.setRange(xRange=(0, data.shape[1]), yRange=(0, data.shape[0]))


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.cube = None
        self.cube_adapter = None
        self.data = None
        self.comparison_cube = None
        self.position_table = None

        # Initial window configuration
        self.setWindowTitle(strings.WINDOW_TITLE)
        self.setup_window_size()

        # Central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Menu
        self.setup_menu()

        # Notebook (QTabWidget)
        self.tab_widget = QTabWidget()
        main_layout.addWidget(self.tab_widget)

        # Tab
        self.create_first_tab()

        # Menu actions
        self.connect_menu_actions()

        # Mostrar la ventana maximizada
        self.showMaximized()

    def setup_window_size(self):
        if platform.system() == "Windows":
            self.showMaximized()
        else:
            screen = QApplication.primaryScreen().geometry()
            self.setGeometry(0, 0, screen.width(), screen.height())

    def setup_menu(self):
        menubar = self.menuBar()

        # Menú Archivo
        file_menu = menubar.addMenu(strings.FILE_MENU)
        file_menu.setObjectName(strings.FILE_MENU)

        open_action = QAction(strings.OPEN_FITS, self)
        open_action.setShortcut('Ctrl+O')
        file_menu.addAction(open_action)

        new_workspace_action = QAction(strings.NEW_WORKSPACE, self)
        new_workspace_action.setShortcut('Ctrl+N')
        file_menu.addAction(new_workspace_action)

        save_action = QAction(strings.SAVE_SPECTRUM, self)
        save_action.setShortcut('Ctrl+S')
        file_menu.addAction(save_action)

        # Menú Administrador
        manager_menu = menubar.addMenu(strings.MANAGER_MENU)
        manager_menu.setObjectName(strings.MANAGER_MENU)

        window_manager_action = QAction(strings.WINDOW_MANAGER, self)
        window_manager_action.setShortcut('Ctrl+W')
        manager_menu.addAction(window_manager_action)

        lambda_limits_action = QAction(strings.LAMBDA_LIMITS, self)
        lambda_limits_action.setShortcut('Ctrl+L')
        manager_menu.addAction(lambda_limits_action)

        # Menú Herramientas
        tools_menu = menubar.addMenu(strings.TOOLS_MENU)
        tools_menu.setObjectName(strings.TOOLS_MENU)

        sonification_action = QAction(strings.SONIFICATION, self)
        sonification_action.setShortcut('Ctrl+M')
        tools_menu.addAction(sonification_action)

        residual_action = QAction(strings.VIEW_RESIDUALS, self)
        residual_action.setShortcut('Ctrl+R')
        tools_menu.addAction(residual_action)

        fit_action = QAction(strings.FIT_SPECTRUM, self)
        fit_action.setShortcut('Ctrl+F')
        tools_menu.addAction(fit_action)

    def create_first_tab(self):
        tab = QWidget()
        self.tab_widget.addTab(tab, strings.DEFAULT_TAB_NAME)

        # Main layout of the tab
        tab_layout = QHBoxLayout(tab)

        # Configuration's frame
        config_group = QGroupBox(strings.CONFIGURATION_FRAME)
        config_layout = QVBoxLayout(config_group)
        self.setup_config_panel(config_layout)  # components

        # Plot's frame
        plots_group = QGroupBox(strings.PLOTS_FRAME)
        workspace_layout = QVBoxLayout(plots_group)
        self.setup_workspace_panel(workspace_layout)  # components

        # Add frames to main layout
        tab_layout.addWidget(config_group, 1)
        tab_layout.addWidget(plots_group, 3)

    def setup_config_panel(self, layout):
        # FITS File Path
        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel(strings.LABEL_FILE_PATH))
        self.btn_search_fits = QPushButton(strings.BUTTON_FILE_PATH)
        self.btn_search_fits.clicked.connect(self.on_search_fits)
        file_layout.addWidget(self.btn_search_fits)
        
        # Mostrar ruta del archivo seleccionado
        self.fits_path_label = QLabel("")
        file_layout.addWidget(self.fits_path_label)
        
        layout.addLayout(file_layout)
        
        # Tabla de posiciones
        pos_table_layout = QHBoxLayout()
        pos_table_layout.addWidget(QLabel('Tabla de posiciones externa para RSS Viewer'))
        self.btn_search_table = QPushButton("Buscar tabla")
        self.btn_search_table.clicked.connect(self.on_search_table)
        pos_table_layout.addWidget(self.btn_search_table)
        layout.addLayout(pos_table_layout)
        
        # Fichero FITS de comparación
        comp_file_layout = QHBoxLayout()
        comp_file_layout.addWidget(QLabel('Fichero FITS de comparación'))
        self.btn_search_comp = QPushButton("Buscar archivo")
        self.btn_search_comp.clicked.connect(self.on_search_comparison)
        comp_file_layout.addWidget(self.btn_search_comp)
        layout.addLayout(comp_file_layout)
        
        # DATA extension
        data_ext_layout = QHBoxLayout()
        data_ext_layout.addWidget(QLabel('DATA extension (default: None)'))
        self.data_ext_spin = QSpinBox()
        self.data_ext_spin.setRange(0, 100)
        data_ext_layout.addWidget(self.data_ext_spin)
        layout.addLayout(data_ext_layout)
        
        # ERROR extension
        error_ext_layout = QHBoxLayout()
        error_ext_layout.addWidget(QLabel('ERROR extension (default: None)'))
        self.error_ext_spin = QSpinBox()
        self.error_ext_spin.setRange(0, 100)
        error_ext_layout.addWidget(self.error_ext_spin)
        layout.addLayout(error_ext_layout)
        
        # FLAG/MASK extension
        flag_ext_layout = QHBoxLayout()
        flag_ext_layout.addWidget(QLabel('FLAG/MASK extension (default: None)'))
        self.flag_ext_spin = QSpinBox()
        self.flag_ext_spin.setRange(0, 100)
        flag_ext_layout.addWidget(self.flag_ext_spin)
        layout.addLayout(flag_ext_layout)
        
        # HEADER extension
        header_ext_layout = QHBoxLayout()
        header_ext_layout.addWidget(QLabel('HEADER extension (default: 0)'))
        self.header_ext_spin = QSpinBox()
        self.header_ext_spin.setRange(0, 100)
        self.header_ext_spin.setValue(0)
        header_ext_layout.addWidget(self.header_ext_spin)
        layout.addLayout(header_ext_layout)
        
        # Factor de multiplicación original
        fo_factor_layout = QHBoxLayout()
        fo_factor_layout.addWidget(QLabel('Multiplicative factor for original file'))
        self.fo_factor_spin = QDoubleSpinBox()
        self.fo_factor_spin.setRange(0, 100)
        self.fo_factor_spin.setSingleStep(0.1)
        self.fo_factor_spin.setValue(1.0)
        fo_factor_layout.addWidget(self.fo_factor_spin)
        layout.addLayout(fo_factor_layout)
        
        # Factor de multiplicación comparación
        fc_factor_layout = QHBoxLayout()
        fc_factor_layout.addWidget(QLabel('Multiplicative factor for comparison file'))
        self.fc_factor_spin = QDoubleSpinBox()
        self.fc_factor_spin.setRange(0, 100)
        self.fc_factor_spin.setSingleStep(0.1)
        self.fc_factor_spin.setValue(1.0)
        fc_factor_layout.addWidget(self.fc_factor_spin)
        layout.addLayout(fc_factor_layout)
        
        # IVAR to error checkbox
        ivar_layout = QHBoxLayout()
        ivar_layout.addWidget(QLabel('Conversion from IVAR to error'))
        self.ivar_combo = QComboBox()
        self.ivar_combo.addItems(["False", "True"])
        ivar_layout.addWidget(self.ivar_combo)
        layout.addLayout(ivar_layout)
        
        # Botón de carga
        load_layout = QHBoxLayout()
        self.btn_load = QPushButton("Cargar")
        self.btn_load.clicked.connect(self.on_load_clicked)
        self.btn_load.setEnabled(False)  # Deshabilitado hasta que se seleccione un archivo
        load_layout.addWidget(self.btn_load)
        layout.addLayout(load_layout)

    def setup_workspace_panel(self, layout):
        """Configura el panel de trabajo"""
        # Área de gráficos
        plots_layout = QHBoxLayout()
        
        # Frame izquierdo (Spaxel)
        spaxel_frame = QFrame()
        spaxel_frame.setFrameStyle(QFrame.Box | QFrame.Raised)
        spaxel_layout = QVBoxLayout(spaxel_frame)
        
        # Frame derecho (Espectro)
        spectrum_frame = QFrame()
        spectrum_frame.setFrameStyle(QFrame.Box | QFrame.Raised)
        spectrum_layout = QVBoxLayout(spectrum_frame)
        
        # Añadir los frames al layout principal
        plots_layout.addWidget(spaxel_frame)
        plots_layout.addWidget(spectrum_frame)
        
        # Añadir el layout de gráficos al layout principal
        layout.addLayout(plots_layout)
        
        # Configurar las proporciones del layout de gráficos
        plots_layout.setStretch(0, 1)  # Spaxel
        plots_layout.setStretch(1, 1)  # Spectrum
        
        # Guardar referencias a los layouts para usar más tarde
        self.spaxel_layout = spaxel_layout
        self.spectrum_layout = spectrum_layout

    def connect_menu_actions(self):
        """Conecta las acciones del menú con sus funciones correspondientes"""
        # Obtener las acciones del menú
        menubar = self.menuBar()
        file_menu = menubar.findChild(QMenu, 'Archivo')
        manager_menu = menubar.findChild(QMenu, 'Administrador')
        tools_menu = menubar.findChild(QMenu, 'Herramientas')

        # Conectar acciones del menú Archivo
        for action in file_menu.actions():
            if action.text() == 'Abrir archivo FITS':
                action.triggered.connect(self.on_search_fits)
            elif action.text() == 'Crear nuevo espacio de trabajo':
                action.triggered.connect(self.create_new_workspace)
            elif action.text() == 'Guardar espectro':
                action.triggered.connect(self.save_spectrum)

        # Conectar acciones del menú Administrador
        for action in manager_menu.actions():
            if action.text() == 'Administrador de ventanas':
                action.triggered.connect(self.show_window_manager)
            elif action.text() == 'Límites Lambda':
                action.triggered.connect(self.show_lambda_limits)

        # Conectar acciones del menú Herramientas
        for action in tools_menu.actions():
            if action.text() == 'Sonificación':
                action.triggered.connect(self.show_sonification)
            elif action.text() == 'Ver residuos':
                action.triggered.connect(self.show_residuals)
            elif action.text() == 'Ajustar espectro':
                action.triggered.connect(self.fit_spectrum)

    def create_new_workspace(self):
        """Crea un nuevo espacio de trabajo"""
        # Limpiar datos existentes
        self.cube = None
        self.cube_adapter = None
        self.data = None
        self.comparison_cube = None
        self.position_table = None

        # Limpiar gráficos
        if self.cube_adapter:
            self.cube_adapter.spaxel_widget.clear()
            self.cube_adapter.spectrum_widget.clear()

        # Resetear título
        self.setWindowTitle(strings.WINDOW_TITLE)

    def show_window_manager(self):
        """Muestra el administrador de ventanas"""
        if self.cube:
            self.cube.WindowManager()

    def show_lambda_limits(self):
        """Muestra el diálogo de límites lambda"""
        if self.cube_adapter:
            self.cube_adapter.show_lambda_limits_dialog()
        else:
            QMessageBox.warning(self, "Error", "Primero debe cargar un archivo FITS")

    def on_search_fits(self):
        """Maneja el evento de búsqueda de archivo FITS"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Seleccionar archivo FITS",
            "",
            "FITS files (*.fits *.fit)"
        )
        if file_path:
            self.selected_fits_file = file_path
            self.fits_path_label.setText(os.path.basename(file_path))
            self.btn_load.setEnabled(True)
    
    def on_load_clicked(self):
        """Maneja el evento de clic en el botón de carga"""
        if hasattr(self, 'selected_fits_file'):
            try:
                # Crear el diccionario de kwargs con la configuración actual
                kwargs = {
                    'exdata': self.data_ext_spin.value() if self.data_ext_spin.value() > 0 else None,
                    'exhdr': self.header_ext_spin.value(),
                    'exerror': self.error_ext_spin.value() if self.error_ext_spin.value() > 0 else None,
                    'exflag': self.flag_ext_spin.value() if self.flag_ext_spin.value() > 0 else None,
                    'fo': self.fo_factor_spin.value(),
                    'fc': self.fc_factor_spin.value(),
                    'ivar': self.ivar_combo.currentText() == "True"
                }
                
                # Crear instancia de CubeViewer con los parámetros actuales
                self.cube = cv.CubeViewer(
                    name_fits=self.selected_fits_file,
                    ptable=self.position_table,
                    fitscom=self.comparison_cube,
                    **kwargs
                )
                
                # Limpiar los layouts existentes
                for i in reversed(range(self.spaxel_layout.count())): 
                    widget = self.spaxel_layout.itemAt(i).widget()
                    if widget is not None:
                        widget.setParent(None)
                
                for i in reversed(range(self.spectrum_layout.count())):
                    widget = self.spectrum_layout.itemAt(i).widget()
                    if widget is not None:
                        widget.setParent(None)
                
                # Crear el adaptador para PyQtGraph
                self.cube_adapter = CubeViewerAdapter(self.cube)
                
                # Añadir los nuevos widgets
                self.spaxel_layout.addWidget(self.cube_adapter.spaxel_widget)
                self.spectrum_layout.addWidget(self.cube_adapter.spectrum_widget)
                
                # Actualizar título de la ventana
                self.setWindowTitle(f"{strings.WINDOW_TITLE} - {os.path.basename(self.selected_fits_file)}")
                
                # Mostrar información del archivo
                try:
                    info = [
                        fits_file_info(self.selected_fits_file),
                        "",
                        "Información del cubo:",
                        self.cube.get_info()
                    ]
                    QMessageBox.information(self, "Información del archivo", "\n".join(info))
                except Exception as e:
                    print(f"Error mostrando información: {e}")
                
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error al cargar el archivo:\n{str(e)}")
                self.data = None
                self.cube = None
                self.cube_adapter = None
                self.setWindowTitle(strings.WINDOW_TITLE)
        else:
            QMessageBox.warning(self, "Error", "Por favor, seleccione un archivo FITS primero")

    def on_search_table(self):
        """Maneja el evento de búsqueda de tabla de posiciones"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Seleccionar tabla de posiciones",
            "",
            "All Files (*.*)"
        )
        if file_path:
            try:
                # Cargar la tabla de posiciones
                self.position_table = file_path

                # Si hay un cubo cargado, actualizarlo con la nueva tabla
                if self.cube:
                    self.cube.ptable = self.position_table
                    self.update_visualizations()

                QMessageBox.information(self, "Éxito", "Tabla de posiciones cargada correctamente")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error al cargar la tabla de posiciones: {str(e)}")

    def on_search_comparison(self):
        """Maneja el evento de búsqueda de archivo FITS de comparación"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Seleccionar archivo FITS de comparación",
            "",
            "FITS files (*.fits *.fit)"
        )
        if file_path:
            try:
                # Guardar la referencia al archivo de comparación
                self.comparison_cube = file_path
                
                # Si hay un cubo cargado, actualizarlo con el nuevo archivo de comparación
                if self.cube:
                    self.cube.fitscom = self.comparison_cube
                    # Recargar el cubo para actualizar la comparación
                    self.load_fits_file(self.cube.name_fits)
                
                QMessageBox.information(self, "Éxito", "Archivo de comparación cargado correctamente")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error al cargar el archivo de comparación: {str(e)}")
                self.comparison_cube = None

    def update_cube_parameters(self):
        """Actualiza los parámetros del cubo cuando cambian los controles"""
        if self.cube:
            try:
                # Actualizar extensiones
                self.cube.kwargs['exdata'] = self.data_ext_spin.value() if self.data_ext_spin.value() > 0 else None
                self.cube.kwargs['exhdr'] = self.header_ext_spin.value()
                self.cube.kwargs['exerror'] = self.error_ext_spin.value() if self.error_ext_spin.value() > 0 else None
                self.cube.kwargs['exflag'] = self.flag_ext_spin.value() if self.flag_ext_spin.value() > 0 else None

                # Actualizar factores de multiplicación
                self.cube.kwargs['fo'] = self.fo_factor_spin.value()
                self.cube.kwargs['fc'] = self.fc_factor_spin.value()

                # Actualizar otros parámetros
                self.cube.kwargs['ivar'] = self.ivar_combo.currentText() == "True"

                # Recargar el cubo con los nuevos parámetros
                self.load_fits_file(self.cube.name_fits)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Error actualizando parámetros: {str(e)}")

    def load_fits_file(self, file_path):
        """Carga un archivo FITS y actualiza las visualizaciones"""
        try:
            # Crear instancia de CubeViewer con los parámetros de la interfaz
            self.cube = cv.CubeViewer(
                name_fits=file_path,
                ptable=self.position_table,
                fitscom=self.comparison_cube,
                exdata=self.data_ext_spin.value() if self.data_ext_spin.value() > 0 else None,
                exhdr=self.header_ext_spin.value(),
                exerror=self.error_ext_spin.value() if self.error_ext_spin.value() > 0 else None,
                exflag=self.flag_ext_spin.value() if self.flag_ext_spin.value() > 0 else None,
                ivar=self.ivar_combo.currentText() == "True"
            )
            
            # Limpiar los layouts existentes
            for i in reversed(range(self.spaxel_layout.count())): 
                widget = self.spaxel_layout.itemAt(i).widget()
                if widget is not None:
                    widget.setParent(None)
            
            for i in reversed(range(self.spectrum_layout.count())):
                widget = self.spectrum_layout.itemAt(i).widget()
                if widget is not None:
                    widget.setParent(None)
            
            # Crear el adaptador para PyQtGraph
            self.cube_adapter = CubeViewerAdapter(self.cube)
            
            # Añadir los nuevos widgets
            self.spaxel_layout.addWidget(self.cube_adapter.spaxel_widget)
            self.spectrum_layout.addWidget(self.cube_adapter.spectrum_widget)
            
            # Actualizar título de la ventana
            self.setWindowTitle(f"{strings.WINDOW_TITLE} - {os.path.basename(file_path)}")
            
            # Mostrar información del archivo
            try:
                info = [
                    fits_file_info(file_path),
                    "",
                    "Información del cubo:",
                    self.cube.get_info()
                ]
                QMessageBox.information(self, "Información del archivo", "\n".join(info))
            except Exception as e:
                print(f"Error mostrando información: {e}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error al cargar el archivo:\n{str(e)}")
            self.data = None
            self.cube = None
            self.cube_adapter = None
            self.setWindowTitle(strings.WINDOW_TITLE)

    def update_visualizations(self):
        """Actualiza las visualizaciones de spaxel y espectro"""
        if self.cube_adapter:
            self.cube_adapter.update_visualizations()

    def save_spectrum(self):
        """Guarda el espectro actual"""
        if self.cube:
            self.cube.SaveFile()

    def show_sonification(self):
        """Muestra la interfaz de sonificación"""
        if self.cube:
            self.cube.Sonification()

    def show_residuals(self):
        """Muestra el visor de residuos"""
        if self.cube_adapter and self.comparison_cube:
            self.cube_adapter.show_residuals()
        else:
            QMessageBox.warning(self, "Error",
                              "Para ver residuos, primero debe cargar un archivo de comparación")

    def fit_spectrum(self):
        """Muestra la interfaz de ajuste de espectro"""
        if self.cube:
            self.cube.FitSpec(None)


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main() 