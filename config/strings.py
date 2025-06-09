"""Text configurations for GUI"""
import viewcube.version as version

# General GUI text
WINDOW_TITLE = "ViewCube " + version.__version__
DEFAULT_TAB_NAME = "Workspace Tab"
SPECTRAL_FRAME_NAME = "Spectral Viewer"
SPAXEL_FRAME_NAME = "Spaxel Viewer"

# Menu bar text
FILE_MENU = "Archivo"
TOOLS_MENU = "Herramientas"
MANAGER_MENU = "Administrador"

# Menu items
OPEN_FITS = "Abrir archivo FITS"
NEW_WORKSPACE = "Crear nuevo espacio de trabajo"
SAVE_SPECTRUM = "Guardar espectro"
WINDOW_MANAGER = "Administrador de ventanas"
LAMBDA_LIMITS = "Límites Lambda"
SONIFICATION = "Sonificación"
VIEW_RESIDUALS = "Ver residuos"
FIT_SPECTRUM = "Ajustar espectro"

# Frames and layout
CONFIGURATION_FRAME = "Configuración"
PLOTS_FRAME = "Gráficas"

# Configuration widgets
LABEL_FILE_PATH = "FITS File Path"
BUTTON_FILE_PATH = "Buscar fichero"
BUTTON_SELECT_FILES = "Seleccionar archivos"
LABEL_NO_FILE = "Archivo actual: ninguno"

# Messages
ERROR_NO_DATA = "No se pudieron leer datos del archivo FITS"
ERROR_NOT_3D = "¡Este archivo FITS no es un cubo 3D!"
ERROR_LOADING = "Error al cargar el archivo: {}"
ERROR_NO_COMPARISON = "Para ver residuos, primero debe cargar un archivo de comparación"

# Plot labels
SPAXEL_TITLE = "Spaxel (Canal {})"
SPECTRUM_TITLE = "Espectro en ({}, {})"
SPECTRUM_SELECT = "Seleccione un punto en el spaxel"
SPECTRUM_REGION = "Espectro promedio de región ({},{}) a ({},{})"
LABEL_CHANNEL = "Canal"
LABEL_INTENSITY = "Intensidad"

