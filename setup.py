'''
First created on Jul 12, 2013

@author: RGB
'''

import sys, os
#sys.path.insert(0, os.path.abspath('./src'))
from viewcube.version import __version__

from viewcube.config import WriteConfigFile, viewcuberc
import shutil

configfile = os.path.join(os.environ['HOME'],'.'+viewcuberc)
WriteConfigFile(filerc=configfile)

os.chmod(configfile, 0o666)

from setuptools import setup, find_packages

setup(name='ViewCube',
      version=__version__,
      description='Datacube visualization and sonification',
      author='Ruben Garcia-Benito',
      author_email='rgb@iaa.es',
      license='MIT',
      url='https://github.com/rgbIAA/viewcube/',
      project_urls={
        "Documentation": "https://viewcube.readthedocs.io",
        "Source Code": "https://github.com/rgbIAA/viewcube/",
      },
      packages=find_packages(), #['viewcube'],
      scripts=['scripts/ViewCube.py'],
      install_requires=['numpy', 'astropy'],
      extras_require={
          "sonicube": ["hashlib",
                       "tensorflow",
                       "pythonosc",
                       "ctcsound",
                       "librosa"],
      },
      keywords=['Scientific/Astrophysics/Spectroscopy/Sonification/'],
      classifiers=[
                   "Programming Language :: Python :: 3",
                   "License :: OSI Approved :: MIT License",
                  ],
      python_requires=">=3.6",   # Minimum Python version
     )
