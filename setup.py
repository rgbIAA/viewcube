'''
Created on Jul 12, 2013

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

#pyexe = sys.executable

from setuptools import setup

setup(name='ViewCube',
      version=__version__,
      description='ViewCube',
      author='Ruben Garcia-Benito',
      author_email='rgb@iaa.es',
      license='MIT',
      url='http://rgb.iaa.es',
      #download_url=None
      packages=['viewcube'],
      provides=['viewcube'],
      scripts=['scripts/ViewCube.py'],
      requires=['numpy', 'astropy'],
      keywords=['Scientific/Engineering'],
      classifiers=[
                   "Development Status :: 4 - Beta",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: MIT License",
                  ],
     )
