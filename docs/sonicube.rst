.. _sonicube:

********
SoniCube
********

ViewCube has also implemented a 
`sonification <https://rgb.iaa.es/que-es-la-sonificacion/>`_, 
module, named **SoniCube**, which provides an open, 
comprehensive, and general-purpose multimodal tool for IFS analysis, supporting the development 
of future sonification techniques focused on specific potential features within datacubes. 

The aim of the SoniCube interface is to offer a diverse "palette of sonifications", akin to the 
range of color palettes available for visual 2D maps in ViewCube. This variety of sonification 
options will allow users to extract or enhance different data characteristics by selecting 
specific sonification methods, much like how a color palette reveals visual details. 

A comprehensive description of the sonification module, along with its evaluation, is available 
in the following paper:

Adrián García Riber, Rubén García-Benito, Francisco Serradilla. *Interactive multimodal Integral 
Field Spectroscopy*, RAS Techniques and Instruments, 2024, 
`https://doi.org/10.1093/rasti/rzae049 <https://doi.org/10.1093/rasti/rzae049>`_ 


.. vimeo:: 1005208084
   :width: 100%


Installing Dependencies
=======================

Two libraries are required to make the SoniCube module functional.

1. **Csound**: An open-source sound and music computing system for audio synthesis and signal 
processing. You can download the version for your operating system from the official webpage:

   `https://csound.com/download.html <https://csound.com/download.html>`_

   Csound is available for all major operating systems.

2. **Python Libraries**: Ensure the following Python libraries are installed (if not already):

.. code-block:: bash

   hashlib
   tensorflow
   pythonosc
   ctcsound
   librosa

These libraries can be installed using `pip`.

CALIFA sonicubes
================

In this version, we introduce the first sonification method integrated into SoniCube's toolset. 
This initial implementation is intended to replicate the functionality of the visual ViewCube, 
offering a quick qualitative overview of the datacube.

SoniCube’s sonification utilizes an autoencoder architecture to accurately transform the spectral 
information of CALIFA survey COMBO datacubes into sound.

Configuration
^^^^^^^^^^^^^

To access the sonification, users need an original datacube from CALIFA Data Release 3 
(in COMBO grating) along with its corresponding autoencoder-generated sonification files.

Download the CALIFA DR3 COMBO datacubes from: 

`https: //califa.caha.es/FTP-PUB/reduced/COMB/reduced_v2.2/ <https: //califa.caha.es/FTP-PUB/reduced/COMB/reduced_v2.2/>`_

and the corresponding sonicubes from:

`https://zenodo.org/records/10570065 <https://zenodo.org/records/10570065>`_ 

Create a folder (e.g. ``sonicube``) at any location of your choice. Inside this folder, 
deploy the contents from the `sonicube` folder of the 
`ViewCube GitHub repository <https://github.com/rgbIAA/viewcube/>`_. 
The directory should contain three subfolders for sound configuration files, along 
with a data folder for the CALIFA sonicubes.

For example, if your folder is named ``sonicube``, it should contain the following structure:

.. code-block:: bash

   autoencoder/  binaural/  csound/  data/

The subfolders ``autoencoder``, ``binaural``, and ``csound`` are obtained from 
the ``sonicube`` directory in the ViewCube GitHub repository. The ``data`` folder, 
on the other hand, should contain the CALIFA sonicubes, which can be downloaded 
from the `Zenodo repository <https://zenodo.org/records/10570065>`_. 

For instance, your ``data`` folder might look like:

.. code-block:: bash

   ...
   NGC5732.COMB.rscube.fits/
   NGC5784.COMB.rscube.fits/
   NGC5794.COMB.rscube.fits/
   NGC5797.COMB.rscube.fits/
   ...

and the complete contents within each directory might be:

.. code-block:: bash

  ...
  NGC5732.COMB.rscube.fits:
  NGC5732.COMB.rscube.fits_Reference.npy
  NGC5732.COMB.rscube.fits_Weights.data-00000-of-00001
  NGC5732.COMB.rscube.fits_Weights.index
  NGC5732.COMB.rscube.fits_learning_rate.png
  checkpoint
  
  NGC5784.COMB.rscube.fits:
  NGC5784.COMB.rscube.fits_Reference.npy
  NGC5784.COMB.rscube.fits_Weights.data-00000-of-00001
  NGC5784.COMB.rscube.fits_Weights.index
  NGC5784.COMB.rscube.fits_learning_rate.png
  checkpoint
  
  NGC5794.COMB.rscube.fits:
  NGC5794.COMB.rscube.fits_Reference.npy
  NGC5794.COMB.rscube.fits_Weights.data-00000-of-00001
  NGC5794.COMB.rscube.fits_Weights.index
  NGC5794.COMB.rscube.fits_learning_rate.png
  checkpoint
  
  NGC5797.COMB.rscube.fits:
  NGC5797.COMB.rscube.fits_Reference.npy
  NGC5797.COMB.rscube.fits_Weights.data-00000-of-00001
  NGC5797.COMB.rscube.fits_Weights.index
  NGC5797.COMB.rscube.fits_learning_rate.png
  checkpoint
  ...

You do not need to download all of the sonicubes—just those that are relevant to your project.

Now in the ViewCube configuration file ``.viewcuberc`` uncomment the keyword ``dsoni`` and 
write the absolute path of the ``sonicube`` directoy:


.. code-block:: bash
   
   dsoni : "/my/absolute/path/to/sonicube/"

Sounding [Data/Soni]Cubes
^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have installed the dependencies and configured ViewCube, you need at least one
original CALIFA datacube and its corresponding sonicube. The datacube can be located in 
any directory, while the sonicube should be placed as explained in the previous 
configuration section.

Open the datacube in the standard way using ViewCube.

Ensure you are wearing your headphones correctly (right headphone on the right ear and
left on the left). Adjust the volume to a lower or medium setting to avoid high volumes 
at the beginning. You can increase the volume later, once you are familiar with the sound 
of that particular galaxy.

To activate sonification mode, press the ``h`` key. The first time you enable this mode, 
it may take 5-8 seconds to load the necessary libraries (``tensorflow`` is known to take 
some time to import).

Explore as usual in ViewCube by moving the mouse over the spaxel window. If the mouse
moves outside the axis or window, no sound will be produced.

In standard mode, the volume corresponds to the median intensity of the spectrum. Dim 
regions of the galaxy will have a lower volume, and the sky will be practically silent. 
In contrast, high-intensity HII regions or the galaxy's center will produce louder sounds.

If you want to deactivate the volume intensity linkage and have a uniform volume for all
spaxels, press the ``j`` key. This is useful if there is an interesting region with lower 
flux that you want to listen to carefully but its volume is relatively low.

Press the same key again to restore intensity-sensitive volume.
