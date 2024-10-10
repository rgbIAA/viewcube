#!/usr/bin/env python
# coding: utf-8
# Adrián García Riber 2023
# Rubén García Benito 2023 (class adaptation to ViewCube)

from astropy.io import fits
import numpy as np
import hashlib

import tensorflow as tf
from tensorflow import keras

from pythonosc import udp_client

import ctcsound

import librosa
import math

import os

# import session_info
# session_info.show()


class SoniCube(object):
    """Califa survey Sonification"""

    def __init__(
        self,
        fig,
        file=None,
        base_dir=None,
        verbose=False,
        port=9970,
        ref=None,
        data=None,
        flux_sensitive=True,
        **kwargs
    ):
        self.verbose = verbose
        self.flux_sensitive = flux_sensitive

        self.is_cube = False
        self.autoencoder = None

        self.base_dir = os.path.dirname(__file__) if base_dir is None else base_dir

        self.set_interactive(fig)
        self.set_sound(port=port)
        self.open_file(file, ref=ref, data=data)
        self.start_sound()

    def set_interactive(self, fig):
        fig.canvas.mpl_connect("key_press_event", self.PressKey)

    def PressKey(self, event):
        if event.key == "j":
            self.flux_sensitive = not self.flux_sensitive

    def open_file(self, file=None, ref=None, data=None):
        if isinstance(file, str):
            self.check_datacube(file)

        if not self.is_cube:
            return

        self.open_autoencoder()
        self.open_weights()
        self.preprocessing(file, ref=ref, data=data)

    def set_sound(self, port=9970):
        self.port = port  # OSC port
        self.cs = None
        self.pt = None

        self.auto_path = os.path.join(self.base_dir, "autoencoder/")
        self.weights_path = os.path.join(self.base_dir, "data/")

        self.hrtf_left = os.path.join(self.base_dir, "binaural", "hrtf-48000-left.dat")
        self.hrtf_right = os.path.join(self.base_dir, "binaural", "hrtf-48000-right.dat")

        self.csound_file = os.path.join(self.base_dir, "csound", "SoniCube.csd")

    def load_csound(self):
        if not os.path.exists(self.csound_file):
            print(">>> CSound file NOT found! [%s]" % self.csound_file)
            return
        # Open CSound
        self.cs = ctcsound.Csound()
        self.cs.compileCsd(self.csound_file)
        self.cs.start()

        # Sends the hrtf paths to CSound
        self.cs.setStringChannel("hrtf_L", self.hrtf_left)
        self.cs.setStringChannel("hrtf_R", self.hrtf_right)

    def start_thread(self):
        if self.cs is None:
            return
        # Starts the synthesizer thread
        self.pt = ctcsound.CsoundPerformanceThread(self.cs.csound())
        self.pt.play()

    def start_sound(self):
        if self.is_cube:
            self.load_csound()
            self.start_thread()

    def check_datacube(self, file):
        if not os.path.exists(file):
            print(">>> File does NOT exist!! [%s]" % file)
            return

        self.bfile = os.path.basename(file)
        self.root = self.bfile.split(".rscube")[0]

        if self.bfile.endswith(".gz"):
            print(">>> Please, use the uncompressed .fits file cube.")
            return
        else:
            self.name = self.bfile

        full_path = os.path.join(self.weights_path, self.name, "%s_Reference.npy" % self.name)

        if not os.path.exists(full_path):
            print(full_path)
            print(">>> No encoded file found. This cube is not processed for sonification yet.")
            return

        reference = np.load(full_path)
        self.lat_dim = reference[0][2]
        self.checksum = self.md5(file)
        self.is_cube = self.checksum == reference[0][1]

        if not self.is_cube:
            print("This file doesn't match the encoded cube: ", file)
        else:
            print(
                "IMPORTANT ADVICE: Some spectra could produce high sound pressure levels. Perform a first scan of the cube to adapt the volume of your headphones."
            )
            self.set_labels()
            self.check_mask()

    def check_mask(self):
        self.mask_path = os.path.join(self.base_dir, "data/", self.name, "%s.mask.fits" % self.root)
        if os.path.exists(self.mask_path):
            mask = fits.getdata(self.mask_path)
            self.mask = np.where(mask > 0, True, False)
        else:
            self.mask = None

    def apply_mask(self, data):
        if self.mask is not None:
            try:
                data[:, self.mask] = 0.0
            except:
                print(">>> Mask found NOT compatible with the datacube!")
                print(self.mask_path)

        return data

    def open_autoencoder(self):
        if not self.is_cube:
            print(">>> No encoder available!")
            return

        # Autoencoder initialization: Importing the model
        self.autoencoder = tf.keras.models.load_model(self.auto_path)

    def open_weights(self):
        if not self.is_cube or self.autoencoder is None:
            print(">>> No encoder available!")
            return

        # Importing the weights
        weights_path = os.path.join(self.weights_path, self.name, "%s_Weights" % self.name)
        self.weights = self.autoencoder.load_weights(weights_path)

    def md5(self, fname):
        """Calculates the checksum of a file"""

        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

    def preprocessing(self, file, ref=None, data=None):
        """Preprocesses and extracts variables from the cube"""

        if data is None:
            hdulist = fits.open(file)
            data = hdulist[0].data
            hdulist.close()

        # extracting data dimensions
        self.spectrum_dim, self.y_tot, self.x_tot = data.shape

        # Apply mask
        data = self.apply_mask(data)

        median_set = np.nanmedian(data, axis=0)
        self.max_median = np.nanmax(median_set)
        self.min_median = np.nanmin(median_set)

        if ref is None:
            ref = (self.y_tot / 2.0, self.x_tot / 2.0)

        self.y_ref, self.x_ref = ref

    def autoencoding(self, spectrum):
        """Calculates the latent space of a single spectrum"""

        redim_spectrum = np.zeros(
            (1, self.spectrum_dim)
        )  # We need to adapt the dimensions to the autoencoder input
        redim_spectrum[0] = spectrum  # We need to normalize and kill nans before autoencoding
        redim_spectrum[np.isnan(redim_spectrum)] = 0
        max_flux = np.nanmax(redim_spectrum)
        min_flux = np.nanmin(redim_spectrum)

        if np.nanmax(redim_spectrum) != 0:  # MAX-MIN Normalizing
            normalized_spectrum = (redim_spectrum - min_flux) / (max_flux - min_flux)
            encoded_spectrum = self.autoencoder.encoder(normalized_spectrum).numpy()  # Encoding
        else:
            normalized_spectrum = (redim_spectrum - min_flux) / 1e-16
            encoded_spectrum = self.autoencoder.encoder(normalized_spectrum).numpy()  # Encoding

        lat = encoded_spectrum * 10000  # Scaling latent values to audible frequencies

        return lat[0]

    def coord_converter(self, y, x):
        """Calculates azimuth and distance factor (normalized between 0 and 2) for the selected pixel (y,x)"""

        azimuth = (
            90 - math.degrees(math.atan2((y - self.y_ref), (x - self.x_ref))) % 360
        )  # Calculates azimuth
        radial_distance = math.sqrt(
            (y - self.y_ref) ** 2 + (x - self.x_ref) ** 2
        )  # Calculates x,y distance to x_ref, y_ref

        coordinates = [
            (0, 0),
            (0, self.y_tot),
            (self.x_tot, 0),
            (self.x_tot, self.y_tot),
        ]  # Calculates max_distance to normalize
        max_distance = 0
        for k, j in coordinates:
            distance = math.sqrt((k - self.y_ref) ** 2 + (j - self.x_ref) ** 2)
            max_distance = max(max_distance, distance)
        # Calculates a distance factor between 0 and 2 to control Reverb(distance emulator)
        norm_distance = (radial_distance / max_distance) * 2  #

        return azimuth, norm_distance

    def ponderating(
        self, frequencies
    ):  # returns the A-weighted amplification factor to compensate loudness
        """Applies the A-weighted levels curve to the input array of frequencies to return their normalized amplification factors"""

        amp_factors = []
        for freq in frequencies:
            if freq != 0:
                weight = librosa.A_weighting(freq)
                norm_weight = (weight - (-80)) / (
                    (1.27134) - (-80)
                )  # MAX-MIN normalization for 24KHz bandwith
                if (norm_weight) != 0:
                    amp_factor = round(1 / norm_weight, 2)
                else:
                    amp_factor = 0
            else:
                amp_factor = 0
            amp_factors.append(amp_factor)
        return amp_factors

    def set_labels(self):
        # Generating labels for the OSC messages
        self.amps_label = "amp0"
        self.freqs_label = "lat0"
        for dimension in range(1, self.lat_dim):
            self.amps_label += "/amp%i" % dimension
            self.freqs_label += "/lat%i" % dimension

    def sonify(self, y, x, spectrum):
        # if (0 <= x < self.x_tot) and (0 <= y < self.y_tot):
        try:
            # Calculating coordinates and radial distance
            (azimuth, dist) = self.coord_converter(y, x)
            azimuth = round(azimuth, 2)
            dist = round(dist, 2)

            lat = np.round(self.autoencoding(spectrum))

            # Creates the list of frequencies to send
            lat_vect = [round(lat[dimension]) for dimension in range(self.lat_dim)]

            # calculating median of absolut flux for current spectrum
            median_flux = round(np.ma.median(spectrum, axis=0), 5)

            if median_flux != 0:
                if self.flux_sensitive:
                    # MAX-MIN Normalizing fluxes between 1 & 10
                    norm = (median_flux - self.min_median) / (self.max_median - self.min_median)
                    median_norm = 1 + norm * 9
                else:
                    median_norm = 8  # Not too loud
            else:
                median_norm = 1
            log_norm_flux = round(np.log(median_norm), 5)

            if log_norm_flux >= 1:  # Dynamic range compression in two stages
                log_norm_flux = round(log_norm_flux / 2.5, 5)
            if 0.5 <= log_norm_flux < 1:
                log_norm_flux = round(log_norm_flux / 1.2, 5)

            # Ponderating latent vector amplitudes (A-weighted)
            amp_weights = self.ponderating(lat_vect)

            self.send_sound(amp_weights, lat_vect, azimuth, dist, log_norm_flux, median_flux)

        except:
            pass

    def send_sound(self, amp_weights, lat_vect, azimuth, dist, log_norm_flux, median_flux=None):
        # Sending OSC messages to CSound
        client_amps = udp_client.SimpleUDPClient("127.0.0.1", self.port)
        client_amps.send_message(self.amps_label, amp_weights)

        client_freqs = udp_client.SimpleUDPClient("127.0.0.1", self.port + 1)
        client_freqs.send_message(self.freqs_label, lat_vect)

        client_coords = udp_client.SimpleUDPClient("127.0.0.1", self.port + 2)
        client_coords.send_message("az/dist", [azimuth, dist])

        client_coords = udp_client.SimpleUDPClient("127.0.0.1", self.port + 3)
        client_coords.send_message("flux", log_norm_flux)

        if self.verbose:
            print("Amp. factors =", amp_weights)
            print("Latent vector =", lat_vect)
            print("Azimuth =", azimuth, "Distance =", dist)
            if median_flux is not None:
                print("Median absoult flux=", median_flux)
            print("Amplitude factor=", log_norm_flux)

    def stop_sound(self):
        client_freqs = udp_client.SimpleUDPClient("127.0.0.1", self.port + 1)
        client_freqs.send_message(self.freqs_label, (0, 0, 0, 0, 0, 0))
        client_coords = udp_client.SimpleUDPClient("127.0.0.1", self.port + 2)
        client_coords.send_message("az/dist", (0, 0))
        client_coords = udp_client.SimpleUDPClient("127.0.0.1", self.port + 3)
        client_coords.send_message("flux", 0)

    def close_sound(self):
        self.pt.stop()  # Stops the synthesizer thread
        self.cs.reset()  # Stops CSound
