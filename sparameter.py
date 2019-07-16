#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Class for handling S parameter data and touch stone format
Code history:

Xaratustrah 2019
"""

import numpy as np
from scipy.signal import savgol_filter
import sys
import os


class SParameter():
    @staticmethod
    def get_impedance(cspar):
        return (1 + cspar) / (1 - cspar) * 50

    @staticmethod
    def get_admittance(cspar):
        return 1 / SParameter.get_impedance(cspar)

    @staticmethod
    def get_vswr(cspar):
        return (1 + np.abs(cspar)) / (1 - np.abs(cspar))

    @staticmethod
    def get_dbmag_phase(cspar):
        return 20 * np.log10(np.abs(cspar)), np.rad2deg(np.angle(cspar))

    @staticmethod
    def get_linmag_phase(cspar):
        return np.abs(cspar), np.rad2deg(np.angle(cspar))

    @staticmethod
    def rotate_phase(cspar, phase_shift):
        real_cspar = np.abs(cspar) * np.cos(np.angle(cspar) + phase_shift)
        imag_cspar = np.abs(cspar) * np.sin(np.angle(cspar) + phase_shift)
        return np.vectorize(complex)(real_cspar, imag_cspar)

    def __init__(self, filename):
        self.filename = filename
        self.file_basename = os.path.basename(filename)
        self.filename_wo_ext = os.path.splitext(filename)[0]
        self.file_basename_wo_ext = self.file_basename.split('.')[0]
        self.read_touchstone()

    def read_touchstone(self):
        """
        touch stone format:
        freq, real, imaginary
        comma seperated values
        S11 data is centered on the resonance frequency

        Data format:
         HZ   S   RI   R     50.0
        ! Rohde & Schwarz ZVL
        !
        !
        !
        4.021237600000000E8    -6.366180130598420E-1  -7.604446991938126E-1
        ...
        """

        # ignore the first 5 header lines
        sri = np.genfromtxt(self.filename, skip_header=5)

        self.freqs = sri[:, 0] / 1e6
        self.cspar = np.vectorize(complex)(sri[:, 1], sri[:, 2])

    def find_reflection_resonance(self):
        # find resonant frequency
        idx_res = np.argmin(np.abs(self.cspar))
        f_res = self.freqs[idx_res]
        return f_res, idx_res

    def find_transmission_resonance(self):
        # find resonant frequency
        idx_res = np.argmax(np.abs(self.cspar))
        f_res = self.freqs[idx_res]
        return f_res, idx_res

    def get_dsp(self):
        """
        Get detuned short position
        """
        offset = 50
        # determine the rotation angle from the offset for detuned short location
        idx_max = np.argmin(np.abs(self.cspar)) + offset
        idx_min = np.argmin(np.abs(self.cspar)) - offset

        ang_min = np.angle(self.cspar[idx_min])
        ang_max = np.angle(self.cspar[idx_max])

        # make sure angles are positive
        if ang_min < 0 or ang_max < 0:
            ang_min = ang_min + 2 * np.pi
            ang_max = ang_max + 2 * np.pi

        # find the rotation angle
        rotation = (ang_min + ang_max) / 2

        # The rotation angle must always result in 0 deg so that it can
        # be added with 180 degrees in order to get dsp
        rotation = rotation if rotation < 2 * np.pi else rotation - np.pi

        # Ae^(jphi-rot) = Acos(phi-rot) + jAsin(phi-rot)
        real_dsp = np.abs(self.cspar) * \
            np.cos(np.angle(self.cspar) - rotation + np.pi)

        imag_dsp = np.abs(self.cspar) * \
            np.sin(np.angle(self.cspar) - rotation + np.pi)

        real_dsp = savgol_filter(real_dsp, 51, 3)
        imag_dsp = savgol_filter(imag_dsp, 51, 3)
        return np.vectorize(complex)(real_dsp, imag_dsp)

    def get_dop(self):
        """
        Get detuned open position
        """
        cspar_dsp = self.get_dsp()
        real_dop = np.abs(cspar_dsp) * np.cos(np.angle(cspar_dsp) - np.pi)
        imag_dop = np.abs(cspar_dsp) * np.sin(np.angle(cspar_dsp) - np.pi)
        return np.vectorize(complex)(real_dop, imag_dop)

# ------------------------


if __name__ == '__main__':
    filename = sys.argv[1]
    sp = SParameter(filename)
    if len(sys.argv) == 3:
        index = int(sys.argv[2])
        dbmag, phase = SParameter.get_dbmag_phase(sp.cspar)
        print(dbmag[index], phase[index])
        linmag, phase = SParameter.get_linmag_phase(sp.cspar)
        print(linmag[index], phase[index])
        swr = SParameter.get_vswr(sp.cspar)
        print(swr[index])
        z = SParameter.get_impedance(sp.cspar)
        print(np.real(z)[index], np.imag(z)[index])
        y = SParameter.get_admittance(sp.cspar)
        print(np.real(y)[index], np.imag(y)[index])
