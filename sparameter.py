import matplotlib.pyplot as plt
import numpy as np
from smithplot import SmithAxes
from scipy.signal import savgol_filter
import sys
import os

"""
Get Q values from S11 touchstone format
or any other complex vector format

Based on the paper by F. Caspers
http://arxiv.org/abs/1201.4068

More info on Q measurement also in:
Microwave Measurements, E. Ginzton, McGraw-Hill 1957


Code history:

Xaratustrah 2019
"""


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
        return idx_res, f_res

    def find_transmission_resonance(self):
        # find resonant frequency
        idx_res = np.argmax(np.abs(self.cspar))
        f_res = self.freqs[idx_res]
        return idx_res, f_res

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

    def get_q_values(self):
        """
        Calculate Q values
        """
        # Calculating unloaded Q
        # find the index of the point where Re(Z) = |Im(Z)| left of the resonance and call it f5

        idx_res, f_res = self.find_reflection_resonance()
        cplx_dsp = self.get_dsp()
        impz_dsp = SParameter.get_impedance(cplx_dsp)

        offset = 100
        f5 = self.freqs[np.argmin(np.abs(
            np.abs(np.real(impz_dsp)[idx_res - offset: idx_res]) -
            np.abs(np.imag(impz_dsp)[idx_res - offset: idx_res])
        )
        )
            + idx_res - offset]

        # find the index of the point where Re(Z) = |Im(Z)| right of the resonance and call it f6
        f6 = self.freqs[np.argmin(np.abs(
            np.abs(np.real(impz_dsp)[idx_res: idx_res + offset]) -
            np.abs(np.imag(impz_dsp)[idx_res: idx_res + offset])
        )
        )
            + idx_res]

        delta_f_u = np.abs(f5 - f6)
        Qu = f_res / delta_f_u
        print('Qu: ', Qu, '∆fu: ', delta_f_u, 'MHz')

        # Calculating the loaded Q
        f1 = self.freqs[np.argmax(np.imag(cplx_dsp))]
        f2 = self.freqs[np.argmin(np.imag(cplx_dsp))]
        delta_f_l = np.abs(f1 - f2)
        Ql = f_res / delta_f_l
        print('Ql: ', Ql, '∆fl: ', delta_f_l, 'MHz')

        # Calculate the external Q from the other two Qs
        Qext_calc = 1 / (1 / Ql - 1 / Qu)
        print('Qext_calc: ', Qext_calc)

        # Calculate the coupling factor
        beta_calc = Qu / Qext_calc
        print('beta_calc', beta_calc)

        # Alternative calculation of the external Q using the impedance ~ ±j
        # you need to use the detuned open position for this
        # find the first zero crossing from beginning to half of the vector
        # https://stackoverflow.com/a/3843124/5177935
        #

        # impz_dop = self.get_dop()
        #
        # f3 = self.freqs[np.where(
        #     np.diff(np.sign(np.imag(impz_dop[:idx_res]))))[0][0]]
        # # now on the second half of the array
        # f4 = self.freqs[np.where(np.diff(np.sign(np.imag(impz_dop[idx_res:]))))[
        #     0][0] + idx_res]
        # delta_f_ext = np.abs(f3 - f4)
        # Qext_meas = f_res / delta_f_ext
        # print('Qext_meas: ', Qext_meas)
        #
        # # Calculate the coupling factor
        # beta_meas = Qu / Qext_meas
        # print('beta_meas', beta_meas)
        return Qu, Ql, delta_f_u, f_res, beta_calc

    def plot_smith(self, Qu, Ql, delta_f_u, f_res, beta_calc):
        # Plot section
        plt.figure()  # figsize=(6, 6))
        ax = plt.subplot(1, 1, 1, projection='smith', grid_minor_enable=True)
        plt.plot(self.cspar, markevery=10, label='Measurement data',
                 datatype=SmithAxes.S_PARAMETER)
        plt.plot(self.get_dsp(), markevery=10, label='Detuned short position',
                 datatype=SmithAxes.S_PARAMETER)
        plt.plot(self.get_dop(), markevery=10, label='Detuned open position',
                 datatype=SmithAxes.S_PARAMETER)
        plt.legend(loc="lower right", fontsize=8)
        plt.title(self.file_basename_wo_ext)
        print('Plotting to file {}.png'.format(self.filename_wo_ext))
        txt = 'Qu={:0.0f}, Ql={:0.0f}\n∆fu={:0.3f} [MHz]\nf_res={:0.3f} [MHz]\nbeta={:.3f}'.format(
            Qu, Ql, delta_f_u, f_res, beta_calc)
        plt.text(0, 0.6, txt, size=9, rotation=0,
                 ha="left", va="top",
                 bbox=dict(boxstyle="square",
                           ec=(1., 0.5, 0.5),
                           fc=(1., 0.8, 0.8),
                           )
                 )

        plt.savefig("{}.png".format(self.filename_wo_ext),
                    format="png",
                    dpi=600,
                    bbox_inches="tight",)

        plt.clf()

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

    Qu, Ql, delta_f_u, f_res, beta_calc = sp.get_q_values()
    sp.plot_smith(Qu, Ql, delta_f_u, f_res, beta_calc)
