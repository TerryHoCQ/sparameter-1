# -*- coding: utf-8 -*-

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

from sparameter import SParameter
import numpy as np


def get_unloaded_q(sparam):
    """
    Calculate Q values
    """
    # Calculating unloaded Q
    # find the index of the point where Re(Z) = |Im(Z)| left of the resonance and call it f5

    f_res, idx_res = sparam.find_reflection_resonance()
    cplx_dsp = sparam.get_dsp()
    impz_dsp = SParameter.get_impedance(cplx_dsp)

    offset = 100
    f5 = sparam.freqs[np.argmin(np.abs(
        np.abs(np.real(impz_dsp)[idx_res - offset: idx_res]) -
        np.abs(np.imag(impz_dsp)[idx_res - offset: idx_res])
    )
    )
        + idx_res - offset]

    # find the index of the point where Re(Z) = |Im(Z)| right of the resonance and call it f6
    f6 = sparam.freqs[np.argmin(np.abs(
        np.abs(np.real(impz_dsp)[idx_res: idx_res + offset]) -
        np.abs(np.imag(impz_dsp)[idx_res: idx_res + offset])
    )
    )
        + idx_res]

    delta_f_u = np.abs(f5 - f6)
    return f_res / delta_f_u, delta_f_u
    #print('Qu: ', Qu, '∆fu: ', delta_f_u, 'MHz')


def get_loaded_q(sparam):
    # Calculating the loaded Q
    f_res, _ = sparam.find_reflection_resonance()
    cplx_dsp = sparam.get_dsp()
    f1 = sparam.freqs[np.argmax(np.imag(cplx_dsp))]
    f2 = sparam.freqs[np.argmin(np.imag(cplx_dsp))]
    delta_f_l = np.abs(f1 - f2)
    return f_res / delta_f_l, delta_f_l

    #print('Ql: ', Ql, '∆fl: ', delta_f_l, 'MHz')


def calc_ext_q(Qu, Ql):
    # Calculate the external Q from the other two Qs
    return 1 / (1 / Ql - 1 / Qu)
    # print('Qext_calc: ', Qext_calc)


def calc_beta(Qu, Qext):
    # Calculate the coupling factor
    return Qu / Qext
    # print('beta_calc', beta_calc)


def get_info_txt(Qu, delta_f_u, Ql, delta_f_l, f_res, beta_calc):
    return 'Qu={:0.0f}, ∆fu={:0.3f} [MHz]\nQl={:0.0f}, ∆fl={:0.3f} [MHz]\nf_res={:0.3f} [MHz], beta={:.3f}\n'.format(
        Qu, delta_f_u, Ql, delta_f_l, f_res, beta_calc)

# def get_ext_q(Ql, Qu):
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
