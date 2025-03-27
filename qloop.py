#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Loop to get Q values on many measurements

Xaratustrah 2019
"""
import sys
from sparameter import SParameter
from sparameter_tools import *
from plotters import *

with open('results.txt', 'a') as the_file:
    the_file.write('#name, Qu, delta_f_u, Ql, delta_f_l, f_res, beta_calc\n')
    for filename in sys.argv[1:]:
        print('Filename: \n', filename)
        sp = SParameter(filename)
        Qu, delta_f_u = get_unloaded_q(sp)
        Ql, delta_f_l = get_loaded_q(sp)
        Qext = calc_ext_q(Qu=Qu, Ql=Ql)
        f_res, _ = sp.find_reflection_resonance()
        beta_calc = calc_beta(Qu, Qext)
        txt = get_info_txt(Qu, delta_f_u, Ql, delta_f_l, f_res, beta_calc)
        print(txt)
        the_file.write(','.join(str(s)
                                for s in (Qu, delta_f_u, Ql, delta_f_l, f_res, beta_calc)) + '\n')
        print('Plotting to file {}.png\n'.format(sp.filename_wo_ext))
        plot_smith(sp, info_txt=txt)
