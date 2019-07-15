#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Loop to get Q values on many measurements

Xaratustrah 2019
"""
import sys
from sparameter import SParameter

with open('results.txt', 'a') as the_file:
    the_file.write('#file, Qu, Ql, delta_f_u, f_res, beta_calc\n')
    for filename in sys.argv[1:]:
        print('Filename: ', filename)
        sp = SParameter(filename)
        Qu, Ql, delta_f_u, f_res, beta_calc = sp.get_q_values()
        sp.plot_smith(Qu, Ql, delta_f_u, f_res, beta_calc)
        the_file.write(','.join(str(s)
                                for s in (Qu, Ql, delta_f_u, f_res, beta_calc)) + '\n')
