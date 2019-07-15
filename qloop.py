#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Loop to get Q values on many measurements

Xaratustrah 2019
"""
import sys
from sparameter import SParameter

with open('results.txt', 'a') as the_file:
    the_file.write('#file, Qu, Ql, delta_f_u, f_res, beta\n')
    for filename in sys.argv[1:]:
        print('Filename: ', filename)
        sp = SParameter(filename)
        the_file.write(','.join(str(s)
                                for s in sp.get_q_values()) + '\n')
