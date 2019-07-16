# -*- coding: utf-8 -*-
"""
Plotters for different S parameters operations

"""
import matplotlib.pyplot as plt
from smithplot import SmithAxes


def plot_smith(sparam, info_txt=None):
    # Plot section
    plt.figure()  # figsize=(6, 6))
    ax = plt.subplot(1, 1, 1, projection='smith', grid_minor_enable=True)
    plt.plot(sparam.cspar, markevery=10, label='Measurement data',
             datatype=SmithAxes.S_PARAMETER)
    plt.plot(sparam.get_dsp(), markevery=10, label='Detuned short position',
             datatype=SmithAxes.S_PARAMETER)
    plt.plot(sparam.get_dop(), markevery=10, label='Detuned open position',
             datatype=SmithAxes.S_PARAMETER)
    plt.legend(loc="lower right", fontsize=8)
    plt.title(sparam.file_basename_wo_ext)
    plt.text(0, 0.6, info_txt, size=9, rotation=0,
             ha="left", va="top",
             bbox=dict(boxstyle="square",
                       ec=(1., 0.5, 0.5),
                       fc=(1., 0.8, 0.8),
                       )
             )

    plt.savefig("{}.png".format(sparam.filename_wo_ext),
                format="png",
                dpi=600,
                bbox_inches="tight",)

    plt.clf()
