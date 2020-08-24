from __future__ import print_function
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.cm as cmx

import sys
import os

import pixel2radec as p2rd
import select_from_awedb as slct
from Astrom import AstromEuclid


# -----------------------------------------------------------------------------------------------
# DEFINITION USED TO BE ABLE TO SET LABELS FOR PLOTS
# This definition is made by Peyton Murray, to fix Matplotslib's Scientific Notation from  ...
# https://peytondmurray.github.io/coding/fixing-matplotlibs-scientific-notation/
# -----------------------------------------------------------------------------------------------


def label_offset(ax, axis="y"):
    if axis == "y":
        fmt = ax.yaxis.get_major_formatter()
        ax.yaxis.offsetText.set_visible(False)
        set_label = ax.set_ylabel
        label = ax.get_ylabel()

    elif axis == "x":
        fmt = ax.xaxis.get_major_formatter()
        ax.xaxis.offsetText.set_visible(False)
        set_label = ax.set_xlabel
        label = ax.get_xlabel()

    def update_label(event_axes):
        offset = fmt.get_offset()
        if offset == '':
            set_label("{}".format(label))
        else:
            set_label(u"{} ({})".format(label, offset))
        return

    ax.callbacks.connect("ylim_changed", update_label)
    ax.callbacks.connect("xlim_changed", update_label)
    ax.figure.canvas.draw()
    update_label(None)
    return


# -----------------------------------------------------------------------------------------------
# DEFINITIONS USED TO BE ABLE TO VISUALIZE KIDS SIMULATION DATA COMING FROM THE EUCLID DATABASE
# It allows to show dither information, as well as the infuence of astrometrical solutions
# -----------------------------------------------------------------------------------------------








