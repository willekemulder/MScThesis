from __future__ import print_function

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
display(HTML("<style> div.prompt {display:true} </style>"))
display(HTML("<style>.output_png {display: table-cell; text-align: center; vertical-align: middle;} </style>"))

try:
    from common.database.Context import context
    context.set_project('EUCLID')
    context.set_privileges(1)
    import pyfits as pfits
except ImportError as err:
    print("ImportError: Module Context not imported, go to AW environment")
    pass
except NameError as err:
    print("NameError: Module Context not imported, go to AW environment")
    pass


# ---------> In notebook <---------
#gaiasl_stage2 = (SourceList.SLID == 38382691)[0]   # stage2_gaia_v0.2
#gaiasl_stage2.info()
#gaiasl_sc3 = (SourceList.SLID == 40569211)[0]   # GAIA_SC3_pilot

import datetime
import numpy as np
import matplotlib.pyplot as plt

import sys
import os


from Astrom import AstromEuclid

import pixel2radec as p2rd
import select_from_awedb as slct
import plotutils_awe as pltut
