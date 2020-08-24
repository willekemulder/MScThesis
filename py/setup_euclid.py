from __future__ import print_function

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
display(HTML("<style> div.prompt {display:true} </style>"))
display(HTML("<style>.output_png {display: table-cell; text-align: center; vertical-align: middle;} </style>"))


from euclid.config import startup
from euclid.config.Environment import Env
from euclid.main.aweimports import *
from common.database.Context import context
context.set_project('TEST');

from astropy.io import fits
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord

import datetime
import numpy as np
import matplotlib.pyplot as plt

import sys
import os



from Astrom import AstromEuclid
from select_from_dssdb import DpdExtObsCol, Get


import pixel2radec as p2rd
import select_from_dssdb as slct
import plotutils_euclid as pltut
