from __future__ import print_function

from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
display(HTML("<style> div.prompt {display:true} </style>"))
display(HTML("<style>.output_png {display: table-cell; text-align: center; vertical-align: middle;} </style>"))

import warnings
warnings.filterwarnings('ignore')

import astroquery
from astroquery.gaia import Gaia
from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord

import pandas as pd
import requests


import datetime
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.gridspec import GridSpec


import sys
import os

from Astrom import AstromEuclid

import pixel2radec as p2rd
import select_from_gaiadb as slct
import plotutils_gaia as pltut

