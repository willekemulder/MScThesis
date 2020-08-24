
# transform_CV & pix2radec_cv from http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/wirwolf/docs/CD_PV_keywords.pdf
# transform_cvtan & transform_pvtan & pix2radec_pvtan from     # https://lost-contact.mit.edu/afs//mpe.mpg.de/i386_linux26/scisoft7.4/lib/python2.5/site-packages/stscidocs/docs/pytools/pytools.wcsutil-pysrc.html  WCSObject.xy2rd   'imgtools.xy2rd'
# --> Transforming from pixel coordinates to celestial coordinates [ra,dec]

import numpy as np


def transform_cv(data_CD, X, Y):
    ra = data_CD.CRVAL1 + data_CD.CD1_1 * (X - data_CD.CRPIX1) + data_CD.CD1_2 * (Y - data_CD.CRPIX2)
    dec = data_CD.CRVAL2 + data_CD.CD2_1 * (X - data_CD.CRPIX1) + data_CD.CD2_2 * (Y - data_CD.CRPIX2)
    return ra, dec

def transform_cvtan(data_CD, X, Y):
    x = data_CD.CD1_1 * (X - data_CD.CRPIX1) + data_CD.CD1_2 * (Y - data_CD.CRPIX2)
    y = data_CD.CD2_1 * (X - data_CD.CRPIX1) + data_CD.CD2_2 * (Y - data_CD.CRPIX2)
    xi, eta = np.deg2rad(x), np.deg2rad(y) # we need radians to work with goniometry
    alpha0 = np.deg2rad(data_CD.CRVAL1) # /180)*pi # convert to units of unit circle to use numpy modules
    delta0 = np.deg2rad(data_CD.CRVAL2) # /180)*pi # convert to units of unit circle to use numpy modules
    ra = np.arctan((xi / (np.cos(delta0)-eta*np.sin(delta0)))) + alpha0
    dec = np.arctan( ((eta*np.cos(delta0)+np.sin(delta0)) / (np.sqrt((np.cos(delta0)-eta*np.sin(delta0))**2 + xi**2))) )
    ra, dec=  np.rad2deg(ra), np.rad2deg(dec) # from radians back to degree
    return ra, dec

def transform_pvtan(data_CD, data_PV, X, Y):
    x = data_CD.CD1_1 * (X - data_CD.CRPIX1) + data_CD.CD1_2 * (Y - data_CD.CRPIX2)
    y = data_CD.CD2_1 * (X - data_CD.CRPIX1) + data_CD.CD2_2 * (Y - data_CD.CRPIX2)
    r = np.sqrt(x**2 + y**2)
    xi = data_PV.PV1_0 + (data_PV.PV1_1 * x) + (data_PV.PV1_2 * y) + (data_PV.PV1_3 * r) + (data_PV.PV1_4 * x**2) + (data_PV.PV1_5 * x * y) + (data_PV.PV1_6 * y**2) + (data_PV.PV1_7 * x**3) + (data_PV.PV1_8 * x**2*y) + (data_PV.PV1_9 * x*y**2) + (data_PV.PV1_10 * y**3)
    eta = data_PV.PV2_0 + (data_PV.PV2_1 * y) + (data_PV.PV2_2 * x) + (data_PV.PV2_3 * r) + (data_PV.PV2_4 * y**2) + (data_PV.PV2_5 * x * y) + (data_PV.PV2_6 * x**2) + (data_PV.PV2_7 * y**3) + (data_PV.PV2_8 * x*y**2) + (data_PV.PV2_9 * x**2*y) + (data_PV.PV2_10 * x**3)
    xi, eta = np.deg2rad(xi), np.deg2rad(eta) # we need radians to work with goniometry
    alpha0 = np.deg2rad(data_CD.CRVAL1) # /180)*pi # convert to units of unit circle to use numpy modules
    delta0 = np.deg2rad(data_CD.CRVAL2) # /180)*pi # convert to units of unit circle to use numpy modules
    ra = np.arctan((xi / (np.cos(delta0)-eta*np.sin(delta0)))) + alpha0
    dec = np.arctan( ((eta*np.cos(delta0)+np.sin(delta0)) / (np.sqrt((np.cos(delta0)-eta*np.sin(delta0))**2 + xi**2))) )
    ra, dec=  np.rad2deg(ra), np.rad2deg(dec) # from radians back to degree
    return ra, dec


def pix2radec_cd(data_CD, X, Y):

    try:
        data_CD.CD1_1
    except AttributeError:
        raise AttributeError \
            ("Check your data_CD input containing CD elements for using the pix2radec_pvtan function.")

    RA = data_CD.CRVAL1 + data_CD.CD1_1 * (X - data_CD.CRPIX1) + data_CD.CD1_2 * (Y - data_CD.CRPIX2)
    DEC = data_CD.CRVAL2 + data_CD.CD2_1 * (X - data_CD.CRPIX1) + data_CD.CD2_2 * (Y - data_CD.CRPIX2)
    return RA, DEC


def pix2radec_cdtan(data_CD, X, Y):

    try:
        data_CD.CD1_1
    except AttributeError:
        raise AttributeError \
            ("Check your data_CD input containing CD elements for using the pix2radec_pvtan function.")

    x = data_CD.CD1_1 * (X - data_CD.CRPIX1) + data_CD.CD1_2 * (Y - data_CD.CRPIX2)
    y = data_CD.CD2_1 * (X - data_CD.CRPIX1) + data_CD.CD2_2 * (Y - data_CD.CRPIX2)
    xi, eta = np.deg2rad(x), np.deg2rad(y)  # we need radians to work with goniometry
    alpha0 = np.deg2rad(data_CD.CRVAL1)  # /180)*pi # convert to units of unit circle to use numpy modules
    delta0 = np.deg2rad(data_CD.CRVAL2)  # /180)*pi # convert to units of unit circle to use numpy modules
    ra = np.arctan((xi / (np.cos(delta0) - eta * np.sin(delta0)))) + alpha0
    dec = np.arctan(((eta * np.cos(delta0) + np.sin(delta0)) / (np.sqrt((np.cos(delta0) - eta * np.sin(delta0)) ** 2 + xi ** 2))))
    RA, DEC = np.rad2deg(ra), np.rad2deg(dec)  # from radians back to degree
    return RA, DEC


def pix2radec_pvtan(data_CD, data_PV, X, Y):

    try:
        data_CD.CD1_1
    except AttributeError:
        raise AttributeError\
        ("Check your data_CD input containing CD elements for using the pix2radec_pvtan function.")
    try:
        data_PV.PV1_0
    except AttributeError:
        raise AttributeError\
        ("Check your data_PV input containing PV elements for using the pix2radec_pvtan function.")

    x = data_CD.CD1_1 * (X - data_CD.CRPIX1) + data_CD.CD1_2 * (Y - data_CD.CRPIX2)
    y = data_CD.CD2_1 * (X - data_CD.CRPIX1) + data_CD.CD2_2 * (Y - data_CD.CRPIX2)
    r = np.sqrt(x ** 2 + y ** 2)
    xi = data_PV.PV1_0 + (data_PV.PV1_1 * x) + (data_PV.PV1_2 * y) + (data_PV.PV1_3 * r) + (data_PV.PV1_4 * x ** 2) + (data_PV.PV1_5 * x * y) + (data_PV.PV1_6 * y ** 2) + (data_PV.PV1_7 * x ** 3) + (data_PV.PV1_8 * x ** 2 * y) + (data_PV.PV1_9 * x * y ** 2) + (data_PV.PV1_10 * y ** 3)
    eta = data_PV.PV2_0 + (data_PV.PV2_1 * y) + (data_PV.PV2_2 * x) + (data_PV.PV2_3 * r) + (data_PV.PV2_4 * y ** 2) + (data_PV.PV2_5 * x * y) + (data_PV.PV2_6 * x ** 2) + (data_PV.PV2_7 * y ** 3) + (data_PV.PV2_8 * x * y ** 2) + (data_PV.PV2_9 * x ** 2 * y) + (data_PV.PV2_10 * x ** 3)
    xi, eta = np.deg2rad(xi), np.deg2rad(eta)  # we need radians to work with goniometry
    alpha0 = np.deg2rad(data_CD.CRVAL1)  # /180)*pi # convert to units of unit circle to use numpy modules
    delta0 = np.deg2rad(data_CD.CRVAL2)  # /180)*pi # convert to units of unit circle to use numpy modules
    ra = np.arctan((xi / (np.cos(delta0) - eta * np.sin(delta0)))) + alpha0
    dec = np.arctan(((eta * np.cos(delta0) + np.sin(delta0)) / (np.sqrt((np.cos(delta0) - eta * np.sin(delta0)) ** 2 + xi ** 2))))
    RA, DEC = np.rad2deg(ra), np.rad2deg(dec)  # from radians back to degree
    return RA, DEC







# ----------------------------------------------------------------------------------
#            Calculating RA and DEC coordinates from pixel information
# ----------------------------------------------------------------------------------

#def calc_radec(data_CD, X, Y, data_PV=None, method="pvtan"):
#    print("Converting input pixel coordinates to RA and DEC coordinates")
#    if method == "pvtan":
#        RA, DEC = pix2radec_pvtan(data_CD, data_PV, X, Y)
#    if method == "cdtan":
#        RA, DEC = pix2radec_cdtan(data_CD, X, Y)
#    if method == "cd":
#       RA, DEC = pix2radec_cd(data_CD, X, Y)
#    return RA, DEC

