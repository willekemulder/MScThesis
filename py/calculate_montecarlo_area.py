from __future__ import print_function

import astropy.units as u
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import colors
import matplotlib.cm as cmx

#from matplotlib.axes._axes import _log as matplotlib_axes_logger
#matplotlib_axes_logger.setLevel('ERROR')

import numpy as np
import pandas as pd

import time


class EuclidSky:
    def __init__(self):
        self.galactic_latitude = 30  # Galactic plane: |b| < 30 deg.
        self.celestial_latitude = 15  # ecliptic plane: |beta| < 15 deg.
        gal_l = np.linspace(-180, 180, 10000)

        sky_ra = np.linspace(0, 360, 360)
        sky_dec = np.linspace(-90, 90, 180)
        sky_ra_grid, sky_dec_grid = np.meshgrid(sky_ra, sky_dec)

        self.coords_north = SkyCoord(l=gal_l * u.deg, b=self.galactic_latitude * u.deg, frame='galactic')
        self.coords_south = SkyCoord(l=gal_l * u.deg, b=-1 * self.galactic_latitude * u.deg, frame='galactic')
        self.coords_sky = SkyCoord(ra= sky_ra_grid, dec= sky_dec_grid, frame='icrs', unit=u.deg)
        self.icrs_coords_north = self.coords_north.icrs
        self.icrs_coords_south = self.coords_south.icrs
        self.icrs_coords_sky = self.coords_sky.icrs

    def north_patch_ra(self):
        return np.linspace(min(self.icrs_coords_north.ra.deg), max(self.icrs_coords_north.ra.deg), num=10,
                           endpoint=True)

    def north_patch_dec(self):
        return np.linspace(min(self.icrs_coords_north.dec.deg), max(self.icrs_coords_north.dec.deg), num=10,
                           endpoint=True)

    def south_patch_right_ra(self):
        mean_icrs_ra = np.mean(self.icrs_coords_south.ra.deg)
        right_ra_lim = min(self.icrs_coords_south[self.icrs_coords_south.ra.deg > mean_icrs_ra].ra.deg)
        return np.linspace(right_ra_lim, max(self.icrs_coords_south.ra.deg), num=5, endpoint=True)

    def south_patch_left_ra(self):
        mean_icrs_ra = np.mean(self.icrs_coords_south.ra.deg)
        left_ra_lim = max(self.icrs_coords_south[self.icrs_coords_south.ra.deg < mean_icrs_ra].ra.deg)
        return np.linspace(min(self.icrs_coords_south.ra.deg), left_ra_lim, num=6, endpoint=True)

    def south_patch_dec(self):
        return np.linspace(min(self.icrs_coords_south.dec.deg), max(self.icrs_coords_south.dec.deg), num=10,
                           endpoint=True)

    def fill_entire_skyPerSdeg_deg(self):
        return self.icrs_coords_sky

    def fill_entire_skyPerSdeg_gal(self):
        return self.icrs_coords_sky.galactic

    def fill_entire_skyRandom_gal_ecl(self, n=10000):
        # generate values between 0 and 1 for convertion to random radians
        ra, dec = np.random.random(2*n).reshape(2, -1)

        # use the random numbers range [0.0;1.0] to determine random radians
        Random_ra_rad = (2*np.pi) * (ra - 0.5)
        Random_dec_rad = np.arcsin(2.*(dec-0.5))

        # transform radians to degree
        to_deg = 180/np.pi
        skyRandom_ra = (Random_ra_rad*to_deg)+180
        skyRandom_dec = Random_dec_rad*to_deg

        skyRandom_ra_grid, skyRandom_dec_grid = np.meshgrid(skyRandom_ra, skyRandom_dec)
        coords_skyRandom = SkyCoord(ra=skyRandom_ra_grid, dec=skyRandom_dec_grid, frame='icrs', unit=u.deg)
        icrs_coords_skyRandom = coords_skyRandom.icrs
        return icrs_coords_skyRandom.galactic, icrs_coords_skyRandom.barycentricmeanecliptic



def MC_calcarea(n_sample = 250):
    ''' Calculate the surface areas of irregular shapes on the sky that are too difficult to fit. Done using the Monte
    Carlo Method which will generate random dummy sources and checks in which of the irregular shapes they show up.

    (1) Generating a known number of points, at n_sample random declinations and n_sample of right ascension: Using
        EuclidSky().fill_entire_skyRandom_gal_ecl(n_sample), it generates n_sample^2 randomly located dummy sources.
    (2) Counting the number of dummy sources that lie inside the different irregular shapes: The shapes are set by
        exclusion zones for both galactic and ecliptic coordinates. The galactic exclusion zone is easily filtered out
        since we define this entire zone as a single 'irregular area'. The ecliptic exlusion zone affects the entire
        sky and therefore all the irregular areas. To take this area into account, we count and store the number of
        dummy sources lying in |beta| < 15 separately as 'ecliptic_exclusion'
    (3) Calculating the proportionality [sources in patch]/[total generated sources]
    (4) Return lists containing # dummy sources, proportionality per area and the sources excluded by the ecliptic.

        Parameters
        ----------
        n_sample : int [Default = 250]
            The number of random dummy sources you want to create [CAUTION: it will be n_sample^2]

        Returns
        -------
        sources
            List of lists; one for every irregular area with the number of dummy sources in its specific region

        percentage
            List of lists; one for every irregular area with [sources in patch]/[total generated sources] per region

        ecliptic_exclusion
            integer number revealing the number of sources that where generated lying in |beta| < 15

    ------------------------------- User example --------------------------------
    ########### Using MC_calcarea() to create calculate area ratios #############

        sources, percentage, ecliptic_exclusion = MC_calcarea(n_sample = 10)

    '''
    slice_options = [-90.0, -60.0, -45.0, -30.0, 0.0, 30.0, 45.0, 60.0, 90.0]

    arr_area = np.array(slice_options)

    # (1) Generating n_sample^2 random dummy sources using EuclidSky()
    E = EuclidSky()
    C_gal, C_ecl = E.fill_entire_skyRandom_gal_ecl(n_sample)

    areas = [[], [], [], [], [], [], [], []]

    ecliptic_exclusion = 0

    # (2) Counting the number of dummy sources that lie inside the different irregular shapes
    for idx_dummy, dummy in enumerate(C_gal.flatten()):
        b_coord, l_coord = dummy.b.deg, dummy.l.deg # Find galactic lat & long: b, l = Skycoord.b.deg, Skycoord.l.deg
        beta_coord = C_ecl.flatten()[idx_dummy].lat.deg # Find ecliptic latitude for the dummy: beta = Skycoord.lat.deg
        if abs(beta_coord) > 15.0:
            # |beta| > 15.0: we find the galactic lat to look to its boarders
            find_area = np.argwhere([arr_area <= b_coord][0] == True)
            idx_area = np.max(find_area)
            areas[idx_area].append(dummy)
        else:
            # |beta| < 15.0: we only count the number of dummies, since it lies in the ecliptic exclusion zone
            ecliptic_exclusion += 1

    # (3) Calculating the proportionality [sources in patch]/[total generated sources]
    sources = [float(len(areas[i])) for i in range(len(slice_options)-1)] # number of sources in patch
    percentage = [i / (float(len(C_gal))**2) for i in sources] # [sources in patch]/[total sources]
    # (4) Return lists containing # dummy sources, proportionality per area and the sources excluded by the ecliptic
    return sources, percentage, ecliptic_exclusion


def _use_array(iterations = None, n_sample = 250):
    ''' Obtain lists containing information about the surface areas of irregular shapes on the sky. Done using lists.
    (1) Set up two initial empty lists called MC_s and MC_p,
    (2) Creating n_sample^2 dummy sources and checks in which patch they are located and,
    (3) Return MC_s and MC_p containing the number of sources per patch and their ratios.

        Parameters
        ----------
        n_sample : int [Default = 250]
            The number of random dummy sources you want to create [CAUTION: it will be n_sample^2]

        iterations : int [Default = None]
            The amount of iterations you want to have for the creation of the random dummy sources

        Returns
        -------
        MC_s
            List containing counting of the number of dummy sources in specific regions

        MC_p
            List containing the ratio [dummy sources found in patch] / [n_sample^2 generated dummy sources]

    ------------------------------- User example --------------------------------
    ######### Using _use_array() to print table containing area ratios #########

        MC_s, MC_p = _use_array(iterations = 10, n_sample = 10) # sources, percentages

        header = "| {:^12s} | {:^12s} | {:^12s} | {:^12s} | {:^12s} | {:^12s} | {:^12s} | {:^12s} | {:^12s} |".format(
        "Galactic latitude", '-90$^\circ$ < b < -60$^\circ$', '-60$^\circ$ < b < -45$^\circ$',
        '-45$^\circ$ < b < -30$^\circ$', '-30$^\circ$ < b < 0$^\circ$', '0$^\circ$ < b < 30$^\circ$',
        '30$^\circ$ < b < 45$^\circ$', '45$^\circ$ < b < 60$^\circ$', '60$^\circ$ < b < 90$^\circ$')
        if print_table == "yes":
        print(header)
        print("|:-----------:|:------------:|:------------:|:------------:|:------------:|:------------:|:------------:|:-----------:|:------------:|")
        s_perarea = "| {:^12s} | {:5.0f} | {:5.0f} | {:5.0f} | {:5.0f} | {:5.0f} | {:5.0f} | {:5.0f} | {:5.0f} | ".format(
            "# dummy sources", sources[0], sources[1], sources[2], sources[3], sources[4], sources[5], sources[6],
            sources[7])
        s_percentage = "| {:^12s} | {:0.2f} | {:0.2f} | {:0.2f} | {:0.2f} | {:0.2f} | {:0.2f} | {:0.2f} | {:0.2f} | ".format(
            "ratio", percentage[0], percentage[1], percentage[2], percentage[3], percentage[4], percentage[5],
            percentage[6], percentage[7])
        print(s_perarea, '\n', s_percentage)

    '''
    # (1) Set up two initial lists
    MC_s = []
    MC_p = []
    if iterations == None:
        iterations = int(input('Number of iterations? ..\n'))

    for i in range(itterations):
        # (2) Creating n_sample^2 dummy sources and checks in which patch they are located
        sources, percentage, ecl_excl = MC_calcarea(n_sample)
        MC_s.append(sources)
        MC_p.append(percentage)
        print('{0:3.1f}/100.0  with {1:1.0f} sources in ecliptic exclusion zone'.format((float(i+1)/iterations)*100, ecl_excl), end="\r")
    # (3) Return the arrays containing the number of sources per patch and their ratios
    return MC_s, MC_p

def _use_pandas(iterations = None, n_sample = 250, get_excl = 'no'):
    ''' Obtain table containing the surface areas of irregular shapes on the sky. Done using Pythons pandas module.
    (1) Set up an initial DataFrame
    (2) Creating n_sample^2 dummy sources and checks in which patch they are located
    (3) Store the table as function Output and .csv file.

        Parameters
        ----------
        n_sample : int [Default = 250]
            The number of random dummy sources you want to create [CAUTION: it will be n_sample^2]

        iterations : int [Default = None]
            The amount of iterations you want to have for the creation of the random dummy sources

        get_excl : str [Default = 'no']
            Option to obtain list of length [iterations] showing the # sources in the ecliptic excl zone per iteration

        Returns
        -------
        df_table : pandas.DataFrame
            Contains the number of dummy sources and its ratio compaired to all sources for specific regions

        {time}_areacalc_i{iterations}n{n_sample}.csv:
            Stored table as a .csv file which can be opened by using e.g.:
            table = pd.read_csv("../data/gaia/area_calc/20200410_095611_areacalc_i100n300.csv")

        total_ecl_excl : List [optional]
            A list of length [iterations] showing the number of sources in the ecliptic exclusion zone per iteration

    ------------------------------- User example --------------------------------
    ####### Using _use_pandas() to create table containing area ratios #########

        table = _use_pandas(iterations = 10, n_sample = 10)


    ######## Using _use_pandas() for sanity check on area calculations #########

        n_sample = 10; iterations = 10
        table, total_ecl_excl = _use_pandas(iterations = 10, n_sample = n_sample, get_excl = 'yes')
        sCols = ['S_9060', 'S_6045', 'S_4530', 'S_3000', 'S_0030', 'S_3045', 'S_4560', 'S_6090'] # obtaining sources columns
        pCols = ['P_9060', 'P_6045', 'P_4530', 'P_3000', 'P_0030', 'P_3045', 'P_4560', 'P_6090'] # obtaining sources columns

        A_hemisphere = (2/np.pi) * 360 * 90 # calculating the surface area of the entire sky
        A_sky = 2 * A_hemisphere

        check = np.copy(iterations)
        for i in range(iterations):
            total_s_patch = 0
            areas_patch = []
            for idx_col, Col in enumerate(sCols):
                total_s_patch += table[Col][i]
                areas_patch.append(table[pCols[idx_col]][i] * A_sky)
            areas_patches = np.sum(areas_patch)
            total_s_patches = np.sum(total_s_patch)
            area_excl = areas_patches * total_ecl_excl[i] / total_s_patches
            if int(areas_patches+area_excl) == int(A_sky):
                print("SANITY CHECK COMPLETE, total area {}/{} as expected".format(i, iterations), end='\r')
            else:
                print("POSSIBLE WRONG CALCULATION, total area {}/{} not as expected".format(i, iterations), end='\r')
                check -= 1
        print("Correct area calculations for {}/{} iterations".format(check, iterations), end='\r')

    '''
    # (1) Set up an initial DataFrame, areas are defined as 9060 meaning -90 < b < -60 where 6090 means 60 < b < 90
    df_table = pd.DataFrame({'idx': pd.Series([], dtype='int'),
                        'S_9060': pd.Series([], dtype='int'),
                        'P_9060': pd.Series([], dtype='float'),
                        'S_6045': pd.Series([], dtype='int'),
                        'P_6045': pd.Series([], dtype='float'),
                        'S_4530': pd.Series([], dtype='int'),
                        'P_4530': pd.Series([], dtype='float'),
                        'S_3000': pd.Series([], dtype='int'),
                        'P_3000': pd.Series([], dtype='float'),
                        'S_0030': pd.Series([], dtype='int'),
                        'P_0030': pd.Series([], dtype='float'),
                        'S_3045': pd.Series([], dtype='int'),
                        'P_3045': pd.Series([], dtype='float'),
                        'S_4560': pd.Series([], dtype='int'),
                        'P_4560': pd.Series([], dtype='float'),
                        'S_6090': pd.Series([], dtype='int'),
                        'P_6090': pd.Series([], dtype='float') })
    if iterations == None:
        iterations = int(input('Number of iterations? ..\n'))

    total_ecl_excl = []
    for i in range(iterations):
        # (2) Creating n_sample^2 dummy sources and checks in which patch they are located
        sources, percentage, ecl_excl = MC_calcarea(n_sample)
        data = [{'idx': i,
                 'S_9060': sources[0], 'P_9060': percentage[0], 'S_6045': sources[1], 'P_6045': percentage[1],
                 'S_4530': sources[2], 'P_4530': percentage[2], 'S_3000': sources[3], 'P_3000': percentage[3],
                 'S_0030': sources[4], 'P_0030': percentage[4], 'S_3045': sources[5], 'P_3045': percentage[5],
                 'S_4560': sources[6], 'P_4560': percentage[6], 'S_6090': sources[7], 'P_6090': percentage[7] }]
        df_table = df_table.append(data, ignore_index=True, sort=False)
        count = (float(i+1)/iterations)*100
        print('{0:3.1f}/100.0  with {1:1.0f} sources in ecliptic exclusion zone'.format(count, ecl_excl), end="\r")
        total_ecl_excl.append(ecl_excl)
    # (3) Store the table by using the time of the saving as filename to prevent overwriting older tables
    timestr = time.strftime("%Y%m%d_%H%M%S")
    filename = "{}_areacalc_i{}n{}.csv".format(timestr, iterations, n_sample)
    df_table.to_csv(filename, index=False)
    if get_excl != 'no':
        return df_table, total_ecl_excl
    return df_table


table = _use_pandas(iterations = 1, n_sample = 3)
# ------------------- Load csv to obtain table [optional] ---------------------
# table = pd.read_csv("../data/gaia/area_calc/20200410_095611_areacalc_i100n300.csv")


# ------------------ Visualise calculation of areas by MC ---------------------
# ------------------ Visualise calculation of areas by MC ---------------------
def plot_surfacearea_latitude(table, slice_option='b_excl', n_sample=250, get_sigmas='no'):
    ''' Visualisation for the Monte Carlo method used to obtain the surface areas between the lines dividing galactic
    latitudes. It calculates the area in squared degrees, by taking the proportions obtained by '_use_pandas()' and
    divide it by the surface area of the entire sky. By repeating the latter, we call it an Monte Carlo method since
    more iterations will eventually lead to an accurate estimation of the 'real surface area'. Each time we calculate
    the average of the previous iterations together with the new one. The function returns a plot showing an estimate
    for the surface areas in units of degrees.

    (1) Retrieving iteration information from table
    (2) Calculating and plot the means of the ratio [sources in patch]/[total generated sources] for each single patch
    (3) Calculating and plot the surface area determined by multiplying ratio with total sky area for single patches
    (4) Saving, showing and closing a figure showing the surface area
    (5) [optional] Returning deviation [sigmas] in sources per iteration giving sigma for surface area calculation

        Parameters
        ----------
        table: pandas.DataFrame
            Data Frame obtained by use of _use_pandas()

        slice_option: str [Default = 'b_excl']
            Either 'b_excl', or 'b_incl', for excluding or including the galactic exclusion area |b| < 30

        n_sample : int [Default = 250]
            The number of random dummy sources you want to create [CAUTION: it will be n_sample^2]

        get_sigmas : str [Default = 'no']
            Option to obtain list of length [iterations] showing sigma for our calculation of the surface area


        Returns
        -------
        sigmas [optional]
            List of deviation in sources per iteration giving sigma for our calculation of the surface area

        MC_excl_i{iterations}n{n_sample}.pdf:
            Plot stored as a .pdf file which is stored in the directory "plots"

    ------------------------------- User example --------------------------------
    ###### Using plot_surfacearea_latitude() to create plot showing areas #######

        sigmas = plot_surfacearea_latitude(table, n_sample = 250)

    '''
    # predefined column names that we are using in our analysis
    if slice_option == 'b_excl':
        perc = ['P_9060', 'P_6045', 'P_4530', 'P_3045', 'P_4560',
                'P_6090']  # obtaining data for all b patches excluding the galactic plane area, like Euclid
        slices = ['$-90^{\circ} < b < -60^{\circ}$', '$-60^{\circ} < b <-45^{\circ}$',
                  '$-45^{\circ} < b < -30^{\circ}$', '   $30^{\circ} < b <$  $45^{\circ}$',
                  '   $45^{\circ} < b <$  $60^{\circ}$',
                  '   $60^{\circ} < b <$  $90^{\circ}$']  # set to perc_all if plot should contain 0 < |b| < 30
    else:
        perc = ['P_9060', 'P_6045', 'P_4530', 'P_3000', 'P_0030', 'P_3045', 'P_4560',
                'P_6090']  # obtaining data for all b patches from table
        slices = ['$-90^{\circ} < b < -60^{\circ}$', '$-60^{\circ} < b < -45^{\circ}$',
                  '$-45^{\circ} < b < -30^{\circ}$',
                  '$-30^{\circ} < b < 0^{\circ}$', '$0^{\circ} < b < 30^{\circ}$', '$30^{\circ} < b < 45^{\circ}$',
                  '$45^{\circ} < b < 60^{\circ}$',
                  '$60^{\circ} < b < 90^{\circ}$']  # set to slices_all if plot should contain 0 < |b| < 30

    # (1) Retrieving iteration information from table
    iterations = len(table['idx'])  # counting the amount of iterations that were used to create this table

    # setting up the colors used in the plot. Each area will have its own color, given by the colormap 'viridis'
    viridis = plt.get_cmap('viridis')
    cNorm = colors.Normalize(vmin=0, vmax=len(perc))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=viridis)
    colorVal = [scalarMap.to_rgba(idx_color) for idx_color in range(iterations)]

    # calculating the surface area of the entire sky, needed to obtain the surface area per patch
    A_hemisphere = (2 / np.pi) * 360 * 90
    A_sky = 2 * A_hemisphere

    # ------------------------------- Graph ----------------------------------
    # PLOT SHOWING THE SURFACE RATIO [PATCH/SKY] AS PATCHES ARE DEFINES BY LATITUDE COORDINATES IN slice_option

    # setting up the base of the figure using matplotlib
    fig, ax = plt.subplots(figsize=(8, 8), ncols=1, nrows=1)
    ax.set_title(r"Surface area (A) of galactic latitude patches (A$_b$) ", fontsize=16)
    ax.xaxis.set_label_text(r"iterations", fontsize=14)
    ax.yaxis.set_label_text(r"ratio A$_b$/A$_{sky}$", fontsize=14)

    # (2) Calculating and plot the means of the ratio [sources in patch]/[total generated sources] for single patches
    for idx_p, p in enumerate(perc):
        ratio = [table[p][:i].mean() for i in range(iterations)]
        ax.plot(range(iterations), ratio, color=colorVal[idx_p], linestyle='--', label="MC")

    # (3) Calculating and plot the surface area determined by multiplying ratio with total sky area for single patches
    ax1 = ax.twinx()  # instantiate a second axes that shares the same x-axis as the ratio plot
    sigmas = []
    for idx_p, p in enumerate(perc):
        A_patch = [(table[p][:i].mean()) * (A_sky) for i in range(iterations)]
        ax1.plot(range(iterations), A_patch, color=colorVal[idx_p], linestyle='--', label="MC")
        sigmas.append(np.std(np.array(A_patch)[~np.isnan(A_patch)]))  # calculating 1 sigma uncertainty for A_b means
    ax1.yaxis.labelpad = 20
    ax1.set_ylabel(r"A $[deg^2]$", fontsize=14, rotation=270)

    # defining and adding a legend to help the reader understand the differences in areas per galactic latitude range
    if len(perc) == 8:
        legend_elements = [Line2D([0], [0], color=colorVal[0], label=slices[0] + r"   $\sim$ {0:4.0f} deg$^2$".format(
            (table[perc[0]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[1], label=slices[1] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[1]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[2], label=slices[2] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[2]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[3], label=slices[3] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[3]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[4], label=slices[4] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[4]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[5], label=slices[5] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[5]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[6], label=slices[6] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[6]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[7], label=slices[7] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[7]][:].mean()) * (A_sky))), ]
    else:  # if  len(perc) != 8, it shoud be len(perc) == 6, as defined by slice_option
        legend_elements = [Line2D([0], [0], color=colorVal[0], label=slices[0] + r"   $\sim$ {0:4.0f} deg$^2$".format(
            (table[perc[0]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[1], label=slices[1] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[1]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[2], label=slices[2] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[2]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[3], label=slices[3] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[3]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[4], label=slices[4] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[4]][:].mean()) * (A_sky))),
                           Line2D([0], [0], color=colorVal[5], label=slices[5] + r"   $\sim$ {0:4.0f} deg$^2$".format(
                               (table[perc[5]][:].mean()) * (A_sky))), ]
    ax.legend(handles=legend_elements)

    # (4) Saving, showing and closing a figure showing the surface area
    plt.savefig("./plots/thesis_MC_excl_i{}n{}.pdf".format(iterations, n_sample))
    plt.show()
    plt.close(fig)
    if get_sigmas != 'no':
        # (5) [optinal] Returning deviation [sigmas] in sources per iteration giving sigma for our calculation of the surface area
        return sigmas
    return None



