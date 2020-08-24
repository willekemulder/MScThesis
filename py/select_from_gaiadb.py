
from __future__ import print_function

from astroquery.gaia import Gaia
from astropy.io import fits, ascii

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np

import os



# ----------------------------------------------------------------------
#                          GAIA Database, TAP+
# ----------------------------------------------------------------------


class EuclidSky:
    def __init__(self):
        self.galactic_latitude = 30  # Galactic plane: |b| < 30 deg.
        self.celestial_latitude = 15  # ecliptic plane: |beta| < 15 deg.

        gal_l = np.linspace(-180, 180, 10000)

        self.coords_north = SkyCoord(l=gal_l * u.deg, b=self.galactic_latitude * u.deg, frame='galactic')
        self.coords_south = SkyCoord(l=gal_l * u.deg, b=-1 * self.galactic_latitude * u.deg, frame='galactic')
        self.icrs_coords_north = self.coords_north.icrs
        self.icrs_coords_south = self.coords_south.icrs
        # self.mean_icrs_coords_south = np.mean(self.coords_south.icrs)

        self.KIDS_S1_ra = [330.0, 360.0];
        self.KIDS_S_dec = [-35.5, 26.6]
        self.KIDS_S2_ra = [0.0, 52.5];
        self.KIDS_S_dec = [-35.5, 26.6]
        self.KIDS_N1_ra = [155.5, 225.5];
        self.KIDS_N1_dec = [-4.0, 4.0]
        self.KIDS_N2_ra = [225.5, 238.5];
        self.KIDS_N2_dec = [-2.0, 3.0]
        self.KIDS_NW2_ra = [128.5, 141.5];
        self.KIDS_NW2_dec = [-2.0, 3.0]
        self.KIDS_ND2_ra = [149.5, 150.5];
        self.KIDS_ND2_dec = [-1.7, 2.7]

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

    def kids_patch_ra(self):
        S1_ra = np.linspace(self.KIDS_S_ra[0], 360.0, num=10, endpoint=True)
        S2_ra = np.linspace(0.0, self.KIDS_S_ra[1], num=10, endpoint=True)
        N1_ra = np.linspace(self.KIDS_N1_ra[0], self.KIDS_N1_ra[1], num=7, endpoint=True)
        N2_ra = np.linspace(self.KIDS_N2_ra[0], self.KIDS_N2_ra[1], num=1, endpoint=True)
        NW2_ra = np.linspace(self.KIDS_NW2_ra[0], self.KIDS_NW2_ra[1], num=2, endpoint=True)
        ND2_ra = np.linspace(self.KIDS_ND2_ra[0], self.KIDS_ND2_ra[1], num=1, endpoint=True)
        return [S1_ra, S2_ra], N1_ra, N2_ra.NW2_ra, ND2_ra

    def kids_patch_dec(self):
        S_dec = np.linspace(self.KIDS_S_dec[0], self.KIDS_S_dec[1], num=10, endpoint=True)
        N1_dec = np.linspace(self.KIDS_N1_dec[0], self.KIDS_N1_dec[1], num=7, endpoint=True)
        N2_dec = np.linspace(self.KIDS_N2_dec[0], self.KIDS_N2_dec[1], num=1, endpoint=True)
        NW2_dec = np.linspace(self.KIDS_NW2_dec[0], self.KIDS_NW2_dec[1], num=2, endpoint=True)
        ND2_dec = np.linspace(self.KIDS_ND2_dec[0], self.KIDS_ND2_dec[1], num=1, endpoint=True)
        return S_dec, N1_dec, N2_dec.NW2_dec, ND2_dec


class Select_DATA_from_Gaia:
    def __init__(self):
        self.patch_filename = "gaia_euclid_patch"
        self.attributes = "solution_id, designation, source_id, random_index, ref_epoch, ra, ra_error, dec, dec_error, parallax, parallax_error, parallax_over_error, pmra, pmra_error, pmdec, pmdec_error, ra_dec_corr, ra_parallax_corr, ra_pmra_corr, ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, astrometric_n_obs_al, astrometric_n_obs_ac, astrometric_n_good_obs_al, astrometric_n_bad_obs_al, astrometric_gof_al, astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, astrometric_params_solved, astrometric_primary_flag, astrometric_weight_al, astrometric_pseudo_colour, astrometric_pseudo_colour_error, mean_varpi_factor_al, astrometric_matched_observations, visibility_periods_used, astrometric_sigma5d_max, frame_rotator_object_type, matched_observations, duplicated_source, phot_g_n_obs, phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_flux_over_error, phot_g_mean_mag, phot_bp_n_obs, phot_bp_mean_flux, phot_bp_mean_flux_error, phot_bp_mean_flux_over_error, phot_bp_mean_mag, phot_rp_n_obs, phot_rp_mean_flux, phot_rp_mean_flux_error, phot_rp_mean_flux_over_error, phot_rp_mean_mag, phot_bp_rp_excess_factor, phot_proc_mode, bp_rp, bp_g, g_rp, radial_velocity, radial_velocity_error, phot_variable_flag, l, b, ecl_lon, ecl_lat, datalink_url"

        self.galactic_latitude = 30  # Galactic plane: |b| < 30 deg.
        self.celestial_latitude = 15  # ecliptic plane: |beta| < 15 deg.

        if not os.path.exists("../../data/gaia/"):
            print("Trying directory: ../../data/gaia/", end='\n')
            print("New directory ../../data/gaia/")
            os.makedirs("../../data/gaia/")

    def query_gaia_limEUCLIDSKY(self):
        query = "select * from gaiadr2.gaia_source where abs(B)>{} and abs(ecl_lat)>{}".format(self.galactic_latitude, self.celestial_latitude)
        print("Start query_gaia_limEUCLIDSKY()")

        if not os.path.exists("../../data/gaia_euclid.csv"):
            filename = "gaia_euclid"
        else:
            print("PATH ../../data/gaia_euclid.csv already exist \n Overwrite existing file? [y/n]")
            answer = input()
            if answer == "y" or answer == "yes":
                filename = "gaia_euclid"
            elif answer == "n" or answer == "no":
                print("Give another filename (no gaia_euclid) to save the data")
                filename = input()
            else:
                print("Your answer: {} is not valid: Please return y, yes, n or no".format(answer))
                return query_gaia_EUCLIDSKY()
        job = Gaia.launch_job_async(query)
        r = job.get_results()
        ascii.write(r, "../../data/gaia/{}.csv".format(filename), delimiter=',')
        return None

    def query_gaia_EUCLIDSKY_north_pm(self):
        print('Start query_gaia_EUCLIDSKY_north_pm()')
        euclid_patch_ = EuclidSky()
        ra_north = euclid_patch_.north_patch_ra()
        dec_north = euclid_patch_.north_patch_dec()


        for idx_ra, ra in enumerate(ra_north[:-1]):
            for idx_dec, dec in enumerate(dec_north[:-1]):
                query_north = "select {} from gaiadr2.gaia_source where ra<{} and ra>{} and dec<{} and dec>{} and abs(B)>{} and abs(ecl_lat)>{} and parallax IS NOT NULL".format( self.attributes, ra_north[idx_ra + 1], ra, dec_north[idx_dec + 1], dec, self.galactic_latitude, self.celestial_latitude)
                print(query_north)
                job = Gaia.launch_job_async(query_north)
                r = job.get_results()
                file_ra = (ra + ra_north[idx_ra + 1]) / 2
                file_dec = (dec + dec_north[idx_dec + 1]) / 2
                ascii.write(r, "../../data/gaia/{}_{:.3f}_{:.3f}.csv".format(self.patch_filename, file_ra, file_dec),
                            delimiter=',')
        print('Queries for the north part of the Euclid sky survey  are done')
        return None

    def query_gaia_EUCLIDSKY_southleft_pm(self):
        print('Start query_gaia_EUCLIDSKY_southleft_pm()')
        euclid_patch_ = EuclidSky()
        ra_south_left = euclid_patch_.south_patch_left_ra()
        dec_south = euclid_patch_.south_patch_dec()

        for idx_ra, ra in enumerate(ra_south_left[:-1]):
            for idx_dec, dec in enumerate(dec_south[:-1]):
                query_south_left = "select {} from gaiadr2.gaia_source where ra<{} and ra>{} and dec<{} and dec>{} and abs(B)>{} and abs(ecl_lat)>{} and parallax IS NOT NULL".format( self.attributes, ra_south_left[idx_ra + 1], ra, dec_south[idx_dec + 1], dec, self.galactic_latitude, self.celestial_latitude)
                job = Gaia.launch_job_async(query_south_left)
                r = job.get_results()
                file_ra = (ra + ra_south_left[idx_ra + 1]) / 2
                file_dec = (dec + dec_south[idx_dec + 1]) / 2
                ascii.write(r, "../../data/gaia/{}_{:.3f}_{:.3f}.csv".format(self.patch_filename, file_ra, file_dec),
                            delimiter=',')
        print('Queries for the left south part of the Euclid sky survey are done')
        return None

    def query_gaia_EUCLIDSKY_southright_pm(self):
        print('Start query_gaia_EUCLIDSKY_southright_pm()')
        euclid_patch_ = EuclidSky()
        ra_south_right = euclid_patch_.south_patch_right_ra()
        dec_south = euclid_patch_.south_patch_dec()

        for idx_ra, ra in enumerate(ra_south_right[:-1]):
            for idx_dec, dec in enumerate(dec_south[:-1]):
                query_south_right = "select {} from gaiadr2.gaia_source where ra<{} and ra>{} and dec<{} and dec>{} and abs(B)>{} and abs(ecl_lat)>{} and parallax IS NOT NULL".format(self.attributes, ra_south_right[idx_ra + 1], ra, dec_south[idx_dec + 1], dec, self.galactic_latitude, self.celestial_latitude)
                job = Gaia.launch_job_async(query_south_right)
                r = job.get_results()
                file_ra = (ra + ra_south_right[idx_ra + 1]) / 2
                file_dec = (dec + dec_south[idx_dec + 1]) / 2
                ascii.write(r, "../../data/gaia/{}_{:.3f}_{:.3f}.csv".format(self.patch_filename, file_ra, file_dec),
                            delimiter=',')
        print('Queries for the right south part of the Euclid sky survey are done')
        return None

    def query_gaia_EUCLIDSKY_all_pm(self):
        print("Start query_gaia_EUCLIDSKY_all_pm()")
        self.query_gaia_EUCLIDSKY_north_pm()
        self.query_gaia_EUCLIDSKY_southleft_pm()
        self.query_gaia_EUCLIDSKY_southright_pm()
        return None

    def query_gaia_KiDS_testpatch(self, pointing): #pointing = [central_ra, central_dec]
        patch_ra = [pointing[0] - 0.5, pointing[0] + 0.5];
        patch_dec = [pointing[1] - 0.5, pointing[1] + 0.5]
        query_KiDS_testpatch = "select * from gaiadr2.gaia_source where ra<{} and ra>{} and dec<{} and dec>{} and parallax IS NOT NULL".format(patch_ra[1], patch_ra[0], patch_dec[1], patch_dec[0])
        print("Start query_gaia_KiDS_testpatch()")

        if not os.path.exists("../../data/gaia/gaia_KiDS_testpatch[{}_{}].csv".format(pointing[0], pointing[1])):
            filename = "gaia_KiDS_testpatch[{}_{}].csv".format(pointing[0], pointing[1])
        else:
            print("PATH ../../data/gaia_KiDS_testpatch[{}_{}].csv already exist \n Overwrite existing file? [y/n]").format(pointing[0], pointing[1])
            answer = input()
            if answer == "y" or answer == "yes":
                filename = "gaia_euclid"
            elif answer == "n" or answer == "no":
                print("Give another filename (no gaia_KiDS_testpatch) to save the data")
                filename = input()
            else:
                print("Your answer: {} is not valid: Please return y, yes, n or no".format(answer))
                return query_gaia_KiDS_testpatch()
        job = Gaia.launch_job_async(query_KiDS_testpatch)
        r = job.get_results()
        ascii.write(r, "../../data/gaia/{}".format(filename), delimiter=',')
        return None



def check_pc_import_gaiafromcsv():
    pc = input('On kapteyn or appeltje?')
    if pc == 'kapteyn':
        return os.chdir('/net/virgo01/data/users/mulder/MP/masterproject/notebooks/gaia')
    elif pc == 'appeltje':
        return os.chdir(
            '/Users/Willeke/Documents/Study_Astronomy/MSc_Astronomy/2019-2020/MP/masterproject/notebooks/gaia')
    else:
        print('Please import either kapteyn or appeltje to continue..')
        return check_pc_import_gaiafromcsv()


class Gaia_table:
    gaia_pm_columns = ['solution_id', 'source_id', 'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 'parallax_error', 'pmra', 'pmra_error', 'pmdec','pmdec_error']
    def __init__(self, pm_table, columns=gaia_pm_columns):
        if (type(pm_table) == type(pd.DataFrame())) == True:
            self.pm_table = pm_table
        else:
            print("No [pandas.core.frame.DataFrame] as input \n---> retrieving from local data can take a while..")
            check_pc_import_gaiafromcsv()
            cwd = os.getcwd()
            os.chdir("../../data/gaia")
            all_filenames = os.listdir('.')
            self.pm_table = pd.concat([pd.read_csv(f, delimiter=',',  usecols=columns) for f in all_filenames])
            os.chdir(cwd)



def select_AV(tb, av_value, av_attibute = "Dist", dev = 0.5):
    cwd = os.getcwd()
    check_pc_import_gaiafromcsv()

    try:
        len(av_value)
        for av in av_value:
            temp_table = tb.loc[ (tb[av_attibute] > av-dev) & (tb[av_attibute] < av+dev) ]
            temp_table.to_csv("../../data/gaia/allpatch/gaia_euclid_allpatch_pm_{}-{}.csv".format(av-dev, av+dev), sep=',', index=False)
            print('Table has been stored in ../../data/gaia/allpatch/gaia_euclid_allpatch_pm_{}-{}.csv'.format(av-dev, av+dev))
        os.chdir(cwd)
        return None
    except TypeError as err:
        print('Handling TypeError:', err, '\n --> Querying table for {} = {}-{} mas/year'.format(av_attibute, av_value-dev, av_value+dev))
        temp_table = tb.loc[(tb[av_attibute] > av_value-dev) & (tb[av_attibute] < av_value+dev)]
        temp_table.to_csv("../../data/gaia/allpatch/gaia_euclid_allpatch_pm_{}-{}.csv".format(av_value-dev, av_value+dev), sep=',', index=False)
        os.chdir(cwd)
        return 'Table has been stored in ../../data/gaia/allpatch/gaia_euclid_allpatch_pm_{}-{}.csv'.format(av_value-dev, av_value+dev)



def select_area(tb, areaname=['KIDS', 'KIDS_183.5_-2.5', 'DES', 'PSTARS']):
    print("Start select_area({})".format(areaname))
    if areaname == 'KIDS':
        S1_ra  = [330.0, 360.0]; S1_dec = [-35.5, -26.6]  #  RA range [330.0, 52.5]
        S2_ra = [0.0, 52.5]; S2_dec = [-35.5, -26.6]  #  RA range [330.0, 52.5]
        N1_ra = [155.5, 225.5]; N1_dec = [-4.0, 4.0]
        N2_ra = [225.5, 238.5]; N2_dec = [-2.0, 3.0]
        NW2_ra = [128.5, 141.5]; NW2_dec = [-2.0, 3.0]
        ND2_ra = [149.5, 150.5]; ND2_dec = [-1.7, 2.7]
        
        S1_area = tb[ (tb['ra'] >= np.float(S1_ra[0])) & (tb['ra'] < np.float(S1_ra[1])) & (tb['dec'] >= np.float(S1_dec[0])) & (tb['dec'] < np.float(S1_dec[1])) ]
        S2_area = tb[ (tb['ra'] >= np.float(S2_ra[0])) & (tb['ra'] < np.float(S2_ra[1])) & (tb['dec'] >= np.float(S2_dec[0])) & (tb['dec'] < np.float(S2_dec[1]))]
        print("{} sources in South KiDS area selected".format(len(S1_area)+len(S2_area)))
        N1_area = tb[ (tb['ra'] >= np.float(N1_ra[0])) & (tb['ra'] < np.float(N1_ra[1])) & (tb['dec'] >= np.float(N1_dec[0])) & (tb['dec'] < np.float(N1_dec[1])) ]
        print("{} sources in North1 KiDS area selected".format(len(N1_area)))
        N2_area = tb[ (tb['ra'] >= np.float(N2_ra[0])) & (tb['ra'] < np.float(N2_ra[1])) & (tb['dec'] >= np.float(N2_dec[0])) & (tb['dec'] < np.float(N2_dec[1])) ]
        print("{} sources in North2 KiDS area selected".format(len(N2_area)))
        NW2_area = tb[ (tb['ra'] >= np.float(NW2_ra[0])) & (tb['ra'] < np.float(NW2_ra[1])) & (tb['dec'] >= np.float(NW2_dec[0])) & (tb['dec'] < np.float(NW2_dec[1])) ]
        print("{} sources in NorthW2 KiDS area selected".format(len(NW2_area)))
        ND2_area = tb[ (tb['ra'] >= np.float(ND2_ra[0])) & (tb['ra'] < np.float(ND2_ra[1])) & (tb['dec'] >= np.float(ND2_dec[0])) & (tb['dec'] < np.float(ND2_dec[1])) ]
        print("{} sources in NorthD2 KiDS area selected".format(len(ND2_area)))
        
        KIDS_area = pd.DataFrame()
        for table in [S1_area, S2_area, N1_area, N2_area, NW2_area, ND2_area]:
            KIDS_area = KIDS_area.append(table, ignore_index=True)
        select = KIDS_area

    elif areaname == 'KIDS_183.5_-2.5':
        central_pointing = 'KIDS_183.5_-2.5'
        central_ra = 183.5; central_dec = -2.5
        patch_ra = [central_ra-0.5, central_ra+0.5] ; patch_dec = [central_dec-0.5, central_dec+0.5]
        patch_area = tb[ (tb['ra'] >= np.float(patch_ra[0])) & (tb['ra'] < np.float(patch_ra[1])) & (tb['dec'] >= np.float(patch_dec[0])) & (tb['dec'] < np.float(patch_dec[1])) ]
        print("{} sources in South KiDS area selected".format(len(patch_area)))
        select = patch_area


    elif areaname == 'PSTARS':
        pi_ra = [0.0, 360.0]; pi_dec = [-30.0, 0]
        pi_area = tb[(tb['ra'] >= np.float(pi_ra[0])) & (tb['ra'] < np.float(pi_ra[1])) & (tb['dec'] >= np.float(pi_dec[0])) & (tb['dec'] < np.float(pi_dec[1]))]
        print("{} sources in North1 KiDS area selected".format(len(pi_area)))

        PSTARS_area = pd.DataFrame()
        select = PSTARS_area

    else:
        try:
            central_ra, central_dec = areaname
            patch_ra = [central_ra - 0.5, central_ra + 0.5];
            patch_dec = [central_dec - 0.5, central_dec + 0.5]
            patch_area = tb[(tb['ra'] >= np.float(patch_ra[0])) & (tb['ra'] < np.float(patch_ra[1])) & (tb['dec'] >= np.float(patch_dec[0])) & (tb['dec'] < np.float(patch_dec[1]))]
            print("{} sources in South KiDS area selected".format(len(patch_area)))
            select = patch_area
        except:
            print('Enter a value for areaname=[KIDS, KIDS_183.5_-2.5, DES, PSTARS], or [RA, DEC] from central pointing')


    if len(select) == 0:
        print('Nothing found for the selected area ({}) in the input table'.format(areaname))
    else:
        return select
