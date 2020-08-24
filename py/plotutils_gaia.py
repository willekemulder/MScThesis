from __future__ import print_function
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import astropy.units as u
from astropy.coordinates import SkyCoord


import sys
import os



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



class plot_Gaia:
    def __init__(self, pm_table):
        self.pm_table = pm_table

        if 'Dist' not in self.pm_table.columns:
            dist = np.sqrt(self.pm_table["pmdec"] ** 2 + self.pm_table["pmra"] ** 2)
            self.pm_table['Dist'] = dist
        else:
            print('Column [Dist] already present in the imported pandas table')



    def SourcesEuclidSky_PM_Hist(self, info=None):
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

        ax1.set_title("Proper motion [RA direction]")
        hist1 = self.pm_table["pmra"].plot.hist(bins=1000, alpha=0.5, ax=ax1)  # .plot(ax=ax1)
        ax1.yaxis.set_label_text("Angular Velocity[mas/year]")
        ax1.set_xlim((-100, 100))

        ax2.set_title("Proper motion [DEC direction]")
        hist2 = self.pm_table["pmdec"].plot.hist(bins=1000, alpha=0.5, ax=ax2)  # .plot(ax=ax2)
        ax2.yaxis.set_label_text("Angular Velocity[mas/year]")
        ax2.set_xlim((-100, 100))

        ax3.set_title("Proper motion")

        if 'Dist' not in self.pm_table.columns:
            dist = np.sqrt(self.pm_table["pmdec"] ** 2 + self.pm_table["pmra"] ** 2)
            self.pm_table['Dist'] = dist
        else:
            dist = self.pm_table['Dist']
            print('Column [Dist] already present in the imported pandas table')

        hist3 = (dist).plot.hist(bins=1000, alpha=0.5, ax=ax3)
        ax3.yaxis.set_label_text("Angular Velocity (mas/year)")
        ax3.set_xlim((-10, 150))

        fig.subplots_adjust(hspace=0.6)

        if not os.path.exists("./plots/SourcesEuclidSky_PM_Hist/"):
            print("Trying directory: {}".format("./plots/SourcesEuclidSky_PM_Hist", end='\n'))
            print("New directory {} has been created".format("./plots/SourcesEuclidSky_PM_Hist"))
            os.makedirs("./plots/SourcesEuclidSky_PM_Hist")

        if info == None:
            plt.savefig("./plots/SourcesEuclidSky_PM_Hist/EuclidSky_pmhist.pdf", bbox_inches='tight', pad_inches=0)
        else:
            plt.savefig("./plots/SourcesEuclidSky_PM_Hist/EuclidSky_pmhist-{}.pdf".format(info), bbox_inches='tight',
                        pad_inches=0)
        plt.show()
        plt.close()
        return None


    def SourcesEuclidSky_PM_Map(self, info=None):
        gal_l = np.linspace(-180, 180, 10000)

        coords = SkyCoord(l=gal_l * u.deg, b=30 * u.deg, frame='galactic')
        coords_l = SkyCoord(l=gal_l * u.deg, b=-30 * u.deg, frame='galactic')
        icrs_coords = coords.icrs
        icrs_coords_l = coords_l.icrs

        c_sources = SkyCoord(ra=10.68458 * u.degree, dec=41.26917 * u.degree, frame='icrs')

        fig, frame = plt.subplots(1, 1, figsize=(20, 12))
        frame.scatter(icrs_coords.ra, icrs_coords.dec)
        frame.scatter(icrs_coords_l.ra, icrs_coords_l.dec)

        jet = plt.get_cmap('jet')
        cNorm = colors.Normalize(vmin=0, vmax=len(range(80)) - 1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        for i in range(80):
            colorVal = scalarMap.to_rgba(i)
            dist_bin = self.pm_table[(self.pm_table['Dist'] > i) & (self.pm_table['Dist'] < i + 1)]
            if i in np.linspace(1, 79, 4):
                frame.scatter(dist_bin["ra"], dist_bin["dec"], color=colorVal, label='{}-{} [mas/year]'.format(i, i + 1),
                              rasterized=True)
            else:
                print('try')
                frame.scatter(dist_bin["ra"], dist_bin["dec"], color=colorVal, rasterized=True)
            print('{}, check'.format(i), end="\r")

        frame.set_title('Proper motion [Angular velocity] of Gaia sources in Euclid footprint', fontsize=20, y=1.05)
        plt.xlabel('RA (J2000)', fontsize=16)
        plt.ylabel('DEC (J2000)', fontsize=16)
        print('done labels')
        plt.gca().invert_xaxis()
        plt.legend(loc='center left', bbox_to_anchor=(0.90, 0.80))
        print('done legend')

        if not os.path.exists("./plots/SourcesEuclidSky_PM_Map/"):
            print("Trying directory: {}".format("./plots/SourcesEuclidSky_PM_Map", end='\n'))
            print("New directory {} has been created".format("./plots/SourcesEuclidSky_PM_Map"))
            os.makedirs("./plots/SourcesEuclidSky_PM_Map")

        if info == None:
            plt.savefig("./plots/SourcesEuclidSky_PM_Map/EuclidSky_ra_dec.pdf", bbox_inches='tight', pad_inches=0)
        else:
            print('try saving')
            plt.savefig("./plots/SourcesEuclidSky_PM_Map/EuclidSky_ra_dec-{}.pdf".format(info), bbox_inches='tight',
                        pad_inches=0)
        plt.show()
        plt.close()
        return None


    def SourcesEuclidSky_PMperAV_Map(self, info=None):
        gal_l = np.linspace(-180, 180, 10000)

        coords = SkyCoord(l=gal_l * u.deg, b=30 * u.deg, frame='galactic')
        coords_l = SkyCoord(l=gal_l * u.deg, b=-30 * u.deg, frame='galactic')
        icrs_coords = coords.icrs
        icrs_coords_l = coords_l.icrs

        c_sources = SkyCoord(ra=10.68458 * u.degree, dec=41.26917 * u.degree, frame='icrs')

        jet = plt.get_cmap('jet')
        cNorm = colors.Normalize(vmin=0, vmax=len(range(80)) - 1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

        if not os.path.exists("./plots/SourcesEuclidSky_PMperAV_Map/"):
            print("Trying directory: {}".format("./plots/SourcesEuclidSky_PMperAV_Map", end='\n'))
            print("New directory {} has been created".format("./plots/SourcesEuclidSky_PMperAV_Map"))
            os.makedirs("./plots/SourcesEuclidSky_PMperAV_Map")

        for i in range(80):
            fig, frame = plt.subplots(1, 1, figsize=(20, 12))
            frame.scatter(icrs_coords.ra, icrs_coords.dec, rasterized=True)
            frame.scatter(icrs_coords_l.ra, icrs_coords_l.dec, rasterized=True)

            frame.set_title('Proper motion [AV: {}-{} (mas/year)] of Gaia sources in Euclid footprint'.format(i, i + 1),
                            fontsize=20, y=1.05)
            plt.xlabel('RA (J2000)', fontsize=16)
            plt.ylabel('DEC (J2000)', fontsize=16)

            plt.gca().invert_xaxis()
            colorVal = scalarMap.to_rgba(i)
            dist_bin = self.pm_table[(self.pm_table['Dist'] > i) & (self.pm_table['Dist'] < i + 1)]
            frame.scatter(dist_bin["ra"], dist_bin["dec"], color=colorVal, rasterized=True)
            print('{}, check'.format(i), end="\r")
            if info == None:
                plt.savefig("./plots/SourcesEuclidSky_PMperAV_Map/EuclidSky_PMamplitude_{}-{}_ra_dec.pdf".format(i, i + 1),
                            bbox_inches='tight', pad_inches=0)
            else:
                plt.savefig(
                    "./plots/SourcesEuclidSky_PMperAV_Map/EuclidSky_PMamplitude_{}-{}_ra_dec-{}.pdf".format(i, i + 1, info),
                    bbox_inches='tight', pad_inches=0)
            plt.show()
            plt.close()
        return None

