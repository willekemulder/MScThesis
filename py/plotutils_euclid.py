from __future__ import print_function
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as colors
import matplotlib.cm as cmx

import astropy.units as u
from astropy.coordinates import SkyCoord

# import pandas as pd

import sys
import os


import pixel2radec as p2rd
import select_from_dssdb as slct
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


def plot_dithers(Astrom_obj, dithersize=5, save=None):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    colors = ['b', 'g', 'r', 'm', 'y']

    for d in range(dithersize):
        data_cd = Astrom_obj.dithers[d].Data
        data_pv = data_cd.NonLinearCoeffs.NonLinearTPVAstromCoeffs

        ra_center, dec_center = p2rd.pix2radec_pvtan(data_cd, data_pv, 1024, 2050)
        ax.scatter(ra_center, dec_center, color=colors[d])

        # Calculating the coordinates in RA, DEC of the corners of the CCD
        ul = p2rd.pix2radec_pvtan(data_cd, data_pv, 0, 4100) # upperleft
        ur = p2rd.pix2radec_pvtan(data_cd, data_pv, 2048, 4100) # upperright
        ll = p2rd.pix2radec_pvtan(data_cd, data_pv, 0, 0) # lowerleft
        lr = p2rd.pix2radec_pvtan(data_cd, data_pv, 2048, 0) # lowerright

        # Defining the line segments to draw
        x_lineleft, y_lineleft = [ll[0], ul[0]], [ll[1], ul[1]]
        x_lineright, y_lineright = [lr[0], ur[0]], [lr[1], ur[1]]
        x_linedown, y_linedown = [ll[0], lr[0]], [ll[1], lr[1]]
        x_lineup, y_lineup = [ul[0], ur[0]], [ul[1], ur[1]]
        plt.plot(x_lineleft, y_lineleft, x_lineright, y_lineright, x_linedown, y_linedown, x_lineup, y_lineup,
                 color=colors[d])

    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('DEC [deg]')
    ax.set_title('Dither observation')

    if save != None:
        plt.savefig("./plots/dithers/{}.pdf".format(save), bbox_inches='tight', pad_inches=0)
    else:
        print("Plots are not saved!")
        plt.show()


def plot_dither_diff_dot(Astrom_obj, dithersize=5):
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    colors = ['b', 'g', 'r', 'm', 'y']

    for d in range(dithersize):
        data_cd = Astrom_obj.dithers[d].Data
        data_pv = data_cd.NonLinearCoeffs.NonLinearTPVAstromCoeffs

        ra_center, dec_center = p2rd.pix2radec_pvtan(data_cd, data_pv, 1024, 2050)
        ax.scatter(ra_center, dec_center, marker='x', color=colors[d], s=6)

        # Calculating the coordinates in RA, DEC of the corners of the CCD
        ul = p2rd.pix2radec_pvtan(data_cd, data_pv, 0, 4100) # upperleft
        ur = p2rd.pix2radec_pvtan(data_cd, data_pv, 2048, 4100) # upperright
        ll = p2rd.pix2radec_pvtan(data_cd, data_pv, 0, 0) # lowerleft
        lr = p2rd.pix2radec_pvtan(data_cd, data_pv, 2048, 0) # lowerright

        # Plotting the distances from the center of the CCD to its corners
        plt.plot([ul[0], ra_center], [ul[1], dec_center], color=colors[d], ls=':')
        plt.plot([ll[0], ra_center], [ll[1], dec_center], color=colors[d], ls='--')
        plt.plot([ra_center, lr[0]], [dec_center, lr[1]], color=colors[d], ls='-.')
        plt.plot([ra_center, ur[0]], [dec_center, ur[1]], color=colors[d], ls='-',
                 label=data_cd.ObservationDateTime.UTCObservationDateTime)

        # Calculating the distances from center to corner
        corners = [ul, ur, ll, lr]
        for corner in corners:
            d = np.sqrt((corner[0] - ra_center) ** 2 + (corner[1] - dec_center) ** 2)
            d = d * 3600  # 1 degree = 3600 arcsec
            plt.text(corner[0], corner[1], '{0:6.3f}'.format(d))

    dither = Astrom_obj.dithers[0]
    obstime = (dither.Data.ObservationDateTime.UTCObservationDateTime).strftime("%m-%d-%Y-%H-%M-%S")

    det_id = dither.Data.Detector.DetectorId
    name_run = dither.Header.PipelineRun
    obj_id = dither.Header.object_id

    plt.axis('equal')
    plt.legend(loc='lower left')

    ax.set_xlabel('RA [deg]')
    ax.set_ylabel('DEC [deg]')
    ax.set_title('Dither observation {}, {}'.format(obj_id, det_id))

    if not os.path.exists("./plots/astromsol_PV/dithers/{}".format(name_run)):
        print("Trying directory: {}".format("./plots/astromsol_PV/dithers/{}".format(name_run)), end='\n')
        print("New directory {} has been created".format("./plots/astromsol_PV/dithers/{}".format(name_run)))
        os.makedirs("./plots/astromsol_PV/dithers/{}".format(name_run))

    plt.savefig("./plots/astromsol_PV/dithers/{}/dithers_{}.pdf".format(name_run, det_id), bbox_inches='tight',
                pad_inches=0)
    plt.close()


def plot_oneDitherQuiver(Astrom_obj, exp_idxDQ_ref=0, exo_idxDQ=1):
    fig, ax = plt.subplots(figsize=(5, 5))

    # calculating the change of the central CCD pixels by canculating the shift of CCD centers in ra, dec coordinates
    ra_center, dec_center = Astrom_obj.calc_radec_center_ccd(exp_idx=exo_idxDQ)
    ra_center_ref, dec_center_ref = Astrom_obj.calc_radec_center_ccd(exp_idx=exp_idxDQ_ref)
    ra_shift, dec_shift = ra_center - ra_center_ref, dec_center - dec_center_ref

    # calculating the location of the corners of the CCD and reference CCD in ra, dec coordinates
    ul, ur, ll, lr = Astrom_obj.calc_radec_corners_ccd(exp_idx=exo_idxDQ)
    ul_ref, ur_ref, ll_ref, lr_ref = Astrom_obj.calc_radec_corners_ccd(exp_idx=exp_idxDQ_ref)

    corners_ref, corners = [ul_ref, ur_ref, ll_ref, lr_ref], [ul, ur, ll, lr]
    dra, ddec = [], []
    for idx in range(len(corners)):
        ra_corner_ref, ra_corner = corners_ref[idx][0], corners[idx][0]
        dec_corner_ref, dec_corner = corners_ref[idx][1], corners[idx][1]
        ra_diff_ref, dec_diff_ref = ra_corner_ref - ra_center_ref, dec_corner_ref - dec_center_ref
        ra_diff, dec_diff = ra_corner - ra_center, dec_corner - dec_center

        dra.extend([ra_diff - ra_diff_ref])
        ddec.extend([dec_diff - dec_diff_ref])

    ax.set_title(
        '{} {}'.format(Astrom_obj.CCD_data.Detector.DetectorId, Astrom_obj.CCD_data.ObservationDateTime.UTCObservationDateTime))

    q = ax.quiver([ra_center_ref, ra_center_ref, ra_center_ref, ra_center_ref],
                  [dec_center_ref, dec_center_ref, dec_center_ref, dec_center_ref], dra, ddec,
                  color=['b', 'g', 'k', 'r'], angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(q, X=0.3, Y=0.9, U=50, label='length = \n 0.18 arcseconds',
                 fontproperties={'weight': 'bold', 'size': 'smaller'}, labelpos='S')
    ax.set_aspect('equal')

    ax.set_xlabel(r'$RA$ (deg)')
    ax.set_ylabel(r'$DEC$ (deg)')

    ax.set_xlim(ra_center_ref - (np.max(dra) * 1.1), ra_center_ref + (np.max(dra) * 1.1))
    ax.set_ylim(dec_center_ref + (np.min(ddec) * 1.1), dec_center_ref - (np.max(ddec) * 1.1))
    plotting.label_offset(ax, "x")
    plotting.label_offset(ax, "y")
    ax.ticklabel_format(style='plain', axis='both', scilimits=(0, 0))

    plt.savefig("plots/astromsol_PV/quiver/diff{}_{}.pdf".format(Astrom_obj.CCD_data.Detector.DetectorId,
                                                                 Astrom_obj.CCD_data.ObservationDateTime.UTCObservationDateTime),
                bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_QuiverDist_Cntr2CrnrOrigin(Astrom_obj, exp_idxDQ_ref=0, dithersize=5):
    fig, axs = plt.subplots(2, 2, figsize=(12, 12), constrained_layout=True)
    layout = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]

    fig.suptitle(
        'Quiverplots Astrometrical Solutions, \n comparing distance in arcsec from central CCD to corner pixels per Exp \n Reference exposure: {}, DetectorId = {}'.format(
            Astrom_obj.dithers[exp_idxDQ_ref].Data.ObservationDateTime.UTCObservationDateTime,
            Astrom_obj.CCD_data.Detector.DetectorId), fontsize=16, y=1.09)

    dra, ddec, ddist = [], [], []

    ra_xaxis_low, ra_xaxis_high = [], []
    dec_xaxis_low, dec_xaxis_high = [], []

    ul_ref, ur_ref, ll_ref, lr_ref = Astrom_obj.calc_radec_corners_ccd(exp_idx=exp_idxDQ_ref)
    ra_center_ref, dec_center_ref = Astrom_obj.calc_radec_center_ccd(exp_idx=exp_idxDQ_ref)

    for exp_idxDQ in range(1, dithersize):

        dra_temp, ddec_temp, ddist_temp = [], [], []

        # calculating the change of the central CCD pixels by canculating the shift of CCD centers in ra, dec coordinates
        ra_center, dec_center = Astrom_obj.calc_radec_center_ccd(exp_idx=exp_idxDQ)

        # calculating the location of the corners of the CCD and reference CCD in ra, dec coordinates
        ul, ur, ll, lr = Astrom_obj.calc_radec_corners_ccd(exp_idx=exp_idxDQ)
        corners_ref, corners = [ul_ref, ur_ref, ll_ref, lr_ref], [ul, ur, ll, lr]

        for idx in range(len(corners)):
            ra_corner_ref, ra_corner = corners_ref[idx][0], corners[idx][0]
            dec_corner_ref, dec_corner = corners_ref[idx][1], corners[idx][1]
            ra_diff_ref, dec_diff_ref = ra_corner_ref - ra_center_ref, dec_corner_ref - dec_center_ref
            ra_diff, dec_diff = ra_corner - ra_center, dec_corner - dec_center

            dra_temp.append([(ra_diff - ra_diff_ref) * 3600])
            ddec_temp.append([(dec_diff - dec_diff_ref) * 3600])

        dra.append(dra_temp)
        ddec.append(ddec_temp)
        ddist.append(ddist_temp)

    ra_xaxis_high, ra_xaxis_low = max(max(dra))[0] * 1.05, min(min(dra))[0] * 1.05
    dec_xaxis_high, dec_xaxis_low = max(max(ddec))[0] * 1.05, min(min(ddec))[0] * 1.05
    center_plot = [[0, 0, 0, 0], [0, 0, 0, 0]]
    quiver_length = 50

    for idx_sub_plot, sub_plot in enumerate(layout):
        q = sub_plot.quiver(center_plot[0], center_plot[1], dra[idx_sub_plot], ddec[idx_sub_plot],
                            color=['b', 'g', 'k', 'r'], angles='xy', scale_units='xy', scale=1)
        sub_plot.quiverkey(q, X=0.3, Y=0.9, U=quiver_length,
                           label='length = \n {0:.3f} arcseconds'.format(quiver_length / 1000.0),
                           fontproperties={'weight': 'bold', 'size': 'smaller'}, labelpos='S')

        sub_plot.set_title(
            r'{}'.format(Astrom_obj.dithers[idx_sub_plot + 1].Data.ObservationDateTime.UTCObservationDateTime))
        sub_plot.set_xlabel(r'$RA$ (arcsec)')
        sub_plot.set_ylabel(r'$DEC$ (arcsec)')

        sub_plot.set_xlim(ra_xaxis_low, ra_xaxis_high)
        sub_plot.set_ylim(dec_xaxis_low, dec_xaxis_high)

    # Put a legend to the right of the last plotted axis
    legend_elements = [Line2D([0], [0], color='b', lw=4, label=r'$\Delta \;\;  D u,l$'),
                       Line2D([0], [0], color='g', lw=4, label=r'$\Delta \;\;  D u,r$'),
                       Line2D([0], [0], color='k', lw=4, label=r'$\Delta \;\;  D l,l$'),
                       Line2D([0], [0], color='r', lw=4, label=r'$\Delta \;\;  D l,r$')]
    plt.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.75, 0.85))

    plt.savefig(
        "plots/astromsol_PV/quiver/plot_QuiverDist_Cntr2Crnr/{}_{}.pdf".format(Astrom_obj.CCD_data.Detector.DetectorId,
                                                                               Astrom_obj.CCD_data.ObservationDateTime.UTCObservationDateTime),
        bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_dithers_difference(Astrom_obj, dithersize=5):
    fig, axs = plt.subplots(3, 1, figsize=(7, 10), constrained_layout=True)
    fig.suptitle('Astrometrical Solutions for {}'.format(Astrom_obj.CCD_data.Detector.DetectorId), fontsize=16, y=1.02)

    dra, ddec, ddist = [], [], []

    ul_ref, ur_ref, ll_ref, lr_ref = Astrom_obj.calc_radec_corners_ccd(exp_idx=0)
    ra_center_ref, dec_center_ref = Astrom_obj.calc_radec_center_ccd(exp_idx=0)

    for exo_idxDQ in range(dithersize):
        data_cd = Astrom_obj.dithers[d].Data

        # calculating the change of the central CCD pixels by canculating the shift of CCD centers in ra, dec coordinates
        ra_center, dec_center = Astrom_obj.calc_radec_center_ccd(exp_idx=exo_idxDQ)

        # calculating the location of the corners of the CCD and reference CCD in ra, dec coordinates
        ul, ur, ll, lr = Astrom_obj.calc_radec_corners_ccd(exp_idx=exo_idxDQ)

        corners_ref, corners = [ul_ref, ur_ref, ll_ref, lr_ref], [ul, ur, ll, lr]

        dra_temp, ddec_temp, ddist_temp = [], [], []
        corner_label = ['ul', 'ur', 'll', 'lr']

        for idx, corner in enumerate(corners):
            ra_corner_ref, ra_corner = corners_ref[idx][0], corner[0]
            dec_corner_ref, dec_corner = corners_ref[idx][1], corner[1]
            ra_diff_ref, dec_diff_ref = ra_corner_ref - ra_center_ref, dec_corner_ref - dec_center_ref
            ra_diff, dec_diff = ra_corner - ra_center, dec_corner - dec_center

            dra_temp.append([ra_diff - ra_diff_ref])
            ddec_temp.append([dec_diff - dec_diff_ref])

            # Calculating the distances from center to corner
            ddist_temp.append(np.sqrt( (corner[0] - ra_center) ** 2 + (corner[1] - dec_center) ** 2) * 3600)  # 1 degree = 3600 arcsec

        dra.append(dra_temp)
        ddec.append(ddec_temp)
        ddist.append(ddist_temp)
        print("{}/{}".format(exo_idxDQ + 1, dithersize), end='\r')

    dither = Astrom_obj.dithers[0]

    obstime_fancy = dither.Data.ObservationDateTime.UTCObservationDateTime
    obstime = (obstime_fancy).strftime("%m-%d-%Y-%H-%M-%S")

    det_id = dither.Data.Detector.DetectorId
    name_run = dither.Header.PipelineRun
    obj_id = dither.Header.object_id

    obstime = slct.get_UniqueObsTime(Astrom_obj.dithers)

    for idx_corner, corner in enumerate(corner_label):
        axs[0].plot(obstime[1:], [ra[idx_corner] for ra in dra][1:], '-')
        axs[1].plot(obstime[1:], [ra[idx_corner] for ra in ddec][1:], '--')
        axs[2].plot(obstime[1:], [ra[idx_corner] for ra in ddist][1:], ':')

    axs[0].set_title(r'$\Delta RA$ per ObsTime')
    axs[0].set_xlabel('Observation date')
    axs[0].set_ylabel(r'$\Delta RA$ (deg)')

    axs[1].set_title(r'$\Delta DEC$ per ObsTime')
    axs[1].set_xlabel('Observation date')
    axs[1].set_ylabel(r'$\Delta DEC$ (deg)')

    axs[2].set_title(r'$\Delta distance$ per ObsTime')
    axs[2].set_xlabel('Observation date')
    axs[2].set_ylabel(r'$\Delta distance$ (arcsec)')

    if not os.path.exists("./plots/astromsol_PV/PlotsGraphsAstroms/{}".format(name_run)):
        print("Trying directory: {}".format("./plots/astromsol_PV/PlotsGraphsAstroms/{}".format(name_run)),
              end='\n')
        print("New directory {} has been created".format(
            "./plots/astromsol_PV/PlotsGraphsAstroms/{}".format(name_run)))
        os.makedirs("./plots/astromsol_PV/PlotsGraphsAstroms/{}".format(name_run))

    plt.savefig("./plots/astromsol_PV/PlotsGraphsAstroms/{}/deltaAstromSol_{}.pdf".format(name_run, det_id),
                bbox_inches='tight', pad_inches=0)
    plt.close()


def plot_Survey_singleExp_PreDefinedFrame(data_exp_ref, data_exp, XpxSteps=8, YpxSteps=20, CCDcolumns=8, CCDrows=4, labels='off', info=None):
    fig, axs = plt.subplots(figsize=(20, 20), constrained_layout=True)

    column_idx = 0
    row_idx = 3

    qx, qy, qu, qv = [], [], [], []

    for ccd_idx, ccd_data in enumerate(data_exp_ref):
        print('Checking  ', ccd_data.Data.Detector.DetectorId,
              ', import at {:03.0f}%, data:'.format(((float(ccd_idx) / (CCDcolumns * CCDrows)) * 100)), ccd_data,
              end="\r")
        AstromObjRef = AstromEuclid(ccd_data)
        AstromObj = AstromEuclid(data_exp[ccd_idx])

        # calculating the location of the central pixels in ra, dec coordinates
        radecCenterRef, radecCenter = AstromObjRef.calc_radec_center_ccd(), AstromObj.calc_radec_center_ccd()
        shift = tuple(np.subtract(radecCenterRef, radecCenter))
        print(shift)

        Xgrid, Ygrid = 2144, 4200

        Y = np.linspace(50, 4150, num=YpxSteps + 1)
        X = np.linspace(48, 2000, num=XpxSteps + 1)

        for x in X:
            for y in Y:
                radecRef = AstromObjRef.calc_radec(x, y, method="pvtan_forframeplotting")
                radec = AstromObj.calc_radec(x, y, method="pvtan_forframeplotting")
                # diff = tuple(np.subtract(shift, tuple(np.subtract(radecRef, radec))))
                diff = tuple(np.subtract(radec, tuple(np.subtract(radecRef, shift))))

                d = tuple([3600 * i for i in diff])

                qx.append(x + column_idx * Xgrid)
                qy.append(y + row_idx * Ygrid)
                qu.append(d[0])
                qv.append(d[1])

        plt.axis('equal')
        if column_idx == CCDcolumns - 1:
            column_idx = 0
            row_idx -= 1
        else:
            column_idx += 1

    q = axs.quiver(qx, qy, qu, qv, angles='xy', scale_units='xy', scale=0.001, linewidths=0.5, headlength=2,
                   headwidth=2, width=0.001)
    qk = axs.quiverkey(q, 0.94, 1.035, (0.05 / np.sqrt(2)), r'length = $0.05$ arcsec', labelpos='S',
                       coordinates='figure')

    # plt.gca().invert_xaxis() # leave commented for correct plotting of NorthSouthEastWest directions
    plt.gca().invert_yaxis()

    if labels == 'off':
        axs.axes.get_xaxis().set_ticks([])
        axs.axes.get_yaxis().set_ticks([])
    elif labels == 'on':
        plt.text(radecCenterRef, radecCenter, '{}'.format(data_exp.Data.Detector.DetectorId),
                 horizontalalignment='center')
    else:
        print('Input for <labels> should be either on/off. Current input for labels=', labels)

    # comparing difference in ra dec for every pixels per Exp
    fig.suptitle('Quiverplot Astrometrical Solutions, \n Exposure = {}, Reference exposure: {} '.format(
        data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
        data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime), fontsize=24, y=1.05)
    plt.xlabel('pixels', fontsize=20)
    plt.ylabel('pixels', fontsize=20)

    if not os.path.exists("./plots/astromsol_PV/quiver/plot_Survey_singleExp_PreDefinedFrame/{}".format(name_run)):
        print("Trying directory: {}".format(
            "./plots/astromsol_PV/quiver/plot_Survey_singleExp_PreDefinedFrame/{}".format(name_run)), end='\n')
        print("New directory {} has been created".format(
            "./plots/astromsol_PV/quiver/plot_Survey_singleExp_PreDefinedFrame/{}".format(name_run)))
        os.makedirs("./plots/astromsol_PV/quiver/plot_Survey_singleExp_PreDefinedFrame/{}".format(name_run))

    if info == None:
        plt.savefig("./plots/astromsol_PV/quiver/plot_Survey_singleExp_PreDefinedFrame/survey_{}_{}.pdf".format(
            data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
            data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime), bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig("./plots/astromsol_PV/quiver/plot_Survey_singleExp_PreDefinedFrame/survey{}_{}_{}.pdf".format(info, data_exp[
            0].Data.ObservationDateTime.UTCObservationDateTime, data_exp_ref[
                                                                                                             0].Data.ObservationDateTime.UTCObservationDateTime),
                    bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_Survey_singleExp_DataFrame(data_exp_ref, data_exp, XpxSteps=8, YpxSteps=20, labels='off', info=None):
    fig, axs = plt.subplots(figsize=(20, 20), constrained_layout=True)

    qx, qy, qu, qv = [], [], [], []

    Xasis, Yaxis = 2048, 4100  # size in pixels Xgrid, Ygrid = 16000, 16000

    for ccd_idx, ccd_data_exp_ref in enumerate(data_exp_ref):
        print('Checking  ', ccd_data_exp_ref.Data.Detector.DetectorId,
              ', import at {:03.0f}%, data:'.format(((float(ccd_idx) / len(data_exp_ref)) * 100)), ccd_data_exp_ref,
              end="\r")
        AstromObjRef = AstromEuclid(ccd_data_exp_ref)
        AstromObj = AstromEuclid(data_exp[ccd_idx])
        CRpix1, CRpix2 = ccd_data_exp_ref.Data.CRPIX1, ccd_data_exp_ref.Data.CRPIX2

        # calculating the location of the central pixels in ra, dec coordinates to determine the shift
        radecCenterRef, radecCenter = AstromObjRef.calc_radec(0, 0, method="pvtan_forframeplotting"), AstromObj.calc_radec(0, 0, method="pvtan_forframeplotting")
        shift = tuple(np.subtract(radecCenterRef, radecCenter))
        X = np.linspace(CRpix1, CRpix1 - Xasis, num=XpxSteps + 1)
        Y = np.linspace(CRpix2, CRpix2 - Yaxis, num=YpxSteps + 1)

        for x in X:
            for y in Y:
                radecRef = AstromObjRef.calc_radec(x, y, method="pvtan_forframeplotting")
                radec = AstromObj.calc_radec(x, y, method="pvtan_forframeplotting")
                # diff = tuple(np.subtract(shift, tuple(np.subtract(radecRef, radec))))
                diff = tuple(np.subtract(radec, tuple(np.subtract(radecRef, shift))))
                d = tuple([3600 * i for i in diff])

                qx.append(x)
                qy.append(y)
                qu.append(d[0])
                qv.append(d[1])
        plt.axis('equal')

        if labels == 'off':
            axs.axes.get_xaxis().set_ticks([])
            axs.axes.get_yaxis().set_ticks([])
        elif labels == 'on':
            plt.text(np.median(X), np.median(Y), '{}'.format(data_exp[ccd_idx].Data.Detector.DetectorId),
                     horizontalalignment='center')
        else:
            print('Input for <labels> should be either on/off. Current input for labels=', labels)

        # print (data_exp[ccd_idx].Data.Detector.DetectorId, '   ra diff= ', d[0], '   dec diff= ', d[1], radec, radecRef, '   difference= ', diff)

    q = axs.quiver(qx, qy, qu, qv, angles='xy', scale_units='xy', scale=0.25, linewidths=0.5, headlength=2, headwidth=2,
                   width=0.001)
    qk = axs.quiverkey(q, 0.94, 1.035, (1.0 / np.sqrt(2)), r'length = $1.0$ arcsec', labelpos='S', coordinates='figure')

    # plt.gca().invert_xaxis() # leave commented for correct plotting of NorthSouthEastWest directions
    plt.gca().invert_yaxis()

    # comparing difference in ra dec for every pixels per Exp
    fig.suptitle('Quiverplot Astrometrical Solutions, \n Exposure = {}, Reference exposure: {} '.format(
        data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
        data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime), fontsize=24, y=1.05)
    plt.xlabel('pixels', fontsize=20)
    plt.ylabel('pixels', fontsize=20)

    if not os.path.exists("./plots/astromsol_PV/quiver/plot_Survey_singleExp_DataFrame/"):
        print("Trying directory: {}".format("./plots/astromsol_PV/quiver/plot_Survey_singleExp_DataFrame", end='\n'))
        print("New directory {} has been created".format("./plots/astromsol_PV/quiver/plot_Survey_singleExp_DataFrame"))
        os.makedirs("./plots/astromsol_PV/quiver/plot_Survey_singleExp_DataFrame")

    if info == None:
        plt.savefig("./plots/astromsol_PV/quiver/plot_Survey_singleExp_DataFrame/survey_{}_{}.pdf".format(
            data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
            data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime), bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig("./plots/astromsol_PV/quiver/plot_Survey_singleExp_DataFrame/survey{}_{}_{}.pdf".format( info,
            data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
            data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime), bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_CCD_diffExp_DataFrame(data_exp_ref, data_exp, XpxSteps=4, YpxSteps=10, labels='off', info=None):
    fig, axs = plt.subplots(figsize=(8, 8), constrained_layout=True)

    qx, qy, qu, qv = [], [], [], []

    Xasis, Yaxis = 2048, 4100  # size in pixels Xgrid, Ygrid = 16000, 16000

    AstromObjRef = AstromEuclid(data_exp_ref[0])
    AstromObj = AstromEuclid(data_exp[0])
    CRpix1, CRpix2 = data_exp_ref[0].Data.CRPIX1, data_exp_ref[0].Data.CRPIX2

    # calculating the location of the central pixels in ra, dec coordinates to determine the shift
    radecCenterRef, radecCenter = AstromObjRef.calc_radec(0, 0, method="pvtan_forframeplotting"), AstromObj.calc_radec(
        0, 0, method="pvtan_forframeplotting")
    shift = tuple(np.subtract(radecCenterRef, radecCenter))

    X = np.linspace(CRpix1, CRpix1 - Xasis, num=XpxSteps + 1)
    Y = np.linspace(CRpix2, CRpix2 - Yaxis, num=YpxSteps + 1)

    for x in X:
        for y in Y:
            radecRef = AstromObjRef.calc_radec(x, y, method="pvtan_forframeplotting")
            radec = AstromObj.calc_radec(x, y, method="pvtan_forframeplotting")
            diff = tuple(np.subtract(shift, tuple(np.subtract(radecRef, radec))))
            d = tuple([3600 * i for i in diff])
            qx.append(x)
            qy.append(y)
            qu.append(d[0])
            qv.append(d[1])
    plt.axis('equal')

    q = axs.quiver(qx, qy, qu, qv, angles='xy', scale_units='xy', scale=0.25, linewidths=0.5, headlength=2,
                   headwidth=2, width=0.001)
    qk = axs.quiverkey(q, 0.94, 1.035, (1.0 / np.sqrt(2)), r'length = $1.0$ arcsec', labelpos='S',
                       coordinates='figure')

    # plt.gca().invert_xaxis() # leave commented for correct plotting of NorthSouthEastWest directions
    plt.gca().invert_yaxis()

    if labels == 'off':
        axs.axes.get_xaxis().set_ticks([])
        axs.axes.get_yaxis().set_ticks([])
    elif labels == 'on':
        plt.text(radecCenterRef, radecCenter, '{}'.format(data_exp.Data.Detector.DetectorId), horizontalalignment='center')
    else:
        print('Input for <labels> should be either on/off. Current input for labels=', labels)


    # comparing difference in ra dec for every pixels per Exp
    fig.suptitle('Quiverplot Astrometrical Solutions, \n Exposure = {}, Reference exposure: {} of {}'.format(
        data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
        data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime,
        data_exp[0].Data.Detector.DetectorId), fontsize=20, y=1.05)
    plt.xlabel('pixels', fontsize=16)
    plt.ylabel('pixels', fontsize=16)

    if not os.path.exists("./plots/astromsol_PV/quiver/plot_CCD_diffExp_DataFrame/"):
        print("Trying directory: {}".format("./plots/astromsol_PV/quiver/plot_CCD_diffExp_DataFrame", end='\n'))
        print("New directory {} has been created".format("./plots/astromsol_PV/quiver/plot_CCD_diffExp_DataFrame"))
        os.makedirs("./plots/astromsol_PV/quiver/plot_CCD_diffExp_DataFrame")

    if info == None:
        plt.savefig("./plots/astromsol_PV/quiver/plot_CCD_diffExp_DataFrame/{}_{}_{}.pdf".format(
            data_exp[0].Data.Detector.DetectorId,
            data_exp[0].Data.ObservationDateTime.UTCObservationDateTime,
            data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime), bbox_inches='tight', pad_inches=0)
    else:
        plt.savefig("./plots/astromsol_PV/quiver/plot_CCD_diffExp_DataFrame/{}_{}_{}_{}.pdf".format(info, data_exp[0].Data.Detector.DetectorId, data_exp[0].Data.ObservationDateTime.UTCObservationDateTime, data_exp_ref[0].Data.ObservationDateTime.UTCObservationDateTime),
                    bbox_inches='tight', pad_inches=0)
    plt.show()
    return


class plot_DpdExtObsCol:
    def SourcesObsCol_XY(self, LDAC_OBJECTS, detector='all', labels='off', info=None):  # detector = 'ESO_CCD_#72' or
        fig = plt.figure(figsize=(4, 8))
        axs = fig.add_axes([0, 0, 1, 1])

        if detector == 'all':
            axs.set_title('Sources from DpdExtObsCollection for all detectors')
        elif detector != 'all':
            tabmask_det = DpdExtObsCol().SourceTableMaskDetID(LDAC_OBJECTS, detector)
            XWIN_IMAGE_DET = LDAC_OBJECTS['XWIN_IMAGE'][tabmask_det]
            YWIN_IMAGE_DET = LDAC_OBJECTS['YWIN_IMAGE'][tabmask_det]
            axs.set_title('Sources from DpdExtObsCollection, marked for {}'.format(detector))
            axs.scatter(XWIN_IMAGE_DET, YWIN_IMAGE_DET, color='r', marker='x', s=2, alpha=0.5)
        XWIN_IMAGE = LDAC_OBJECTS['XWIN_IMAGE']
        YWIN_IMAGE = LDAC_OBJECTS['YWIN_IMAGE']
        axs.scatter(XWIN_IMAGE, YWIN_IMAGE, color='k', s=0.1, alpha=0.3)
        plt.xlabel('XWIN_IMAGE')
        plt.ylabel('YWIN_IMAGE')
        plt.gca().invert_yaxis()

        if labels == 'off':
            axs.axes.get_xaxis().set_ticks([])
            axs.axes.get_yaxis().set_ticks([])
        elif labels == 'on':
            pass
        else:
            print('Input for <labels> should be either on/off. Current input for labels=', labels)

        if not os.path.exists("./plots/DpdExtObsCol/sources/SourcesObsCol_XY"):
            print("Trying directory: {}".format("./plots/DpdExtObsCol/sources/SourcesObsCol_XY"), end='\n')
            print("New directory {} has been created".format("./plots/DpdExtObsCol/sources/SourcesObsCol_XY"))
            os.makedirs("./plots/DpdExtObsCol/sources/SourcesObsCol_XY")

        if info == None:
            plt.savefig("./plots/DpdExtObsCol/sources/SourcesObsCol_XY/LDAC_OBJECTS_{}".format(detector),
                        bbox_inches='tight', pad_inches=0)
        else:
            plt.savefig("./plots/DpdExtObsCol/sources/SourcesObsCol_XY/LDAC_OBJECTS_{}_{}".format(info, detector),
                        bbox_inches='tight', pad_inches=0)
        plt.show()
        plt.close()

    def SourcesObsCol_RD(self, LDAC_OBJECTS_exp1, LDAC_OBJECTS_exp2, DpdEAS1, DpdEAS2, detector='all', info=None):
        # detector = 'ESO_CCD_#72'

        LDAC_OBJECTS = [LDAC_OBJECTS_exp1, LDAC_OBJECTS_exp2]
        DpdExtAsSols = [DpdEAS1, DpdEAS2]

        fig = plt.figure(figsize=(8, 16))
        axs = fig.add_axes([0, 0, 1, 1])

        color = ['k', 'b']
        sourcemarker = ['+', 'x']
        sourcecolor = ['g', 'r']

        for idx_assol, assol in enumerate(DpdExtAsSols):

            AstromObj = AstromEuclid(assol[0])

            datetime = assol[0].Data.ObservationDateTime.UTCObservationDateTime

            CRpix1, CRpix2 = assol[0].Data.CRPIX1, assol[0].Data.CRPIX2
            X_WIN = CRpix1 - LDAC_OBJECTS[idx_assol]["XWIN_IMAGE"]
            Y_WIN = CRpix2 - LDAC_OBJECTS[idx_assol]["YWIN_IMAGE"]

            for idx_src, source in enumerate(LDAC_OBJECTS[idx_assol]):
                radec = AstromObj.calc_radec(X_WIN[idx_src], Y_WIN[idx_src], method="pvtan_forframeplotting")
                axs.scatter(radec[0], radec[1], marker=sourcemarker[idx_assol], color=sourcecolor[idx_assol], s=5)

            ul = AstromObj.calc_radec(CRpix1, CRpix2, method="pvtan_forframeplotting")
            lr = AstromObj.calc_radec(CRpix1 - 2048, CRpix2 - 4100, method="pvtan_forframeplotting")
            ll = AstromObj.calc_radec(CRpix1, CRpix2 - 4100, method="pvtan_forframeplotting")
            ur = AstromObj.calc_radec(CRpix1 - 2048, CRpix2, method="pvtan_forframeplotting")

            plt.plot([ul[0], ur[0]], [ul[1], ur[1]], ls='-.', color=color[idx_assol], label=datetime)
            plt.plot([ll[0], ul[0]], [ll[1], ul[1]], ls='-.', color=color[idx_assol])
            plt.plot([ur[0], lr[0]], [ur[1], lr[1]], ls='-.', color=color[idx_assol])
            plt.plot([ll[0], lr[0]], [ll[1], lr[1]], ls='-.', color=color[idx_assol])

        axs.set_title('ObsCollection Sources for {}'.format(detector))
        axs.set_xlabel('RA (deg)')
        axs.set_ylabel('DEC (deg)')

        plt.gca().invert_yaxis()
        plt.legend(loc='center left', bbox_to_anchor=(0.90, 0.80))

        obstime = (datetime).strftime("%m-%d-%Y-%H-%M-%S")

        if not os.path.exists("./plots/DpdExtObsCol/sources/SourcesObsCol_RD/{}".format(obstime)):
            print("Trying directory: ./plots/DpdExtObsCol/sources/SourcesObsCol_RD/{}".format(obstime), end='\n')
            print("New directory ./plots/DpdExtObsCol/sources/SourcesObsCol_RD/{} has been created".format(obstime))
            os.makedirs("./plots/DpdExtObsCol/sources/SourcesObsCol_RD/{}".format(obstime))

        if info == None:
            plt.savefig("./plots/DpdExtObsCol/sources/SourcesObsCol_RD/{}/sources_{}.pdf".format(obstime, detector), bbox_inches='tight', pad_inches=0)
        else:
            plt.savefig("./plots/DpdExtObsCol/sources/SourcesObsCol_RD/{}/sources_{}_{}.pdf".format(obstime, info, detector), bbox_inches='tight', pad_inches=0)
        plt.show()
        plt.close()


