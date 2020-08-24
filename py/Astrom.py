from __future__ import print_function
import numpy as np
import sys

sys.path.insert(1, "../../pythonscripts/")

import pixel2radec as p2rd


class AstromEuclid:
    def __init__(self, dithers, CCD=None, CCD_ref=None):
        self.dithers = dithers

        if CCD != None and CCD_ref == None:
            self.input = "CCD"
            self.CCD, self.CCD_data = dithers[CCD], dithers[CCD].Data
            self.CCD_PV = dithers[CCD].Data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
            self.input_info = 'Taking AstromEuclid(dithers, CCD=..) as input'

        elif CCD_ref != None and CCD != None and dithers[CCD_ref].Data.Detector.DetectorId != dithers[
            CCD].Data.Detector.DetectorId:
            self.CCD_ref, self.CCD_ref_data = dithers[CCD_ref], dithers[CCD_ref].Data
            self.CCD_ref_PV = self.CCD_ref_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
            self.CCD, self.CCD_data = dithers[CCD], dithers[CCD].Data
            self.CCD_PV = self.CCD_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
            self.input_info = 'Error DetectorId diverge, you are not looking at the same DetectorIds \n Looking at different CCDs in 1DITHER'

        elif CCD_ref != None and CCD != None and dithers[CCD_ref].Data.ObservationDateTime.UTCObservationDateTime != \
                dithers[CCD].Data.ObservationDateTime.UTCObservationDateTime:
            self.input = "CCD&CCD_ref"
            self.CCD_ref, self.CCD_ref_data = dithers[CCD_ref], dithers[CCD_ref].Data
            self.CCD_ref_PV = self.CCD_ref_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
            self.CCD, self.CCD_data = dithers[CCD], dithers[CCD].Data
            self.CCD_PV = self.CCD_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
            self.input_info = 'Error ObservationDates diverge, you are not looking at the same UTCObservationDateTime \n Looking at 1CCD over different DITHER observations'

        else:
            try:
                self.CCD = self.CCD_ref = self.dithers[0]
                self.CCD_data = self.CCD_ref_data = self.dithers[0].Data
                self.CCD_PV = self.CCD_ref_PV = self.dithers[0].Data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
                self.input_info = 'No input given for CCD_ref or CCD \n --> Taking First idx of <dithers> is taken as input for both CCD as CCD_ref'
            except TypeError as err:
                self.CCD = self.CCD_ref = self.dithers
                self.CCD_data = self.CCD_ref_data = self.dithers.Data
                self.CCD_PV = self.CCD_ref_PV = self.dithers.Data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
                self.input_info = 'Handling TypeError:', err, '\n --> Taking single exposure from single dither.'

    def get_info(self, idx=0):
        if idx == 0:
            info = self.dithers
        else:
            info = self.dithers[idx]
        print("PipelineRun: {:50} ,\n Object ID: {:20},\n UTCObservationDateTime: {:20},\n DetectorId: {:20}".format(
            info.Header.PipelineRun, info.Header.object_id,
            (info.Data.ObservationDateTime.UTCObservationDateTime).strftime("%m-%d-%Y-%H-%M-%S"),
            info.Data.Detector.DetectorId))

    def calc_radec(self, X, Y, method="pvtan", exp_idx=None):
        if method == "pvtan_forframeplotting":
            if exp_idx == None:
                exp_data = self.CCD_data
                exp_data_pv = self.CCD_PV
            else:
                exp_data = self.dithers[exp_idx].Data
                exp_data_pv = exp_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
            RA, DEC = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, X, Y)
        if method == "pvtan":
            print("Converting input pixel coordinates to RA and DEC coordinates")
            RA, DEC = p2rd.pix2radec_pvtan(self.CCD_data, self.CCD_PV, X, Y)
        if method == "cdtan":
            RA, DEC = p2rd.pix2radec_cdtan(self.CCD_data, X, Y)
        if method == "cd":
            RA, DEC = p2rd.pix2radec_cd(self.CCD_data, X, Y)
        return RA, DEC

    def calc_sep_radec(self, p1, p2, method="pvtan"):  # http://spiff.rit.edu/classes/phys373/lectures/astrom/astrom.html
        RA_p1, DEC_p1 = self.calc_radec(p1[0], p1[1], method)
        RA_p2, DEC_p2 = self.calc_radec(p2[0], p2[1], method)
        delta_ra, delta_dec = RA_p1 - RA_p2, DEC_p1 - DEC_p2
        delta_RA = delta_ra * (
                    DEC_p1 + DEC_p2) / 2  # project the spherical sky onto a flat plane ---> degree to arscec ? --->    *3600[arcsec/degree]
        delta_DEC = delta_dec  # no correction is needed
        distance_p1p2 = np.sqrt((delta_RA) ** 2 + (delta_DEC) ** 2)
        return distance_p1p2, delta_RA, delta_DEC

    def calc_radec_center_ccd(self, exp_idx=None):  # idx=None): # return RA, DEC from center pixel at [1024, 2050]
        if exp_idx == None:
            exp_data = self.CCD_data
            exp_data_pv = self.CCD_PV
        else:
            exp_data = self.dithers[exp_idx].Data
            exp_data_pv = exp_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        c_ra, c_dec = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 1024, 2050)
        return c_ra, c_dec

    def calc_radec_corners_ccd(self, exp_idx=None):  # idx=None):
        if exp_idx == None:
            exp_data = self.CCD_data
            exp_data_pv = self.CCD_PV
        else:
            exp_data = self.dithers[exp_idx].Data
            exp_data_pv = exp_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        ul = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 0, 4100)
        ur = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 2048, 4100)
        ll = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 0, 0)
        lr = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 2048, 0)
        return ul, ur, ll, lr


class AstromAW:
    def __init__(self, dithers, CCD=None, CCD_ref=None):
        self.dithers = dithers

        if CCD != None and CCD_ref != None:
            self.input = "CCD"
            self.CCD_data, self.CCD_ref_data = dithers[CCD], dithers[CCD_ref]
            self.input_info = 'Taking AstromAW(dithers, CCD=..) as input'
        else:
            try:
                self.CCD_data = self.CCD_ref_data = self.dithers[0]
                self.input_info = 'No input given for CCD_ref or CCD \n --> Taking First idx of <dithers> is taken as input for both CCD as CCD_ref'
            except TypeError as err:
                self.CCD_data = self.CCD_ref_data = self.dithers
                self.input_info = 'Handling TypeError:', err, '\n --> Taking single exposure from single dither.'

        # if CCD != None and CCD_ref == None:
        #     self.input = "CCD"
        #     self.CCD, self.CCD_data = dithers[CCD], dithers[CCD].Data
        #     self.CCD_PV = dithers[CCD].Data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        #     self.input_info = 'Taking AstromEuclid(dithers, CCD=..) as input'
        #
        # elif CCD_ref != None and CCD != None and dithers[CCD_ref].Data.Detector.DetectorId != dithers[
        #     CCD].Data.Detector.DetectorId:
        #     self.CCD_ref, self.CCD_ref_data = dithers[CCD_ref], dithers[CCD_ref].Data
        #     self.CCD_ref_PV = self.CCD_ref_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        #     self.CCD, self.CCD_data = dithers[CCD], dithers[CCD].Data
        #     self.CCD_PV = self.CCD_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        #     self.input_info = 'Error DetectorId diverge, you are not looking at the same DetectorIds \n Looking at different CCDs in 1DITHER'
        #
        # elif CCD_ref != None and CCD != None and dithers[CCD_ref].Data.ObservationDateTime.UTCObservationDateTime != \
        #         dithers[CCD].Data.ObservationDateTime.UTCObservationDateTime:
        #     self.input = "CCD&CCD_ref"
        #     self.CCD_ref, self.CCD_ref_data = dithers[CCD_ref], dithers[CCD_ref].Data
        #     self.CCD_ref_PV = self.CCD_ref_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        #     self.CCD, self.CCD_data = dithers[CCD], dithers[CCD].Data
        #     self.CCD_PV = self.CCD_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
        #     self.input_info = 'Error ObservationDates diverge, you are not looking at the same UTCObservationDateTime \n Looking at 1CCD over different DITHER observations'
        #


    def get_info(self, idx=0):
        if idx == 0:
            info = self.dithers
        else:
            info = self.dithers[idx]
        print("Observing_block name: {:50} ,\n Observingblock ID: {:20},\n Chip name: {:20}".format( info.observing_block.name, observing_block.id, info.chip.name))

    def calc_radec(self, X, Y, method="pvtan", exp_idx=None):
        if method == "pvtan_forframeplotting":
            if exp_idx == None:
                exp_data = self.CCD_data
                exp_data_pv = self.CCD_data
            else:
                exp_data = self.dithers[exp_idx]
                exp_data_pv = self.dithers[exp_idx]
            RA, DEC = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, X, Y)
        if method == "pvtan":
            print("Converting input pixel coordinates to RA and DEC coordinates")
            RA, DEC = p2rd.pix2radec_pvtan(self.CCD_data, self.CCD_data, X, Y)
        if method == "cdtan":
            RA, DEC = p2rd.pix2radec_cdtan(self.CCD_data, X, Y)
        if method == "cd":
            RA, DEC = p2rd.pix2radec_cd(self.CCD_data, X, Y)
        return RA, DEC

    # def calc_sep_radec(self, p1, p2, method="pvtan"):  # http://spiff.rit.edu/classes/phys373/lectures/astrom/astrom.html
    #     RA_p1, DEC_p1 = self.calc_radec(p1[0], p1[1], method)
    #     RA_p2, DEC_p2 = self.calc_radec(p2[0], p2[1], method)
    #     delta_ra, delta_dec = RA_p1 - RA_p2, DEC_p1 - DEC_p2
    #     delta_RA = delta_ra * (
    #                 DEC_p1 + DEC_p2) / 2  # project the spherical sky onto a flat plane ---> degree to arscec ? --->    *3600[arcsec/degree]
    #     delta_DEC = delta_dec  # no correction is needed
    #     distance_p1p2 = np.sqrt((delta_RA) ** 2 + (delta_DEC) ** 2)
    #     return distance_p1p2, delta_RA, delta_DEC
    #
    # def calc_radec_center_ccd(self, exp_idx=None):  # idx=None): # return RA, DEC from center pixel at [1024, 2050]
    #     if exp_idx == None:
    #         exp_data = self.CCD_data
    #         exp_data_pv = self.CCD_PV
    #     else:
    #         exp_data = self.dithers[exp_idx].Data
    #         exp_data_pv = exp_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
    #     c_ra, c_dec = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 1024, 2050)
    #     return c_ra, c_dec
    #
    # def calc_radec_corners_ccd(self, exp_idx=None):  # idx=None):
    #     if exp_idx == None:
    #         exp_data = self.CCD_data
    #         exp_data_pv = self.CCD_PV
    #     else:
    #         exp_data = self.dithers[exp_idx].Data
    #         exp_data_pv = exp_data.NonLinearCoeffs.NonLinearTPVAstromCoeffs
    #     ul = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 0, 4100)
    #     ur = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 2048, 4100)
    #     ll = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 0, 0)
    #     lr = p2rd.pix2radec_pvtan(exp_data, exp_data_pv, 2048, 0)
    #     return ul, ur, ll, lr



# -----------------------------------------------------------------------------------------------
# DEFINITION USED TO BE ABLE TO IMPORT CLASS IN THE NOTEBOOKS
# This is done, in order to only change one file, which can be used in every Jupyter notebook ...
# -----------------------------------------------------------------------------------------------