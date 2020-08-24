
from __future__ import print_function

import sys
import os


# ----------------------------------------------------------------------
#  Selecting unique observationtimes and detectorIDs
# ----------------------------------------------------------------------

class Get:
    def __init__(self, query):
        self.query = query

    def UniqueDetectorId(self):
        DetectorIDS = []
        for q in self.query:
            ID = q.Data.Detector.DetectorId
            if ID not in DetectorIDS:
                DetectorIDS.append( ID )
            else:
                pass
        return sorted(DetectorIDS)

    def UniqueObsTime(self):
        ObsTime = []
        for q in self.query:
            OBS = q.Data.ObservationDateTime.UTCObservationDateTime
            if OBS not in ObsTime:
                ObsTime.append( OBS )
            else:
                pass
        return sorted(ObsTime)


# ----------------------------------------------------------------------
#  Selecting Sources from FITS files obtained through DpdExtObsCol
# ----------------------------------------------------------------------


class Select_FITS_from_DpdExtObsCol:
    def __init__(self):
        pass

    def XYlimDetId(self, query, Xaxis=2048, Yaxis=4100):
        detId = query[0].Data.Detector.DetectorId
        CRpix1 = query[0].Data.CRPIX1
        CRpix2 = query[0].Data.CRPIX2
        xlims = [CRpix1, CRpix1 - Xaxis]
        ylims = [CRpix2, CRpix2 - Yaxis]
        return xlims, ylims, detId

    def SourceTableMaskDetID(self, LDAC_OBJECTS, query):
        # Usage:
        # #1 mask_ccd72 = DpdExtObsCol().SourceTableMaskDetID(OC_ImageObjects, 'ESO_CCD_#72')
        # #2 mask_ccd72 = DpdExtObsCol().SourceTableMaskDetID(OC_ImageObjects, query_exp1ccd72_stage2)
        try:
            [xmax, xmin], [ymax, ymin], detId = self.XYlimDetId(query)
        except AttributeError as err:
            'Handling AttributeError:', err, '\n --> Input query =', query, ' was no query.'
            detId = query
        except TypeError as err:
            'Handling TypeError:', err, '\n --> Input query =', query, ' was no query.'
            detId = query
        mask_detId = (LDAC_OBJECTS['DETECTOR'] == detId)
        return mask_detId

    def SourceTableDetID(self, LDAC_OBJECTS, query, get='xy'):
        # Usage:
        # #1 sources_ccd72_#1 = DpdExtObsCol().SourceTableDetID(OC_ImageObjects, 'ESO_CCD_#72')
        # #2 sources_ccd72_#2 = DpdExtObsCol().SourceTableDetID(OC_ImageObjects, query_exp1ccd72_stage2)
        try:
            [xmax, xmin], [ymax, ymin], detId = self.XYlimDetId(query)
        except AttributeError as err:
            'Handling AttributeError:', err, '\n --> Input query =', query, ' was no query.'
            detId = query
        except TypeError as err:
            'Handling TypeError:', err, '\n --> Input query =', query, ' was no query.'
            detId = query
        mask_detId = (LDAC_OBJECTS['DETECTOR'] == detId)
        if get == 'xy':
            XWIN_IMAGE_DET = LDAC_OBJECTS['XWIN_IMAGE'][mask_detId]
            YWIN_IMAGE_DET = LDAC_OBJECTS['YWIN_IMAGE'][mask_detId]
            table = [XWIN_IMAGE_DET, YWIN_IMAGE_DET]
        elif get == 'all':
            table = LDAC_OBJECTS[:][mask_detId]
        else:
            table = LDAC_OBJECTS[get][mask_detId]
        return table


class DpdExtObsCol(Select_FITS_from_DpdExtObsCol):
    def __init__(self):
        Select_FITS_from_DpdExtObsCol.__init__(self)


# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------


def dss_io():
    from euclid.toolbox.ingest_client_sc456 import DataIO
    from euclid.config.Environment import Env
    dss_server = Env['dss_server']
    dss_username = Env['dss_username']
    dss_password = Env['dss_password']
    ds = DataIO(dss_server, username=dss_username, password=dss_password)
    return ds


def retrieve_from_dss(filename):
    dss = dss_io()
    res = dss.get(path=filename, savepath='./%s' % filename)
    return res


