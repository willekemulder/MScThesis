
from __future__ import print_function

import numpy as np


import sys
import os

import pickle



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
# 
# ----------------------------------------------------------------------


def get_KiDS_Gaia_Test_Patch():
    "Obtain data for the KiDS-Gaia-Test-Patch area."
    "Returns dictionary with tilenames and ObservingBlock objects."
    tilename_central_pointing=['KIDS_183.5_-2.5']
    tilenames_full_patch = [
    'KIDS_180.0_1.5', 'KIDS_180.0_0.5', 'KIDS_180.0_-0.5', 'KIDS_180.0_-1.5', 'KIDS_180.5_-2.5',
    'KIDS_181.0_1.5', 'KIDS_181.0_0.5', 'KIDS_181.0_-0.5', 'KIDS_181.0_-1.5', 'KIDS_181.5_-2.5',
    'KIDS_182.0_1.5', 'KIDS_182.0_0.5', 'KIDS_182.0_-0.5', 'KIDS_182.0_-1.5', 'KIDS_182.5_-2.5',
    'KIDS_183.0_1.5', 'KIDS_183.0_0.5', 'KIDS_183.0_-0.5', 'KIDS_183.0_-1.5', 'KIDS_183.5_-2.5',
    'KIDS_184.0_1.5', 'KIDS_184.0_0.5', 'KIDS_184.0_-0.5', 'KIDS_184.0_-1.5', 'KIDS_184.5_-2.5' ]
    tilenames=tilename_central_pointing
    obnames = []
    for filtername in ['r']:
        for tilename in tilenames:
            obnames.append('%s_%s' % (tilename, filtername))
    KiDSGaiaTestPatch={'tilenames':tilenames,'obnames':obnames}
    return KiDSGaiaTestPatch





# ----------------------------------------------------------------------
# 
# ----------------------------------------------------------------------


class cdfs:
    '''
    #Custodian: Gijs Verdoes Kleijn: g.a.verdoes.kleijn@.rug.nl / tel:+31-50-3638326.

    >> os.chdir("../../pythonscripts/")

    >> import select_from_AWEdb
    >> cdfs=select_from_AWEdb.cdfs()
    >> cdfs.get_coadds()

    >> os.chdir("../../pythonscripts/")
    >> os.chdir(cwd)
    '''

    def get_coadds(self):
        '''
        The pickle was created using https://github.com/willekemulder/masterproject/blob/master/notebooks/awe/GijsMasterThesisWillekeAstroWISE.ipynb
        '''
        f=open('coadds_cdfs.pickle','rb')
        objectids=pickle.load(f)
        f.close()
        self.coadds=[CoaddedRegriddedFrame(object_id=objectid) for objectid in objectids]
