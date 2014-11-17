"""
Pulsar Polarization Simulation Results I/O interfaces
"""

import cPickle as pkl
import numpy as np

id2mode=['Total Intensity','Invariant Interval','Matrix Template Matching']
mode2id={'intensity':0, 'invariant':1, 'matrixtemplate':2}

def rmsFromTxt(rmsFiles):
    """Return a dictionary of RMS values (in nanoseconds) from a list of files
    """
    rmsDict={}
    for fn in rmsFiles:
        #parse file name
        #FORMAT: rms<MODE>_sn<SNR>_dj<DJ>_pp<POLPUR>.<CAL/UNCAL>.txt
        calMode=fn.split('.')[-2]
        lfn=fn.split('/')[-1] #chop off any path information
        mode=int(lfn.split('_')[0][-1])
        snr=float(lfn.split('_')[1].split('sn')[-1])
        dj=float(lfn.split('_')[2].split('dj')[-1])
        polPur=lfn.split('_')[3].split('pp')[-1]
        polPur=float(polPur.split('.%s.txt'%calMode)[0])

        fh=open(fn,'r')
        rmsVal=fh.readline().strip()
        if rmsVal=='': rmsDict[(mode,snr,dj,polPur,calMode)]=-1. #No RMS was computed, seg fault during timing
        else: rmsDict[(mode,snr,dj,polPur,calMode)]=float(rmsVal)
        fh.close()
    return rmsDict

def loadTotalRMSDat(rmsFile):
    """Return a dictionary from the RMS.dat total RMS text file
    """
    rmsDict={}
    fh=open(rmsFile)
    fmt=fh.readline().strip()
    #print fmt
    #assumes format: calMode S/N deltaJ PolPur rms0 rms1 rms2
    for line in fh.readlines():
        rmsList=line.strip().split(' ')
        rmsDict[(float(rmsList[1]),float(rmsList[2]),float(rmsList[3]),rmsList[0])]=map(float,rmsList[4:])
    return rmsDict

def loadTotalRMSPkl(rmsFile):
    """Return a dictionary from the RMS.pkl total RMS file
    """
    return pkl.load(open(rmsFile))
    #('cal', 1000.0, 0.29999999999999999, 0.031622776601683791): (58.729, 9.131, 9.175)

def timFromPkl(timFiles):
    """Return a dictionary of [time of arrival (in MJD), sigma ToA (in microseconds)]xNumber of Observations from a list of pickle files
    """
    timDict={}
    for fn in timFiles:
        #parse file name
        #FORMAT: tim<MODE>_sn<SNR>_dj<DJ>_pp<POLPUR>.<CAL/UNCAL>.pkl
        calMode=fn.split('.')[-2]
        lfn=fn.split('/')[-1] #chop off any path information
        mode=int(lfn.split('_')[0][-1])
        snr=float(lfn.split('_')[1].split('sn')[-1])
        dj=float(lfn.split('_')[2].split('dj')[-1])
        polPur=lfn.split('_')[3].split('pp')[-1]
        polPur=float(polPur.split('.%s.pkl'%calMode)[0])

        fh=open(fn)
        timArr=pkl.load(fh)
        timDict[(mode,snr,dj,polPur,calMode)]=timArr
        fh.close()
    return timDict

def timFromTxt(timFiles):
    """Return a dictionary of [time of arrival (in MJD), sigma ToA (in microseconds)]xNumber of Observations from a list of text files
    """
    #line formatexample: simPsr0000.uncal.fits 1272.000 50000.277181767857425   0.992  7
    #                       FITS_NAME           FREQ        OBSERVED MJD        SIGMA_TOA   PARKES_ID
    timDict={}
    for fn in timFiles:
        #parse file name
        #FORMAT: tim<MODE>_sn<SNR>_dj<DJ>_pp<POLPUR>.<CAL/UNCAL>.txt
        calMode=fn.split('.')[-2]
        lfn=fn.split('/')[-1] #chop off any path information
        mode=int(lfn.split('_')[0][-1])
        snr=float(lfn.split('_')[1].split('sn')[-1])
        dj=float(lfn.split('_')[2].split('dj')[-1])
        polPur=lfn.split('_')[3].split('pp')[-1]
        polPur=float(polPur.split('.%s.txt'%calMode)[0])

        timArr=np.loadtxt(fn,usecols=[2,3],skiprows=1)
        if timArr.size==0: timArr=np.array([np.nan,np.nan])
        timDict[(mode,snr,dj,polPur,calMode)]=timArr
    return timDict

def mjdsFromTxt(mjdFile):
    """Return an array of expected MJD values (in days)
    """
    mjdArr=np.loadtxt(mjdFile)
    mjds=mjdArr[:,0]+((mjdArr[:,1]+mjdArr[:,2])/86400.)
    return mjds

