#!/usr/bin/env python

import os,sys
import numpy as np

import simsio

#import matplotlib
#matplotlib.use('Agg')
#import pylab as p
#matplotlib.rc('xtick',labelsize=25)
#matplotlib.rc('ytick',labelsize=25)

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [RMS txt files]')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    rmsFiles=args
    if rmsFiles[0]=='RMS.dat':
        rmsDict=simsio.loadTotalRMSDat(rmsFiles[0])
    elif rmsFiles[0]=='RMS.pkl':
        rmsDict=simsio.loadTotalRMSPkl(rmsFiles[0])
    else:
        rmsDict=simsio.rmsFromTxt(rmsFiles)
    print rmsDict

