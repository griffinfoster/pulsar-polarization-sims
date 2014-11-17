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
    o.set_usage('%prog [options] [TIM files]')
    o.set_description(__doc__)
    o.add_option('--mjd',dest='mjdFile',default=None,
        help='Expected MJD file')
    opts, args = o.parse_args(sys.argv[1:])

    timTxtFiles=[]
    timPklFiles=[]
    for fn in args:
        if fn.endswith('txt'): timTxtFiles.append(fn)
        elif fn.endswith('pkl'): timPklFiles.append(fn)
    timTxtDict=simsio.timFromTxt(timTxtFiles)
    timPklDict=simsio.timFromPkl(timPklFiles)
    print timPklDict.keys()
    print timTxtDict.keys()

    if not(opts.mjdFile is None): mjds=simsio.mjdsFromTxt(opts.mjdFile)
