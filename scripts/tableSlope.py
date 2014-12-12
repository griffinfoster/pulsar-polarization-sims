#!/usr/bin/env python
"""
"""

import os,sys
import numpy as np

import cPickle as pkl

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [pklReduceDict.py DICT]')
    o.set_description(__doc__)
    o.add_option('--snr',dest='snr',default=100,type='int',
        help='SNR value to use (rounds to nearest int value), default: 100')
    o.add_option('--info',dest='info',action='store_true',
        help='Print parameter information in the dictionary and exit')
    o.add_option('--dJ',dest='dJ',default=0.05,type='float',
        help='Calibration error to select out, default: 0.05')
    o.add_option('-c','--cal',dest='calMode',default='cal',
        help='cal mode to use: cal or uncal, default: cal')
    o.add_option('-m','--mode',dest='mode',default='rms',
        help='Data mode: rms, chi2, sigma ; default: rms')
    o.add_option('-r','--rms', dest='rmsMode', default=0, type='int',
        help='Set RMS mode, 0: total intesity, 1: invariant interval, 2: matrix template matching. default: 0')
    opts, args = o.parse_args(sys.argv[1:])

    print 'Loading PKL file'
    reduceDict=pkl.load(open(args[0]))

    if opts.info:
        snrs=[]
        deltaJs=[]
        ixrs=[]
        for key,val in reduceDict.iteritems():
            snrs.append(key[1])
            deltaJs.append(key[2]*100.)
            ixrs.append(10.*np.log10(1./(key[3]**2)))
        snrs=np.array(snrs)
        deltaJs=np.array(deltaJs)
        ixrs=np.array(ixrs)
        print 'SNR:', np.unique(snrs)
        print 'delta J (\%):',np.unique(deltaJs)
        print 'IXR (dB):', np.unique(ixrs)
        exit()
    
    ixrdbs=[]
    vals=[]
    for key,val in reduceDict.iteritems():
        #key: (mode,snr,dJ,IXR,cal/uncal)
        #val keys: ['rms', 'chi2', 'avgSigma', 'obsMJD', 'nobs', 'expMJD', 'sigmas']
        if key[0]==opts.rmsMode and int(key[1])==opts.snr and key[2]==opts.dJ and key[4].startswith(opts.calMode): #timing mode, snr, dJ, cal mode selection
            ixrdb=10.*np.log10(1./(key[3]**2))
            ixrdbs.append(ixrdb)
            if opts.mode.startswith('rms'): vals.append(val['rms'])
            elif opts.mode.startswith('chi'): vals.append(val['chi2'])
            elif opts.mode.startswith('sigma'): vals.append(val['avgSigma'])

    ixrdbs=np.array(ixrdbs)
    vals=np.array(vals)

    idx=np.argsort(ixrdbs)
    print 'IXR',ixrdbs[idx]
    print 'RMS',vals[idx]
    print 'precent',100.*np.diff(vals[idx])/vals[idx][:-1]

