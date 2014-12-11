#!/usr/bin/env python
"""
"""

import os,sys
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import pylab as p

import cPickle as pkl

from scipy import interpolate

matplotlib.rc('xtick',labelsize=25)
matplotlib.rc('ytick',labelsize=25)

modeTitle=['Total Intensity','Invariant Interval','Matrix Template Matching']
fs=27 #fontsize

import numpy

def smooth(x,window_len=31,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [pklReduceDict.py DICT]')
    o.set_description(__doc__)
    o.add_option('-c','--cal',dest='calMode',default='cal',
        help='cal mode to use: cal or uncal, default: cal')
    o.add_option('-m','--mode',dest='mode',default='rms',
        help='Data mode: rms, chi2, sigma ; default: rms')
    o.add_option('-r','--rms', dest='rmsMode', default=0, type='int',
        help='Set RMS mode, 0: total intesity, 1: invariant interval, 2: matrix template matching. default: 0')
    o.add_option('-l','--leak',dest='leak', action='store_true',
        help='Plot in terms of polarization leakage instead of IXR')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save figure in a format based on name extension')
    o.add_option('-S','--show',dest='show', action='store_true',
        help='Show the plot')
    o.add_option('--dJ',dest='djval',default=5.0,type='float',
        help='Polarization calibration error value to use, default: 5.0')
    o.add_option('--info',dest='info',action='store_true',
        help='Print parameter information in the dictionary and exit')
    opts, args = o.parse_args(sys.argv[1:])

    print 'Loading PKL file'
    reduceDict=pkl.load(open(args[0]))

    normSNR=1000. #SNR to normalize other values against

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
    
    dJVal=opts.djval

    ixrdbs=[]
    polLeakdbs=[]
    deltaJs=[]
    simVals=[]
    snrs=[]
    nobs=[]
    print 'Selecting subset'
    for key,val in reduceDict.iteritems():
        if key[0]==opts.rmsMode and key[4].startswith(opts.calMode) and key[2]*100.==dJVal: #RMS, CAL, dJ mode selection
            #data product slection
            #val keys: ['rms', 'chi2', 'avgSigma', 'obsMJD', 'nobs', 'expMJD', 'sigmas']
            deltaJ=key[2]*100.
            polLeakdb=10.*np.log10((key[3]**2))
            ixrdb=10.*np.log10(1./(key[3]**2))
            snrs.append(key[1])
            if np.isnan(val['rms']) or np.isnan(val['avgSigma']) or val['rms']<0.:
                print 'Warning: some selected values are NaN',key
                continue #skip any simulations which returned NaN
            ixrdbs.append(ixrdb)
            deltaJs.append(deltaJ)
            polLeakdbs.append(polLeakdb)
            if opts.mode.startswith('rms'): simVals.append(val['rms'])
            elif opts.mode.startswith('chi'): simVals.append(val['chi2'])
            elif opts.mode.startswith('sigma'): simVals.append(val['avgSigma'])
            nobs.append(val['nobs'])

    simVals=np.array(simVals)
    polLeakdbs=np.array(polLeakdbs)
    snrs=np.array(snrs)/normSNR

    fig=p.figure()
    ax=fig.add_subplot(1,1,1)
    fig.set_size_inches(9.,5.)

    polLeakVals=np.unique(polLeakdbs)
    cNorm=np.min(polLeakVals)
    for pVal in polLeakVals:
        idx=np.argwhere(polLeakdbs==pVal)
        subSnrs=snrs[idx]
        subSimVals=simVals[idx]
        sIdx=np.argsort(subSnrs[:,0])
        rgb=(np.exp(pVal/10.),0.,1.-np.exp(pVal/10.))
        #slightly hardcoded plots and labels
        #if pVal > -0.0001: labelStr='%0.f dB'%(-1.*pVal) #-0 to 0
        if pVal > -0.1: continue
        elif (pVal < -16. and pVal > -30) or pVal < -31.: continue #skip these lines
        else: labelStr='%0.f dB'%(pVal)
        midPoint=int(subSnrs.shape[0]*.5)
        p.plot(subSnrs[sIdx,0],subSimVals[sIdx,0],color=rgb)
        p.text(subSnrs[sIdx,0][midPoint],0.8*subSimVals[sIdx,0][midPoint],labelStr,fontsize=14)

    #lines of constant time (5, 1, .5 us)
    p.hlines([5,1,.5],np.min(snrs),np.max(snrs),linestyles=['dashed','dashdot','dotted'])

    #slightly hardcoded plot limits
    #print np.min(snrs),np.max(snrs)
    #print np.min(simVals),np.max(simVals)
    p.xlim(np.min(snrs),np.max(snrs))
    p.ylim(np.min(simVals),100)
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'fraction of reference $\tau_{int}$', fontsize=23)
    ax.set_ylabel('rms $(\mu s)$', fontsize=23)

    if opts.savefig is None: p.show()
    else: p.savefig(opts.savefig, bbox_inches='tight')

