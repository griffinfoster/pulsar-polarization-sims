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
    o.add_option('-l','--leak',dest='leak', action='store_true',
        help='Plot in terms of polarization leakage instead of IXR')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save figure in a format based on name extension')
    o.add_option('-S','--show',dest='show', action='store_true',
        help='Show the plot')
    o.add_option('--snr',dest='snr',default=100,type='int',
        help='SNR value to use (rounds to nearest int value), default: 100')
    o.add_option('--info',dest='info',action='store_true',
        help='Print parameter information in the dictionary and exit')
    o.add_option('--lines',dest='lines',action='store_true',
        help='Instead of fill plots, plot the individual lines')
    o.add_option('--ixrlines',dest='ixrlines',action='store_true',
        help='Instead of fill plots, plot the individual IXR lines')
    o.add_option('--dJmax',dest='dJmax',default=0.15,type='float',
        help='Maximum calibration error to use, default: 0.15')
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
    
    #Total Intensity
    ixrdbs0=[]
    polLeakdbs0=[]
    deltaJs0=[]
    rmsVals0=[]
    #Invariant Interval
    ixrdbs1=[]
    polLeakdbs1=[]
    deltaJs1=[]
    rmsVals1=[]
    #MTM cal
    ixrdbs2cal=[]
    polLeakdbs2cal=[]
    deltaJs2cal=[]
    rmsVals2cal=[]
    #MTM uncal
    ixrdbs2uncal=[]
    polLeakdbs2uncal=[]
    deltaJs2uncal=[]
    rmsVals2uncal=[]
    for key,val in reduceDict.iteritems():
        if int(key[1])==opts.snr and key[2] < opts.dJmax: #SNR mode and dJ max selection
            deltaJ=key[2]*100.
            polLeakdb=10.*np.log10((key[3]**2))
            ixrdb=10.*np.log10(1./(key[3]**2))
            if key[0]==0 and key[4].startswith('cal'):
                ixrdbs0.append(ixrdb)
                deltaJs0.append(deltaJ)
                polLeakdbs0.append(polLeakdb)
                rmsVals0.append(val['rms'])
                #rmsVals0.append(val['chi2'])
                #rmsVals0.append(val['avgSigma'])
            elif key[0]==1 and key[4].startswith('cal'):
                ixrdbs1.append(ixrdb)
                deltaJs1.append(deltaJ)
                polLeakdbs1.append(polLeakdb)
                rmsVals1.append(val['rms'])
                #rmsVals1.append(val['chi2'])
                #rmsVals1.append(val['avgSigma'])
            elif key[0]==2:
                if key[4].startswith('cal'):
                    ixrdbs2cal.append(ixrdb)
                    deltaJs2cal.append(deltaJ)
                    polLeakdbs2cal.append(polLeakdb)
                    rmsVals2cal.append(val['rms'])
                else:
                    ixrdbs2uncal.append(ixrdb)
                    deltaJs2uncal.append(deltaJ)
                    polLeakdbs2uncal.append(polLeakdb)
                    rmsVals2uncal.append(val['rms'])

    def fillPlotter(ixrdbs,polLeakdbs,deltaJs,rmsVals):
        ixrdbs=np.array(ixrdbs)
        polLeakdbs=np.array(polLeakdbs)
        deltaJs=np.array(deltaJs)
        rmsVals=np.array(rmsVals)
        polLeakVals=np.unique(polLeakdbs)
        sortIdx=np.argsort(polLeakVals)
        avgRMS=np.zeros_like(polLeakVals)
        minRMS=np.zeros_like(polLeakVals)
        maxRMS=np.zeros_like(polLeakVals)
        for pid,pval in enumerate(polLeakVals):
            rmsValsPolLeak=[]
            for pid0,pval0 in enumerate(polLeakdbs):
                if pval==pval0: rmsValsPolLeak.append(rmsVals[pid0])
            avgRMS[pid]=np.average(np.array(rmsValsPolLeak))
            minRMS[pid]=np.min(np.array(rmsValsPolLeak))
            maxRMS[pid]=np.max(np.array(rmsValsPolLeak))
        polLeakVals=polLeakVals[sortIdx]
        maxRMS=maxRMS[sortIdx]
        minRMS=minRMS[sortIdx]
        
        #interpolate to smooth out curves
        fmax=interpolate.interp1d(polLeakVals,maxRMS,kind='linear')
        fmin=interpolate.interp1d(polLeakVals,minRMS,kind='linear')
        interpPolLeak=np.linspace(polLeakVals[0],polLeakVals[-1],400)
        window_len=31
        return interpPolLeak,smooth(fmax(interpPolLeak),window_len=window_len)[:-30],smooth(fmin(interpPolLeak),window_len=window_len)[:-30]

    def linePlotter(ixrdbs,polLeakdbs,deltaJs,rmsVals):
        ixrdbs=np.array(ixrdbs)
        polLeakdbs=np.array(polLeakdbs)
        deltaJs=np.array(deltaJs)
        rmsVals=np.array(rmsVals)
        deltaJVals=np.unique(deltaJs)
        rmsLines=[]
        polLeakLines=[]
        for dJ in deltaJVals:
            idx=np.argwhere(deltaJs==dJ)
            subPolLeak=polLeakdbs[idx][:,0]
            subRmsVals=rmsVals[idx][:,0]
            sortIdx=np.argsort(subPolLeak)
            rmsLines.append(subRmsVals[sortIdx])
            polLeakLines.append(subPolLeak[sortIdx])
        return polLeakLines,rmsLines,deltaJVals

    def ixrLinePlotter(ixrdbs,polLeakdbs,deltaJs,rmsVals):
        ixrdbs=np.array(ixrdbs)
        polLeakdbs=np.array(polLeakdbs)
        deltaJs=np.array(deltaJs)
        rmsVals=np.array(rmsVals)
        polLeakVals=np.unique(polLeakdbs)
        rmsLines=[]
        dJLines=[]
        for pid,pval in enumerate(polLeakVals):
            idx=np.argwhere(polLeakdbs==pval)
            subdJs=deltaJs[idx][:,0]
            subRmsVals=rmsVals[idx][:,0]
            sortIdx=np.argsort(subdJs)
            rmsLines.append(subRmsVals[sortIdx])
            dJLines.append(subdJs[sortIdx])
        return dJLines,rmsLines,polLeakVals

    fig=p.figure()
    fig.set_size_inches(9.,6.)
    ax=fig.add_subplot(1,1,1)

    if opts.lines:
        polLeakLines,rmsLines,deltaJVals=linePlotter(ixrdbs0,polLeakdbs0,deltaJs0,rmsVals0)
        #color setup
        cVals=np.log10(np.array(rmsLines)[:,0])
        cVals-=np.min(cVals)
        cVals/=np.max(cVals)
        textStep=0.
        for pid,lid,did,cid in zip(polLeakLines,rmsLines,deltaJVals,cVals):
            if did==23.0: continue
            p.plot(pid,lid,color=(cid,0.,1-cid))
            if textStep%2==0: p.text(pid[int(pid.shape[0]*.4)]-textStep,lid[int(lid.shape[0]*.4)],'%0.f%%'%did,verticalalignment='center')
            textStep+=1.
        p.xlim(-30,0)

    elif opts.ixrlines:
        dJLines,rmsLines,polLeakVals=ixrLinePlotter(ixrdbs0,polLeakdbs0,deltaJs0,rmsVals0)
        #color setup
        cVals=np.log10(np.array(rmsLines)[:,-1])
        cVals-=np.min(cVals)
        cVals/=np.max(cVals)
        for djid,lid,pid,cid in zip(dJLines,rmsLines,polLeakVals,cVals):
            if pid > -1. and pid < -0.1: #hack
                templid=lid
                templid[-2]=(templid[-3]+templid[-1])/2.
                p.plot(djid,templid,color=(cid,0.,1-cid))
            else: p.plot(djid,lid,color=(cid,0.,1-cid))
            if pid > -0.1: lblStr='%0.f dB'%(-1.*pid)
            elif (pid >-30. and pid < -11) or pid < -31: lblStr=''
            else: lblStr='%0.f dB'%(pid)
            p.text(djid[int(djid.shape[0]*.25)],lid[int(lid.shape[0]*.25)],lblStr,verticalalignment='center')
        p.gca().invert_xaxis()

    else:
        polLeakVals,maxRMS0,minRMS0=fillPlotter(ixrdbs0,polLeakdbs0,deltaJs0,rmsVals0)
        p.fill_between(polLeakVals,maxRMS0,y2=minRMS0,edgecolor='none',facecolor=(1.,0.5,0.5,1.))
    
        polLeakVals,maxRMS1,minRMS1=fillPlotter(ixrdbs1,polLeakdbs1,deltaJs1,rmsVals1)
        p.fill_between(polLeakVals,maxRMS1,y2=minRMS1,edgecolor='none',facecolor=(0.5,0.5,0.5,.9))
    
        polLeakVals,maxRMS2,minRMS2=fillPlotter(ixrdbs2cal,polLeakdbs2cal,deltaJs2cal,rmsVals2cal)
        p.fill_between(polLeakVals,maxRMS2,y2=minRMS2,edgecolor='none',facecolor=(1.,0.73,0.25,.9))
    
        polLeakVals,maxRMS3,minRMS3=fillPlotter(ixrdbs2uncal,polLeakdbs2uncal,deltaJs2uncal,rmsVals2uncal)
        p.fill_between(polLeakVals,maxRMS3,y2=minRMS3,edgecolor='none',facecolor=(0.5,0.63,1.,.9))

        p.xlim(-30,0)

    if opts.leak: p.xlabel('polarization leakage (dB)',fontsize=fs)
    else: p.xlabel('IXR (dB)',fontsize=fs)
    p.ylabel('rms ($\mu$s)',fontsize=fs)

    if opts.ixrlines: p.xlabel('Calibration Error (%)',fontsize=fs)

    ax.set_yscale('log')

    #maxPolLeak=-3.
    #p.xlim(np.min(np.array(polLeakVals)),np.max(np.array(polLeakVals))+maxPolLeak)

    #minIdx=(np.abs(np.array(polLeakVals)-maxPolLeak)).argmin()
    #maxIdx=np.array(polLeakVals).argmin()
    #p.ylim(min(minRMS1[maxIdx],minRMS2[maxIdx],minRMS3[maxIdx])-.1,max(maxRMS1[minIdx],maxRMS2[minIdx],maxRMS3[minIdx]))

    #p.xlim()
    #p.ylim()

    p.gca().invert_xaxis()

    if opts.show: p.show()
    if not(opts.savefig is None):
        p.savefig(opts.savefig, bbox_inches='tight')

