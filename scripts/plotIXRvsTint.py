#!/usr/bin/env python
"""
"""

import os,sys
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import pylab as p

import cPickle as pkl

matplotlib.rc('xtick',labelsize=25)
matplotlib.rc('ytick',labelsize=25)

modeTitle=['Total Intensity','Invariant Interval','Matrix Template Matching']
fs=27 #fontsize

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [pklReduceDict.py DICT]')
    o.set_description(__doc__)
    o.add_option('-c','--cal',dest='calMode',default='cal',
        help='cal mode to use: cal or uncal, default: cal')
    o.add_option('-g','--grid',dest='grid', default='0.01,1,0,40',
        help='Grid limits for the plot: T_int min, T_int max , IXR min (dB), IXR max (dB) default: 0.01,1,0,40 ')
    o.add_option('-m','--mode',dest='mode',default='rms',
        help='Data mode: rms, chi2, sigma ; default: rms')
    o.add_option('-l','--leak',dest='leak', action='store_true',
        help='Plot in terms of polarization leakage instead of IXR')
    o.add_option('-r','--rms', dest='rmsMode', default=0, type='int',
        help='Set RMS mode, 0: total intesity, 1: invariant interval, 2: matrix template matching. default: 0')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save figure in a format based on name extension')
    o.add_option('-S','--show',dest='show', action='store_true',
        help='Show the plot')
    o.add_option('--dJ',dest='dJ',default=0.05,type='float',
        help='Calibration error to select out, default: 0.05')
    o.add_option('--scatter',dest='scatter',action='store_true',
        help='Scatter plot of raw data as opposed to an interpolated contour plot')
    o.add_option('--min',dest='min',default=None,type='float',
        help='Clip values to this minimum')
    o.add_option('--max',dest='max',default=None,type='float',
        help='Clip values to this maximum')
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
        print 'tInt:', (np.unique(snrs)/1000.)**2.
        print 'delta J (\%):',np.unique(deltaJs)
        print 'IXR (dB):', np.unique(ixrs)
        exit()

    limits=map(float,opts.grid.split(','))
    tIntmin=limits[0]
    tIntmax=limits[1]
    IXRmin=limits[2]
    IXRmax=limits[3]

    ixrdbs=[]
    polLeakdbs=[]
    tInts=[]
    simVals=[]
    nobs=[]
    print 'Selecting subset'
    for key,val in reduceDict.iteritems():
        if key[0]==opts.rmsMode and key[2]==opts.dJ and key[4].startswith(opts.calMode): #RMS, CAL and delta J mode selection
            #tInt=np.sqrt(key[1]/1000.)
            tInt=(key[1]/1000.)**2.
            polLeakdb=10.*np.log10((key[3]**2))
            ixrdb=10.*np.log10(1./(key[3]**2))
            #if tInt >= tIntmin and tInt <= tIntmax and ixrdb >= IXRmin and ixrdb <= IXRmax: #min/max delta J and IXR dB selection
            #data product slection
            #val keys: ['rms', 'chi2', 'avgSigma', 'obsMJD', 'nobs', 'expMJD', 'sigmas']
            if np.isnan(val['rms']) or np.isnan(val['avgSigma']) or val['rms']<0.:
                print 'Warning: some selected values are NaN',key
                continue #skip any simulations which returned NaN
            ixrdbs.append(ixrdb)
            tInts.append(tInt)
            polLeakdbs.append(polLeakdb)
            if opts.mode.startswith('rms'): simVals.append(val['rms'])
            elif opts.mode.startswith('chi'): simVals.append(val['chi2'])
            elif opts.mode.startswith('sigma'): simVals.append(val['avgSigma'])
            nobs.append(val['nobs'])

    if opts.leak:
        ixrdbs=np.array(polLeakdbs)
        temp=IXRmax
        IXRmax=-1.*IXRmin
        IXRmin=-1.*temp
    else: ixrdbs=np.array(ixrdbs)
    tInts=np.array(tInts)
    simVals=np.array(simVals)

    tIntmin=np.min(tInts)
    tIntmax=np.max(tInts)

    fig=p.figure()
    ax=fig.add_subplot(1,1,1)
    fig.set_size_inches(9.,5.)

    if opts.scatter:
        #scatter plot
        print 'Generating scatter plot'
        scatSize=20.
        p.scatter(tInts,ixrdbs,simVals*scatSize,c=simVals,edgecolors='none')
        p.colorbar()
    else:
        #contour plot
        print 'Generating interpolated contour plot'
        import scipy.interpolate
        grid_x,grid_y=np.mgrid[tIntmin:tIntmax:1000j, IXRmin:IXRmax:20j*int(abs(IXRmax-IXRmin))]
        grid_z=scipy.interpolate.griddata(( tInts, ixrdbs ), simVals, (grid_x,grid_y), method='linear')
        #set contour levels and labels
        #get the number of orders of magnitude in the dynamic range
        maxOrder=int(np.log10(np.nanmax(simVals)))
        minOrder=int(np.log10(np.nanmin(simVals)))-1
        lvls0=[]
        lvlLbls=[]
        fmt={}
        for mag in range(minOrder,maxOrder-2):
            maglvls=6 #make 6 contour levels for each magnitude
            magstep=1./maglvls
            for i in range(maglvls):
                l0=10.**(mag+(i*magstep))
                lvls0.append(l0)
                #if mag<0: fmt[l0]='%.3f'%(int(l0/(10.**mag))*(10.**mag))
                if mag<0: fmt[l0]='%.1f'%(int(l0/(10.**mag))*(10.**mag))
                else: fmt[l0]='%i'%(int(l0/(10.**mag))*(10.**mag))
        #Select colors
        #print lvls0
        logLvls0=np.log10(lvls0)
        logLvls0-=np.min(logLvls0)
        logLvls0/=np.max(logLvls0)
        #print logLvls0
        rgbs=[]
        for ll in logLvls0: rgbs.append((ll,0.,1.-ll))
        grid_z=np.clip(grid_z,0.,5.)
        CS=p.contour(grid_x,grid_y,grid_z,lvls0,colors=rgbs)
        
        p.xlim(limits[0],limits[1])
        ax.set_xscale('log')

        if np.abs(maxOrder-minOrder)>2:
            p.clabel(CS,CS.levels[maglvls*2::4],inline=0,fontsize=25,colors='black',fmt=fmt) #label every 4th level
            p.clabel(CS,CS.levels[:2*maglvls],inline=0,fontsize=25,colors='black',fmt=fmt) #label every level of the lowest 2 orders
        else: p.clabel(CS,CS.levels,inline=0,fontsize=25,colors='black',fmt=fmt) #label every level if there are feewer than 2 orders

        #p.imshow(np.rot90(grid_z),interpolation='nearest',extent=(dJmin,dJmax,IXRmin,IXRmax))
        #p.colorbar()

    p.xlabel(r'fraction of reference $\tau_{int}$',fontsize=fs)
    p.ylim(IXRmin,IXRmax)

    if opts.leak: p.ylabel('polarization leakage (dB)',fontsize=fs)
    else: p.ylabel('IXR (dB)',fontsize=fs)
    #p.title('%s [%s,%s]'%(modeTitle[opts.rmsMode],opts.mode,opts.calMode),fontsize=fs)


    if opts.show: p.show()
    if not(opts.savefig is None):
        p.savefig(opts.savefig, bbox_inches='tight')

