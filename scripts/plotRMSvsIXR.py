#!/usr/bin/env python
"""Generate an IXR/Polarization Leakage vs. delta J plot for RMS/average sigma_ToA/chi^2 from the dictionary produced with pklReduceDict.py
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
    o.add_option('-g','--grid',dest='grid', default='0,3000,0,40',
        help='Grid limits for the plot: SNR min , SNR max , IXR min (dB), IXR max (dB) default: 0,3000,0,40 ')
    o.add_option('-m','--mode',dest='mode',default='rms',
        help='Data mode: rms, chi2, sigma ; default: rms')
    o.add_option('-l','--leak',dest='leak', action='store_true',
        help='Plot in terms of polarization leakage instead of IXR')
    o.add_option('-r','--rms', dest='rmsMode', default=0, type='int',
        help='Set RMS mode, 0: total intesity, 1: invariant interval, 2: matrix template matching. default: 0')
    o.add_option('--dJ',dest='dJ',default=0.05,type='float',
        help='Calibration error to select out, default: 0.05')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save figure in a format based on name extension')
    o.add_option('-S','--show',dest='show', action='store_true',
        help='Show the plot')
    o.add_option('--scatter',dest='scatter',action='store_true',
        help='Scatter plot of raw data as opposed to an interpolated contour plot')
    o.add_option('--min',dest='min',default=None,type='float',
        help='Clip values to this minimum')
    o.add_option('--max',dest='max',default=None,type='float',
        help='Clip values to this maximum')
    o.add_option('--info',dest='info',action='store_true',
        help='Print parameter information in the dictionary and exit')
    o.add_option('--outliers',dest='filtout',action='store_true',
        help='Filter outliers before plotting')
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

    limits=map(int,opts.grid.split(','))
    SNRmin=limits[0]
    SNRmax=limits[1]
    IXRmin=limits[2]
    IXRmax=limits[3]

    ixrdbs=[]
    polLeakdbs=[]
    deltaJs=[]
    simVals=[]
    nobs=[]
    snrs=[]
    print 'Selecting subset'
    for key,val in reduceDict.iteritems():
        if key[0]==opts.rmsMode and key[4].startswith(opts.calMode) and key[2]==opts.dJ: #RMS, CAL, dJ selection
            deltaJ=key[2]*100.
            polLeakdb=10.*np.log10((key[3]**2))
            ixrdb=10.*np.log10(1./(key[3]**2))
            snr=key[1]
            #if deltaJ >= dJmin and deltaJ <= dJmax and ixrdb >= IXRmin and ixrdb <= IXRmax: #min/max delta J and IXR dB selection
            #data product slection
            #val keys: ['rms', 'chi2', 'avgSigma', 'obsMJD', 'nobs', 'expMJD', 'sigmas']
            if np.isnan(val['rms']) or np.isnan(val['avgSigma']) or val['rms']<0.:
                print 'Warning: some selected values are NaN',key
                continue #skip any simulations which returned NaN
            ixrdbs.append(ixrdb)
            deltaJs.append(deltaJ)
            polLeakdbs.append(polLeakdb)
            snrs.append(snr)
            if opts.mode.startswith('rms'): simVals.append(val['rms'])
            elif opts.mode.startswith('chi'): simVals.append(val['chi2'])
            elif opts.mode.startswith('sigma'): simVals.append(val['avgSigma'])
            nobs.append(val['nobs'])

    if opts.leak:
        ixrdbs=np.array(polLeakdbs)
        #temp=IXRmax
        #IXRmax=-1.*IXRmin
        #IXRmin=-1.*temp
    else: ixrdbs=np.array(ixrdbs)
    simVals=np.array(simVals)
    snrs=np.array(snrs)

    #filter outliers
    if opts.filtout:
        medianVal=np.median(simVals)
        filterFactor=5. #filter out values above filterFactor*medianVal
        idx=np.argwhere(simVals>filterFactor*medianVal)
        simVals[idx]=np.nan

    if opts.scatter:
        #scatter plot
        print 'Generating scatter plot'
        scatSize=20.
        p.scatter(snrs,ixrdbs,simVals*scatSize,c=simVals,edgecolors='none')
        p.colorbar()
        #print np.median(simVals.flatten())
        #print np.mean(simVals.flatten())
        #print np.std(simVals.flatten())
        #print np.max(simVals.flatten())
        #print np.min(simVals.flatten())
        #p.plot(simVals.flatten(),'b.')
    else:
        #contour plot
        print 'Generating interpolated contour plot'
        import scipy.interpolate
        grid_x,grid_y=np.mgrid[SNRmin:SNRmax:20j*int(abs(SNRmax-SNRmin)), IXRmin:IXRmax:20j*int(abs(IXRmax-IXRmin))]
        grid_z=scipy.interpolate.griddata(( deltaJs, ixrdbs ), simVals, (grid_x,grid_y), method='linear')
        #set contour levels and labels
        #get the number of orders of magnitude in the dynamic range
        maxOrder=int(np.log10(np.nanmax(simVals)))
        minOrder=int(np.log10(np.nanmin(simVals)))-1
        lvls0=[]
        lvlLbls=[]
        fmt={}
        for mag in range(minOrder,maxOrder+1):
            maglvls=6 #make 6 contour levels for each magnitude
            magstep=1./maglvls
            for i in range(maglvls):
                l0=10.**(mag+(i*magstep))
                lvls0.append(l0)
                if mag<0: fmt[l0]='%.3f'%(int(l0/(10.**mag))*(10.**mag))
                else: fmt[l0]='%i'%(int(l0/(10.**mag))*(10.**mag))
        CS=p.contour(grid_x,grid_y,grid_z,lvls0)
        if np.abs(maxOrder-minOrder)>2:
            p.clabel(CS,CS.levels[maglvls*2::4],inline=0,fontsize=25,colors='black',fmt=fmt) #label every 4th level
            p.clabel(CS,CS.levels[:2*maglvls],inline=0,fontsize=25,colors='black',fmt=fmt) #label every level of the lowest 2 orders
        else: p.clabel(CS,CS.levels,inline=0,fontsize=25,colors='black',fmt=fmt) #label every level if there are feewer than 2 orders

        #p.imshow(np.rot90(grid_z),interpolation='nearest',extent=(dJmin,dJmax,IXRmin,IXRmax))
        #p.colorbar()

    p.xlabel('calibration error (%)',fontsize=fs)
    #dJbuffer=.1*np.abs((dJmax-dJmin))
    #p.xlim(dJmin-dJbuffer,dJmax+dJbuffer)
    if opts.leak: p.ylabel('polarization leakage (dB)',fontsize=fs)
    else: p.ylabel('IXR (dB)',fontsize=fs)
    #ixrBuffer=.1*np.abs((IXRmax-IXRmin))
    #p.ylim(IXRmin-ixrBuffer,IXRmax+ixrBuffer)
    p.title('%s [%s,%s]'%(modeTitle[opts.rmsMode],opts.mode,opts.calMode),fontsize=fs)

    if opts.show: p.show()
    if not(opts.savefig is None):
        fig=p.gcf()
        fig.set_size_inches(9.,6.)
        p.savefig(opts.savefig, bbox_inches='tight')

