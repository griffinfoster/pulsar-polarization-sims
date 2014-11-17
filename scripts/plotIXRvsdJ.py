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
    o.add_option('-g','--grid',dest='grid', default='0,30,0,40',
        help='Grid limits for the plot: deltaJ min (%), deltaJ max (%), IXR min (dB), IXR max (dB) default: 0,30,0,40 ')
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
    o.add_option('--snr',dest='snr',default=100,type='int',
        help='SNR value to use (rounds to nearest int value), default: 100')
    o.add_option('--scatter',dest='scatter',action='store_true',
        help='Scatter plot of raw data as opposed to an interpolated contour plot')
    o.add_option('--min',dest='min',default=None,type='float',
        help='Clip values to this minimum')
    o.add_option('--max',dest='max',default=None,type='float',
        help='Clip values to this maximum')
    opts, args = o.parse_args(sys.argv[1:])

    print 'Loading PKL file'
    reduceDict=pkl.load(open(args[0]))

    limits=map(int,opts.grid.split(','))
    dJmin=limits[0]
    dJmax=limits[1]
    IXRmin=limits[2]
    IXRmax=limits[3]

    ixrdbs=[]
    polLeakdbs=[]
    deltaJs=[]
    simVals=[]
    nobs=[]
    print 'Selecting subset'
    for key,val in reduceDict.iteritems():
        if key[0]==opts.rmsMode and int(key[1])==opts.snr and key[4].startswith(opts.calMode): #RMS, CAL and SNR mode selection
            deltaJ=key[2]*100.
            polLeakdb=10.*np.log10((key[3]**2))
            ixrdb=10.*np.log10(1./(key[3]**2))
            if deltaJ >= dJmin and deltaJ <= dJmax and ixrdb >= IXRmin and ixrdb <= IXRmax: #min/max delta J and IXR dB selection
                #data product slection
                #val keys: ['rms', 'chi2', 'avgSigma', 'obsMJD', 'nobs', 'expMJD', 'sigmas']
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

    if opts.leak:
        ixrdbs=np.array(polLeakdbs)
        temp=IXRmax
        IXRmax=-1.*IXRmin
        IXRmin=-1.*temp
    else: ixrdbs=np.array(ixrdbs)
    deltaJs=np.array(deltaJs)
    simVals=np.array(simVals)

    if opts.scatter:
        #scatter plot
        print 'Generating scatter plot'
        scatSize=20.
        p.scatter(deltaJs,ixrdbs,simVals*scatSize,c=simVals,edgecolors='none')
        p.colorbar()
    else:
        #contour plot
        print 'Generating interpolated contour plot'
        import scipy.interpolate
        grid_x,grid_y=np.mgrid[dJmin:dJmax:20j*int(abs(dJmax-dJmin)), IXRmin:IXRmax:20j*int(abs(IXRmax-IXRmin))]
        grid_z=scipy.interpolate.griddata(( deltaJs, ixrdbs ), simVals, (grid_x,grid_y), method='linear')
        #set contour levels and labels
        #get the number of orders of magnitude in the dynamic range
        maxOrder=int(np.log10(np.nanmax(simVals)))
        minOrder=int(np.log10(np.nanmin(simVals)))-1
        lvls0=[]
        #make 4 contour levels for each magnitude:
        fmt={}
        for mag in range(minOrder,maxOrder+1):
            maglvls=4
            magstep=1./maglvls
            for i in range(maglvls):
                l0=10.**(mag+(i*magstep))
                lvls0.append(l0)
                if mag<0: fmt[l0]='%.3f'%(int(l0/(10.**mag))*(10.**mag))
                else: fmt[l0]='%i'%(int(l0/(10.**mag))*(10.**mag))
        CS=p.contour(grid_x,grid_y,grid_z,lvls0)
        #CS=p.contourf(grid_x,grid_y,grid_z,lvls0)
        p.clabel(CS,CS.levels[::4],inline=0,fontsize=25,colors='black',fmt=fmt) #label every 4th level

        #p.imshow(np.rot90(grid_z),interpolation='nearest',extent=(dJmin,dJmax,IXRmin,IXRmax))
        #p.colorbar()

    p.xlabel('calibration error (%)',fontsize=fs)
    p.xlim(dJmin,dJmax)
    if opts.leak: p.ylabel('polarization leakage (dB)',fontsize=fs)
    else: p.ylabel('IXR (dB)',fontsize=fs)
    p.ylim(IXRmin,IXRmax)
    p.title(modeTitle[opts.rmsMode],fontsize=fs)

    if opts.show: p.show()
    if not(opts.savefig is None):
        fig=p.gcf()
        fig.set_size_inches(9.,6.)
        p.savefig(opts.savefig, bbox_inches='tight')

