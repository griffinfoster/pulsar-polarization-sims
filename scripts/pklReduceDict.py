#!/usr/bin/env python
"""Combine MJD,TIM,RMS files into a single dictionary
"""

import os,sys
import numpy as np

import simsio
import cPickle as pkl

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] --mjd=MJDFILE [TIM files] [RMS files]')
    o.set_description(__doc__)
    o.add_option('-m','--mode',dest='mode',default='pkl',
        help='File mode: pkl or txt, default: pkl')
    o.add_option('--mjd',dest='mjdFile',default=None,
        help='MJD file')
    o.add_option('-o','--ofn',dest='ofn',default='combined.pkl',
        help='Output pickle file name, default:combined.pkl')
    opts, args = o.parse_args(sys.argv[1:])

    rmsFiles=[]
    timFiles=[]
    for fn in args:
        if fn.split('/')[-1].startswith('tim'): timFiles.append(fn)
        if fn.split('/')[-1].startswith('rms'): rmsFiles.append(fn)
    if opts.mode=='txt': timDict=simsio.timFromTxt(timFiles)
    else: timDict=simsio.timFromPkl(timFiles)
    rmsDict=simsio.rmsFromTxt(rmsFiles)

    expMJDs=simsio.mjdsFromTxt(opts.mjdFile)

    #construct a combined results dictionary and compute the useful metrics
    reduceDict={}
    for key in timDict.iterkeys():
        resDict={}
        #compute average sigma toa
        if timDict[key].size==0: #no solution for the simulation parameters
            print 'Warning simulation has no solution:',key
        else:
            if rmsDict[key] < 0: #no RMS values
                print 'Bad RMS values:', key
            else:
                resDict['rms']=rmsDict[key]
                validIdx=np.argwhere(timDict[key][:,0]>0) #pick out valid MJDs
                resDict['obsMJD']=timDict[key][validIdx,0]
                resDict['avgSigma']=np.mean(timDict[key][validIdx,1])
                resDict['sigmas']=timDict[key][validIdx,1]
                resDict['expMJD']=expMJDs[validIdx]
                resDict['nobs']=resDict['obsMJD'].shape[0]
                #compute chi^2
                chi2=(1./(resDict['obsMJD'].shape[0]-1))*np.sum((((24.*60.*60.*1e6)*(resDict['obsMJD']-resDict['expMJD']))**2.)/(resDict['sigmas']**2.))
                resDict['chi2']=chi2
                reduceDict[key]=resDict

    #save dictionary to pickle file
    pkl.dump(reduceDict,open(opts.ofn,'wb'))

