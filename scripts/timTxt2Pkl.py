#!/usr/bin/env python
"""Simulate Pulsar Timing ToA's using PSARCHIVE abd TEMPO2
"""

import numpy as np
import os, sys
import cPickle as pickle

def readTimFile(fn):
    """Read a TIM output file from PAT, and return a numpy array of the contents"""
    #Some dodgy stuff is done here to maintain the numerical precision, but it does seem to work
    fh=open(fn,'r')
    timData=fh.read()
    fh.close()
    lines=timData.split('\n')
    lines=lines[1:-1]
    arr=[]
    for l in lines:
        if l.startswith('FORMAT'): continue
        splitLine=l.split()
        mjd=splitLine[2].split('.')
        #print mjd[0],mjd[1],'%.14f'%(float(mjd[0])+float(mjd[1][:12])/(1e12))
        ndecimals=len(mjd[1])
        cmjd=float(mjd[0])+float(mjd[1])/(10.**ndecimals)
        arr.append([cmjd,float(splitLine[3])])
    arr=np.array(arr)
    return arr

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] tim*.txt')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    for fn in args:
        arr=readTimFile(fn)
        ofn=fn.split('/')[-1].split('.txt')[0]+'.pkl'
        ofh=open(ofn,'wb')
        pickle.dump(arr,ofh)
        ofh.close()

