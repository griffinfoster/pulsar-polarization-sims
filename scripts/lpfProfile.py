#!/usr/bin/env python
"""
Apply a low pass filter to a pulsar profile
"""

#broaden filter

import pyfits as pf
import numpy as n
import pylab as p
import os
import sys
import shutil
import time
from scipy import signal

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [FITS file]')
    o.set_description(__doc__)
    o.add_option('-c','--cutoff',dest='cutoff', type='float', default=0.1,
        help='Cutoff frequency, franction of Nyquist sampling, range 0 to 1. default: 0.1')
    o.add_option('-r', '--rot', dest='rot', action='store_true',
        help='Rotate the profile by 0.5 of the phase')
    o.add_option('-s','--save',dest='save',action='store_true',
        help='Save the LPF\'d profile to a new fits file')
    #o.add_option('-r', '--rot', dest='rot', action='store_true',
    #    help='Rotate the profile by 0.5 of the phase')
    opts, args = o.parse_args(sys.argv[1:])
    
    hdulist=pf.open(args[0])
    #print hdulist.info()

    primary=hdulist['PRIMARY'].header
    print primary['FITSTYPE']

    #see www.atnf.csiro.au/research/pulsar/psrfists/fitsdef.html section: Subintegration data
    d=hdulist[3].data
    #print d
    offsets=d[0][-3]
    sclFactor=d[0][-2]
    data=d[0][-1]
    #print sclFactor
    #print offsets
    #print data.shape
    if len(data.shape)==1:
        data.shape=(4,1,data.shape[-1]/4)
        #print data.shape

    dout=n.zeros_like(data, dtype=n.float32)
    for sid,stokes in enumerate(sclFactor): dout[sid,0,:]=data[sid,0,:].astype(n.float32)*sclFactor[sid]+offsets[sid]

    xvals=n.arange(dout.shape[2],dtype=n.float32)

    hdulist.close()

    if opts.rot: dout=n.roll(dout, dout.shape[2]/2, axis=2)

    #LOW PASS FILTER
    ntaps=dout.shape[2]
    cutoff=opts.cutoff
    fir=signal.firwin(ntaps,cutoff)
    ifilter=n.convolve(dout[0,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]
    qfilter=n.convolve(dout[1,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]
    ufilter=n.convolve(dout[2,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]
    vfilter=n.convolve(dout[3,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]

    if opts.save:
        dirname,basename=os.path.split(os.path.abspath(args[0]))
        outputname=basename.split('.fits')[0]+'.lpf.fits'
        outputname=dirname+'/'+outputname
        shutil.copy(os.path.abspath(args[0]),outputname)
        time.sleep(.1)

        hdulist=pf.open(outputname,mode='update')
        dwrite=n.zeros_like(dout)
        dwrite[0,0,:]=(ifilter-offsets[0])/sclFactor[0]
        dwrite[1,0,:]=(qfilter-offsets[1])/sclFactor[1]
        dwrite[2,0,:]=(ufilter-offsets[2])/sclFactor[2]
        dwrite[3,0,:]=(vfilter-offsets[3])/sclFactor[3]
        if opts.rot: dwrite=n.roll(dwrite, -dwrite.shape[2]/2, axis=2)
        dwrite=dwrite.flatten()
        dDict=hdulist[3].data
        dDict[0][-1]=dwrite
        hdulist[3].data=dDict
        hdulist.flush()
        hdulist.close()

    p.subplot(221)
    p.plot((ifilter-offsets[0])/sclFactor[0])
    p.plot((dout[0,0,:]-offsets[0])/sclFactor[0])
    
    p.subplot(222)
    p.plot((qfilter-offsets[1])/sclFactor[1])
    p.plot((dout[1,0,:]-offsets[1])/sclFactor[1])
    
    p.subplot(223)
    p.plot((ufilter-offsets[2])/sclFactor[2])
    p.plot((dout[2,0,:]-offsets[2])/sclFactor[2])

    p.subplot(224)
    p.plot((vfilter-offsets[3])/sclFactor[3])
    p.plot((dout[3,0,:]-offsets[3])/sclFactor[3])

    p.show()

