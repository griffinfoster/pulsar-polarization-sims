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

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
    x: the input signal 
    window_len: the dimension of the smoothing window; should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    flat window will produce a moving average smoothing.

    output: the smoothed signal
    
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also: 
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
        scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=n.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=n.ones(window_len,'d')
    else:
        w=eval('n.'+window+'(window_len)')

    y=n.convolve(w/w.sum(),s,mode='valid')
    #return y
    return y[(window_len/2-1):-(window_len/2+1)]
    
if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [FITS file]')
    o.set_description(__doc__)
    o.add_option('-r', '--rot', dest='rot', action='store_true',
        help='Rotate the profile by 0.5 of the phase')
    o.add_option('-w','--win_len',dest='win_len',default=11,type='int',
        help='Window smoothing size, should be odd, default:11')
    o.add_option('-s','--save',dest='save',action='store_true',
        help='Save the smoothed profile to a new fits file')
    o.add_option('-S','--shift',dest='shift',default=0, type='int',
        help='Shift the smoothed profile to the left N values default:0')
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

    ##LOW PASS FILTER
    #ntaps=dout.shape[2]
    #cutoff=opts.cutoff
    #fir=signal.firwin(ntaps,cutoff)
    #ifilter=n.convolve(dout[0,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]
    #qfilter=n.convolve(dout[1,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]
    #ufilter=n.convolve(dout[2,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]
    #vfilter=n.convolve(dout[3,0,:],fir)[int(ntaps/2)-1:-1*int(ntaps/2)]

    #SMOOTHING
    ifilter=smooth(dout[0,0,:],window_len=opts.win_len)
    qfilter=smooth(dout[1,0,:],window_len=opts.win_len)
    ufilter=smooth(dout[2,0,:],window_len=opts.win_len)
    vfilter=smooth(dout[3,0,:],window_len=opts.win_len)

    #SHIFTING
    if not (opts.shift==0):
        shift=-1*opts.shift
        print 'Applying a shift of %i units'%shift
        ifilter=n.roll(ifilter,shift)
        qfilter=n.roll(qfilter,shift)
        ufilter=n.roll(ufilter,shift)
        vfilter=n.roll(vfilter,shift)

    if opts.save:
        dirname,basename=os.path.split(os.path.abspath(args[0]))
        outputname=basename.split('.fits')[0]+'.smooth.fits'
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
        #dwrite=dwrite.flatten()
        dDict=hdulist[3].data
        print dwrite.shape
        dDict[0][-1]=dwrite
        hdulist[3].data=dDict
        hdulist.flush()
        hdulist.close()

    #p.subplot(221)
    #p.plot((ifilter-offsets[0])/sclFactor[0])
    #p.plot((dout[0,0,:]-offsets[0])/sclFactor[0])
    #
    #p.subplot(222)
    #p.plot((qfilter-offsets[1])/sclFactor[1])
    #p.plot((dout[1,0,:]-offsets[1])/sclFactor[1])
    #
    #p.subplot(223)
    #p.plot((ufilter-offsets[2])/sclFactor[2])
    #p.plot((dout[2,0,:]-offsets[2])/sclFactor[2])

    #p.subplot(224)
    #p.plot((vfilter-offsets[3])/sclFactor[3])
    #p.plot((dout[3,0,:]-offsets[3])/sclFactor[3])

    #p.show()

