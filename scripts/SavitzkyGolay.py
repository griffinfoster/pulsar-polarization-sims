#!/usr/bin/env python
"""
Apply a low pass filter to a pulsar profile
"""

import pyfits as pf
import numpy as n
import pylab as p
import os
import sys
import shutil
import time

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    from math import factorial

    try:
        window_size = n.abs(n.int(window_size))
        order = n.abs(n.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = n.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = n.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - n.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + n.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = n.concatenate((firstvals, y, lastvals))
    return n.convolve( m[::-1], y, mode='valid')
    
if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [FITS file]')
    o.set_description(__doc__)
    o.add_option('-r', '--rot', dest='rot', action='store_true',
        help='Rotate the profile by 0.5 of the phase')
    o.add_option('-w','--win_len',dest='win_len',default=51,type='int',
        help='Window smoothing size, should be odd, default:51')
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
        print data.shape

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
    ifilter=savitzky_golay(dout[0,0,:], opts.win_len, 10)
    qfilter=savitzky_golay(dout[1,0,:], opts.win_len, 10)
    ufilter=savitzky_golay(dout[2,0,:], opts.win_len, 10)
    vfilter=savitzky_golay(dout[3,0,:], opts.win_len, 10)

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

