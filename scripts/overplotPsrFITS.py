#!/usr/bin/env python

import pyfits as pf
import numpy as n
import pylab as p
import os
import sys

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [FITS file]')
    o.set_description(__doc__)
    o.add_option('-m', '--mode', dest='mode', default='basic',
        help='Plotting mode: basic, simplified default: basic')
    o.add_option('-r', '--rot', dest='rot', action='store_true',
        help='Rotate the profile by 0.5 of the phase')
    o.add_option('-s', '--savefig', dest='savefig', default=None,
        help='Save figure as a PDF')
    o.add_option('-t', '--title', dest='title', default='Pulsar Profile',
        help='Plot title, default: Pulsar Profile')
    opts, args = o.parse_args(sys.argv[1:])
    
    hdulist=pf.open(args[0])
    hdulist0=pf.open(args[1])
    #print hdulist.info()

    primary=hdulist['PRIMARY'].header
    print primary['FITSTYPE']

    #see www.atnf.csiro.au/research/pulsar/psrfists/fitsdef.html section: Subintegration data
    d=hdulist[3].data
    d0=hdulist0[3].data
    print d
    offsets=d[0][-3]
    sclFactor=d[0][-2]
    data=d[0][-1]
    data0=d0[0][-1]
    print sclFactor
    print offsets
    print data.shape
    if len(data.shape)==1:
        data.shape=(4,1,data.shape[-1]/4)
        print data.shape

    dout=n.zeros_like(data, dtype=n.float32)
    dout0=n.zeros_like(data, dtype=n.float32)
    for sid,stokes in enumerate(sclFactor):
        dout[sid,0,:]=data[sid,0,:].astype(n.float32)*sclFactor[sid]+offsets[sid]
        dout0[sid,0,:]=data0[sid,0,:].astype(n.float32)*sclFactor[sid]+offsets[sid]

    xvals=n.arange(dout.shape[2],dtype=n.float32)
    #xvals to pulse phase
    xvals/=xvals.shape[0]
    #normalize amplitude
    #dout/=n.max(dout)

    if opts.rot:
        dout=n.roll(dout, dout.shape[2]/2, axis=2)
        dout0=n.roll(dout0, dout.shape[2]/2, axis=2)

    #p.plot(xvals, dout[0,0,:], 'k-.', label='I')
    #p.plot(xvals, dout[1,0,:], 'k--', label='Q')
    #p.plot(xvals, dout[2,0,:], 'k:', label='U')
    #p.plot(xvals, dout[3,0,:], 'k-', label='V')
    #p.plot(xvals, dout[1,0,:], 'g--', label='Q')
    #p.plot(xvals, dout[2,0,:], 'r:', label='U')
    #p.plot(xvals, dout[3,0,:], 'c-', label='V')
    p.plot(xvals, dout[0,0,:], 'k-', label='raw')
    p.plot(xvals, dout0[0,0,:], 'r--', label='smooth')

    #p.ylim(ymax=1.2)

    ax=p.gca()

    if opts.mode.startswith('bas'):
        #mode for viewing
        p.xlabel('Pulse Phase')
        p.ylabel('Amplitude')
        p.title(opts.title)
        p.legend()
    elif opts.mode.startswith('simp'):
        #mode for paper
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        p.text(0.05,0.93,opts.title,fontsize=25,verticalalignment='center',transform=ax.transAxes)
        fig=p.gcf()
        fig.set_size_inches(6.,6.)
    elif opts.mode.startswith('fill'):
        #hack to get empty pdf for latex
        fig=p.gcf()
        fig.set_size_inches(6.,6.)
        p.axis('off')

    #save figure
    if opts.savefig is None: p.show()
    else: p.savefig(opts.savefig, bbox_inches='tight')

    hdulist.close()

