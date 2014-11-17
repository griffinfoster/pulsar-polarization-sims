#!/usr/bin/env python
"""
Convert a PSR FITS file to a headerless dat text file of Stokes' values
"""

import pyfits as pf
import numpy as n
import os
import sys

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [FITS file]')
    o.set_description(__doc__)
    o.add_option('-o', '--output', dest='output', default='template.smooth.dat',
        help='Output file name, default:template.smooth.dat')
    opts, args = o.parse_args(sys.argv[1:])

    #see www.atnf.csiro.au/research/pulsar/psrfists/fitsdef.html section: Subintegration data
    hdulist=pf.open(args[0])
    d=hdulist[3].data
    offsets=d[0][-3]
    sclFactor=d[0][-2]
    data=d[0][-1]
    if len(data.shape)==1:
        data.shape=(4,1,data.shape[-1]/4)

    dout=n.zeros_like(data, dtype=n.float32)
    for sid,stokes in enumerate(sclFactor): dout[sid,0,:]=data[sid,0,:].astype(n.float32)*sclFactor[sid]+offsets[sid]

    #n.savetxt(opts.output, dout, fmt='%10.10f')
    fh=open(opts.output,'w')
    outputstr=''
    for i in range(dout.shape[2]):
        outputstr+='%10.10f %10.10f %10.10f %10.10f\n'%(dout[0,0,i],dout[1,0,i],dout[2,0,i],dout[3,0,i])
    fh.write(outputstr)
    fh.close()

