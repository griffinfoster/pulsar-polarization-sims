#!/usr/bin/env python

import os,sys

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options]')
    o.set_description(__doc__)
    o.add_option('-f', '--file', dest='tempo2file', default=None,
        help='tempo2 fake mjd file, no default')
    opts, args = o.parse_args(sys.argv[1:])

    if opts.tempo2file==None:
        print 'tempo2 file not set'
        exit(1)

    secPerDay=24*60*60

    fh=open(opts.tempo2file,'r')
    fhOut=open(opts.tempo2file+'.mjd','w')
    outData=''
    for l in fh:
        lws= ' '.join(l.split())
        if len(lws.split(' '))<3: continue
        mjd=float(lws.split(' ')[2])
        imjd=int(mjd)
        smjd=(mjd-imjd)*secPerDay
        fmjd=smjd-int(smjd)
        smjd=int(smjd)
        outData+="       %i       %i    %f\n"%(imjd,smjd,fmjd)
    fh.close()
    fhOut.write(outData)
    fhOut.close()

