#!/usr/bin/env python
"""Convert an output tempo2 MJD file into a MJD list readable by beam2fits
tempo2 -gr fake -ndobs 1 -nobsd 1 -randha y -start 50000 -end 50500 -rms auto -f <par file>
parseTempo2fakeSim.py -f <simulated MJD file>
"""

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

    secPerDay=24*60*60.

    fh=open(opts.tempo2file,'r')
    fhOut=open(opts.tempo2file+'.mjd','w')
    outData=''
    #NOTE: odd bug with tempo2 fake that generates 10x the number of MJDs, so just take the first 503
    lcount=0
    for l in fh:
        if lcount > 502: continue
        lws= ' '.join(l.split())
        if len(lws.split(' '))<3: continue
        mjd=float(lws.split(' ')[2])
        imjd=int(mjd)
        smjd=(mjd-imjd)*secPerDay
        fmjd=smjd-int(smjd)
        smjd=int(smjd)
        outData+="       %i       %i    %.20f\n"%(imjd,smjd,fmjd)
        lcount+=1
    fh.close()
    fhOut.write(outData)
    fhOut.close()

