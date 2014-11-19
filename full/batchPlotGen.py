#!/usr/bin/env python

import os, sys
import subprocess
import shutil

plotScript='/home/foster/pulsar/simulations/scripts/plotIXRvsdJ.py'
calMode=['cal','uncal']
mode=['rms','chi2','sigma']
rmsMode=[0,1,2]
imgExt=['png','pdf']
snrVal=100

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options]')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    for fn in args:
        for c in calMode:
            for m in mode:
                for r in rmsMode:
                    for ie in imgExt:
                        psrName=fn.split('/')[-1].split('.')[-2]
                        plotFn='%s.%s.%s.r%i.snr%i.%s'%(psrName,c,m,r,snrVal,ie)
                        cmd='%s -c %s -m %s -r %i --snr=%i -l -s %s %s'%(plotScript,c,m,r,snrVal,plotFn,os.path.abspath(fn))
                        print cmd
                        os.system(cmd)

