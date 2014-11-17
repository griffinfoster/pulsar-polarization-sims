#!/usr/bin/env python
"""Make a pickle of the simulated expected MJD times
"""

import os,sys
import numpy as np

import simsio
import cPickle as pkl

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [MJD files]')
    o.set_description(__doc__)
    opts, args = o.parse_args(sys.argv[1:])

    for fn in args:
        mjds=simsio.mjdsFromTxt(fn)
        ofn=fn+'.pkl'
        print 'Writing %s'%ofn
        pkl.dump(mjds,open(ofn,'wb'))

