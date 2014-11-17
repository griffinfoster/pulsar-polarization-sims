#!/usr/bin/env python
"""Make a pickle of the TIM dictionary from either the text files or pickle files
"""

import os,sys
import numpy as np

import simsio
import cPickle as pkl

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [TIM files]')
    o.set_description(__doc__)
    o.add_option('-m','--mode',dest='mode',default='pkl',
        help='File mode: pkl or txt')
    opts, args = o.parse_args(sys.argv[1:])

    timTxtFiles=[]
    timPklFiles=[]
    for fn in args:
        if fn.endswith('txt'): timTxtFiles.append(fn)
        elif fn.endswith('pkl'): timPklFiles.append(fn)
    if opts.mode=='txt':
        timDict=simsio.timFromTxt(timTxtFiles)
    else:
        timDict=simsio.timFromPkl(timPklFiles)
    ofn='TIM.pkl'
    pkl.dump(timDict,open(ofn,'wb'))

