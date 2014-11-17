#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import pyfits as pf
import numpy as n
import pylab as p
import os
import sys

def Jones2Mueller(J):
    A=n.matrix([[1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1j,-1j,0]])
    M=A*n.kron(J,J.conj())*n.linalg.inv(A)
    return n.real(M)

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] [FITS file]')
    o.set_description(__doc__)
    o.add_option('-r', '--rot', dest='rot', action='store_true',
        help='Rotate the profile by 0.5 of the phase')
    opts, args = o.parse_args(sys.argv[1:])

    fs=25
    
    hdulist=pf.open(args[0])
    #print hdulist.info()

    primary=hdulist['PRIMARY'].header
    print primary['FITSTYPE']

    #see www.atnf.csiro.au/research/pulsar/psrfists/fitsdef.html section: Subintegration data
    d=hdulist[3].data
    print d
    offsets=d[0][-3]
    sclFactor=d[0][-2]
    data=d[0][-1]
    print sclFactor
    print offsets
    print data.shape
    if len(data.shape)==1:
        data.shape=(4,1,data.shape[-1]/4)
        print data.shape

    dout=n.zeros_like(data, dtype=n.float32)
    for sid,stokes in enumerate(sclFactor): dout[sid,0,:]=data[sid,0,:].astype(n.float32)*sclFactor[sid]+offsets[sid]

    xvals=n.arange(dout.shape[2],dtype=n.float32)
    #xvals to pulse phase
    xvals/=xvals.shape[0]
    #normalize amplitude
    dout/=n.max(dout)

    if opts.rot: dout=n.roll(dout, dout.shape[2]/2, axis=2)

    IXR=7.
    polPur=1/n.sqrt(IXR)
    CondNum=(1.+polPur)/(1.-polPur)
    Stokes=n.matrix([dout[0,0,:],dout[1,0,:],dout[2,0,:],dout[3,0,:]])

    deltaJAmp=20./100.
    SN=100.

    u,s,v=n.linalg.svd(n.matrix(n.random.randn(2,2))+1j*n.matrix(n.random.randn(2,2)))
    d=n.matrix([[1.,0],[0,1./CondNum]])
    Jtrue=u*d*v.conj().transpose()
    Mtrue=Jones2Mueller(Jtrue)
    deltaJ=deltaJAmp*(n.matrix(n.random.randn(2,2))+1j*n.matrix(n.random.randn(2,2)))
    Mest=Jones2Mueller(Jtrue+deltaJ)

    #a: profile
    #p.subplot(221)
    StokesTrue=Mtrue*Stokes
    StokesTrue0=StokesTrue/n.max(StokesTrue)
    #p.plot(xvals, StokesTrue[0].T, 'k-.', label='I')
    #p.plot(xvals, StokesTrue[1].T, 'k--', label='Q')
    #p.plot(xvals, StokesTrue[2].T, 'k:', label='U')
    #p.plot(xvals, StokesTrue[3].T, 'k-', label='V')
    p.plot(xvals, StokesTrue[1].T, 'g--', label='Q')
    p.plot(xvals, StokesTrue[2].T, 'r:', label='U')
    p.plot(xvals, StokesTrue[3].T, 'c-', label='V')
    p.plot(xvals, StokesTrue[0].T, 'k-', label='I')
    p.ylim((-0.3,1.0))
    ax=p.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    p.text(0.05,0.93,'$\mathbf{M_{sys} e^S}$',fontsize=fs,verticalalignment='center',transform=ax.transAxes)

    fig=p.gcf()
    fig.set_size_inches(6.,6.)
    p.savefig('sim_profile.pdf', bbox_inches='tight')
    p.clf()

    #b: delta J
    StokesdJ=Mest*Stokes
    StokesdJ0=StokesdJ/n.max(StokesdJ)
    #p.subplot(222)
    #p.plot(xvals, StokesdJ[0].T, 'k-.', label='I')
    #p.plot(xvals, StokesdJ[1].T, 'k--', label='Q')
    #p.plot(xvals, StokesdJ[2].T, 'k:', label='U')
    #p.plot(xvals, StokesdJ[3].T, 'k-', label='V')
    p.plot(xvals, StokesdJ[1].T, 'g--', label='Q')
    p.plot(xvals, StokesdJ[2].T, 'r:', label='U')
    p.plot(xvals, StokesdJ[3].T, 'c-', label='V')
    p.plot(xvals, StokesdJ[0].T, 'k-', label='I')
    p.ylim((-0.3,1.0))
    ax=p.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    p.text(0.05,0.93,'$\mathbf{\^M_{sys} e^S}$',fontsize=fs,verticalalignment='center',transform=ax.transAxes)

    fig=p.gcf()
    fig.set_size_inches(6.,6.)
    p.savefig('sim_deltaj.pdf', bbox_inches='tight')
    p.clf()

    nrSamps=n.shape(Stokes)[1]

    #Create noise
    RecNoise=n.matrix(n.zeros((4,nrSamps)))
    for indJ in range(nrSamps):
        RecNoise[0,indJ]=n.random.randn()
        RecNoise[1,indJ]=n.random.randn()
        RecNoise[2,indJ]=n.random.randn()
        RecNoise[3,indJ]=n.random.randn()

    maxI=max(abs(Stokes[0,k]) for k in range(nrSamps))
    Stokes=(SN)*Stokes/maxI

    #Pass pulsar profile through interferometer
    StokesRaw=Mtrue*(Stokes)+RecNoise  #Compute raw signal

    #c: rec noise
    StokesRaw/=SN
    StokesRaw0=StokesRaw/n.max(StokesRaw)
    #p.subplot(223)
    #p.plot(xvals, StokesRaw[0].T, 'k-.', label='I')
    #p.plot(xvals, StokesRaw[1].T, 'k--', label='Q')
    #p.plot(xvals, StokesRaw[2].T, 'k:', label='U')
    #p.plot(xvals, StokesRaw[3].T, 'k-', label='V')
    p.plot(xvals, StokesRaw[1].T, 'g--', label='Q')
    p.plot(xvals, StokesRaw[2].T, 'r:', label='U')
    p.plot(xvals, StokesRaw[3].T, 'c-', label='V')
    p.plot(xvals, StokesRaw[0].T, 'k-', label='I')
    p.ylim((-0.3,1.0))
    ax=p.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    p.text(0.05,0.93,'$\mathbf{\^v_{obs}^S}$',fontsize=fs,verticalalignment='center',transform=ax.transAxes)

    fig=p.gcf()
    fig.set_size_inches(6.,6.)
    p.savefig('sim_snr.pdf', bbox_inches='tight')
    p.clf()

    Mest0=Mest/n.max(Mest)
    invMest=n.linalg.inv(Mest)
    print n.max(invMest)
    print n.max(StokesRaw)
    invMest/=n.max(invMest)
    StokesCalEst=invMest*StokesRaw #Compute calibrated signal
    StokesCalEst=n.real(StokesCalEst)

    #d: inverted M
    #StokesCalEst0=StokesCalEst/n.max(StokesCalEst)
    #p.subplot(224)
    #p.plot(xvals, StokesCalEst[0].T, 'k-.', label='I')
    #p.plot(xvals, StokesCalEst[1].T, 'k--', label='Q')
    #p.plot(xvals, StokesCalEst[2].T, 'k:', label='U')
    #p.plot(xvals, StokesCalEst[3].T, 'k-', label='V')
    p.plot(xvals, StokesCalEst[1].T, 'g--', label='Q')
    p.plot(xvals, StokesCalEst[2].T, 'r:', label='U')
    p.plot(xvals, StokesCalEst[3].T, 'c-', label='V')
    p.plot(xvals, StokesCalEst[0].T, 'k-', label='I')
    p.ylim((-0.3,1.0))
    ax=p.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    p.text(0.05,0.93,'$\mathbf{\^e^S}$',fontsize=fs,verticalalignment='center',transform=ax.transAxes)

    fig=p.gcf()
    fig.set_size_inches(6.,6.)
    p.savefig('sim_estimated.pdf', bbox_inches='tight')

    #e: original profile
    p.clf()
    Stokes0=Stokes/n.max(Stokes)
    p.plot(xvals, Stokes0[1].T, 'g--', label='Q')
    p.plot(xvals, Stokes0[2].T, 'r:', label='U')
    p.plot(xvals, Stokes0[3].T, 'c-', label='V')
    p.plot(xvals, Stokes0[0].T, 'k-', label='I')
    p.ylim((-0.3,1.0))
    ax=p.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    p.text(0.05,0.93,'$\mathbf{e^S}$',fontsize=fs,verticalalignment='center',transform=ax.transAxes)

    fig=p.gcf()
    fig.set_size_inches(6.,6.)
    p.savefig('sim_original.pdf', bbox_inches='tight')

    #p.show()

    hdulist.close()

