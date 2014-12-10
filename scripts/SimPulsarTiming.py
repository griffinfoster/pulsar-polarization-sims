#!/usr/bin/env python
"""Simulate Pulsar Timing ToA's using PSARCHIVE abd TEMPO2
"""
#TODO:include MJDs in tim files

import numpy as np
import os, sys
import subprocess
import shutil
import glob
import cPickle as pickle
import time
import math

invStokes=np.matrix([[1.,0.,0.,1.],[1.,0.,0.,-1.],[0.,1.,1.,0.],[0.,-1.j,1.j,0.]])
Stokes=.5*np.matrix([[1.,1.,0.,0.],[0.,0.,1.,1.j],[0.,0.,1.,-1.j],[1.,-1.,0.,0.]])   #note 1/2

#Get environiment variables for unset/export hacks when running pat and tempo2
#TEMPO_VAR=os.environ['TEMPO']
#TEMPO2_VAR=os.environ['TEMPO2']

def Jones2Mueller(J):
    """Convert a Jones Matrix to a Mueller Matrix"""
    A=np.matrix([[1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,1j,-1j,0]])
    M=invStokes*np.kron(J,J.conj())*Stokes
    return np.real(M)

def readTimFile(fn):
    """Read a TIM output file from PAT, and return a numpy array of the contents"""
    fh=open(fn,'r')
    timData=fh.read()
    fh.close()
    lines=timData.split('\n')
    lines=lines[1:-1]
    arr=[]
    for l in lines:
        if l.startswith('FORMAT'): continue
        splitLine=l.split()
        arr.append([float(splitLine[2]),float(splitLine[3])])
    arr=np.array(arr)
    return arr

def StokesToTextFile(data,outfile):
    """Write a Stokes spectrum to text file"""
    filehandle=open(outfile,'w')
    for i in range(np.shape(data)[0]):
        line=''
        for j in range(np.shape(data)[1]):
            x=data[i,j]
            line=line+repr(x)+' '
        filehandle.write(line+'\n')
    filehandle.close()

def StokesFromTextFile(infile):
    """Return a Stokes spectrum from a text file"""
    I=[]
    Q=[]
    U=[]
    V=[]
    filehandle=open(infile,'r')
    for line in filehandle.readlines():
        l=line.split()
        if l!=[]:
            I.append(float(l[0]))
            Q.append(float(l[1]))
            U.append(float(l[2]))
            V.append(float(l[3]))
	filehandle.close()
    return I,Q,U,V

def SimInterfMeasPuls(Stokes,ofnPrefix,SN,nomPolPur,deltaJAmp): 
    """Simulate an interferometer measurement of the pulsar profile Stokes spectrum, and writes spectrum to text file
    inputs:
        Stokes: template Stokes specturm
        ofnPrefix: output filename prefix, two files are written ofnPrefix+'.cal.dat' and ofnPrefix+'.uncal.dat'
        SN: nominal signal-to-noise ratio
        nomPolPur: polarization purity, which is related to the IXR
        deltaJAmp: calibration error, 1=100%
    """
    #Generate the polarimeter matrices for interferometer
    #nomPolPur is the polarization purity or 1/sqrt(IXR) (gmin is 1)
    #    nomPolCond=(1.0+nomPolPur)/(1.0-nomPolPur) 
    nomPolCond=(1.+nomPolPur)/(1.-nomPolPur)
    
    U,s,Vh=np.linalg.svd(np.matrix(np.random.randn(2,2))+1j*np.matrix(np.random.randn(2,2))) 
    d=np.matrix([[1.,0],[0,1./nomPolCond]])
    Jtrue=np.dot(U,np.dot(d,Vh))
    Mtrue=Jones2Mueller(Jtrue)

    deltaJ=deltaJAmp*(np.matrix(np.random.randn(2,2))+1j*np.matrix(np.random.randn(2,2)))
    Mest=Jones2Mueller(Jtrue+deltaJ)

    #Create noise
    spectrumLen=np.shape(Stokes)[1]
    RecNoise=np.matrix(np.zeros((4,spectrumLen)))
    for indJ in range(spectrumLen):
        RecNoise[0,indJ]=np.random.randn()
        RecNoise[1,indJ]=np.random.randn()
        RecNoise[2,indJ]=np.random.randn()
        RecNoise[3,indJ]=np.random.randn()

    #Scale data to SNR before adding receiver noise
    maxI=np.max(np.abs(Stokes[0]))
    Stokes=(SN)*Stokes/maxI 

    StokesRaw=Mtrue*(Stokes)+RecNoise  #Compute raw signal
    StokesCalEst=np.linalg.inv(Mest)*StokesRaw #Compute calibrated signal
    StokesCalEst=np.real(StokesCalEst)

    StokesToTextFile(StokesCalEst.transpose(),ofnPrefix+'.cal.dat')
    StokesToTextFile(StokesRaw.transpose(),ofnPrefix+'.uncal.dat')

def TimingSimulation(SN_Values,nomPolPur_Values,deltaJAmp_Values,pd1,parFile,mjdFile,rootDir,b2f,timingModes):
    """Simulate the timing of a pulsar across multiple parameter spaces
    SN_Values: signal to noise values to loop over
    nomPolPur_Values: polarization purity values to loop over
    deltaJAmp_Values: calibration error values to loop over
    pd1: the period of the pulsar you are processing, so that the rms values are in reasonable units
    parFile: pulsar PAR file
    mjdFile: file with a list of MJDs of expected pulse ToA
    rootDir: working root directory
    b2f: beam2fits executable name
    timingModes: table of methods and calibration mode to use
    """
    fr=1.400 # fr is the centre frequency in GHz (not critical)

    rmsDict={}
    prefix='simPsr'
    verbose=False
    FNULL=open(os.devnull, 'w')

    #create the template profile
    fhMJD=open(mjdFile,'r')
    mjdline0=fhMJD.readline() #use the first MJD of the file for the template
    fhMJD.close()
    # Creating a fits file from the ascii profile which is template.dat
    if os.path.exists(os.path.abspath('.')+'/template.fits'): os.remove(os.path.abspath('.')+'/template.fits')
    #use original 'noisy' profile as the template for pat, but use the smoothed profile as the starting profile for simulation
    #command='/home/griffin/pulsar/PSRBeam/beam2fits%i template.dat template.hdr template.fits %f %f %s'%(b2f,pd1,fr,mjdline0)
    command='/home/foster/pulsar/PSRBeam/beam2fits%i template.dat template.hdr template.fits %f %f %s'%(b2f,pd1,fr,mjdline0)
    print command
    os.system(command)
    # Creates the invariant interval for timing and set site to Parkes
    command='pam -IF -e ii --site 7 template.fits'
    print command
    os.system(command)
    # modifies the site in place for template.fits
    command='pam -m --site 7 template.fits'
    print command
    os.system(command)
    # fixes the fits file to work with the timing software
    command='psredit -c polc=1 -m template.fits'
    print command
    os.system(command)

    #MJD files have 3 values per line: DAYS SECONDS FRACTION_OF_SECOND, MJD=DAYS+(SECONDS+FRACTION_OF_SECOND)/86400.
    fhMJD=open(mjdFile,'r')
    mjdTriplet=fhMJD.readlines()
    nobs=len(mjdTriplet)
    fhMJD.close()
                
    compt=0
    nsims=len(SN_Values)*len(deltaJAmp_Values)*len(nomPolPur_Values)
    for SN in SN_Values:
        for deltaJAmp in deltaJAmp_Values:
            for nomPolPur in nomPolPur_Values:
                print 'Simulating Stokes spectra...'
                print 'SNR: %f \t dJ: %f \t PolPur: %f'%(SN,deltaJAmp,nomPolPur)
	            #I,Q,U,V=StokesFromTextFile('template.dat')
	            #I,Q,U,V=StokesFromTextFile('template.lpf.dat'
                I,Q,U,V=StokesFromTextFile('template.smooth.dat')
                Stokes=np.matrix([I,Q,U,V])
                for i in range(nobs):
    	            ofnPrefix='%s%04i'%(prefix,i)
    	            SimInterfMeasPuls(Stokes,ofnPrefix,SN,nomPolPur,deltaJAmp)
                print 'done'
                
                #So now there should be nobs prefix*.cal.dat and prefix*.uncal.dat files ; convert to fits:
                print 'Converting Stokes spectrum to FITS files'
                for datFile in sorted(glob.glob('%s*.dat'%prefix)):
                    if verbose: print 'Converting '+datFile+' to fits'
                    idx=int(datFile.split('.')[0].split(prefix)[-1])
                    mjdline=mjdTriplet[idx]
                    fitsFile=datFile.split('.dat')[0]+'.fits'
                    
                    #command='/home/griffin/pulsar/PSRBeam/beam2fits%i %s template.hdr %s %f %f %s'%(b2f,datFile,fitsFile,pd1,fr,mjdline)
                    command='/home/foster/pulsar/PSRBeam/beam2fits%i %s template.hdr %s %f %f %s'%(b2f,datFile,fitsFile,pd1,fr,mjdline)
                    if verbose: print command
                    return_code=subprocess.call(command, shell=True, stdout=FNULL)
                    
                    command='pam -IF -e ii --site 7 %s > /dev/null'%fitsFile
                    if verbose: print command
                    os.system(command)
                    
                    command='pam -m --site 7 %s > /dev/null'%fitsFile
                    if verbose: print command
                    os.system(command)
                    
                    command='rm -f %s'%datFile
                    if verbose: print command
                    os.system(command)
                print 'done'

                #Now we are left with nobs FITS files; time them in one go
                for mid,mode in enumerate(['cal','uncal']):

                    #Total Intensity
                    if timingModes['ti'][mid] is True:
                        tim0TxtFile='tim0_sn%i_dj%.4f_pp%.4f.%s.txt'%(int(SN),deltaJAmp,nomPolPur,mode)
                        tim0PklFile='tim0_sn%i_dj%.4f_pp%.4f.%s.pkl'%(int(SN),deltaJAmp,nomPolPur,mode)
                        rms0TxtFile='rms0_sn%i_dj%.4f_pp%.4f.%s.txt'%(int(SN),deltaJAmp,nomPolPur,mode)
                        command='pat -f tempo2 -F -s template.fits %s*.%s.fits > %s'%(prefix,mode,tim0TxtFile)
                        print command
                        os.system(command)
                        
                        command="tempo2 -f %s %s |grep \"residual\" |awk '{print $11}' > %s"%(parFile,tim0TxtFile,rms0TxtFile)
                        print command
                        os.system(command)

                        tim0Data=readTimFile(tim0TxtFile)
                        tim0pkl=open(tim0PklFile,'wb')
                        pickle.dump(tim0Data,tim0pkl)
                        tim0pkl.close()

                        #RMS values
                        fh=open('%s/%s'%(rootDir,rms0TxtFile),'r')
                        try: rms0=float(fh.readline())
                        #except ValueError: rms0=float('NaN')
                        except ValueError: rms0=-1.
                        fh.close()
                    else: rms0=-1
                    
                    #Invariant Interval
                    if timingModes['ti'][mid] is True:
                        tim1TxtFile='tim1_sn%i_dj%.4f_pp%.4f.%s.txt'%(int(SN),deltaJAmp,nomPolPur,mode)
                        tim1PklFile='tim1_sn%i_dj%.4f_pp%.4f.%s.pkl'%(int(SN),deltaJAmp,nomPolPur,mode)
                        rms1TxtFile='rms1_sn%i_dj%.4f_pp%.4f.%s.txt'%(int(SN),deltaJAmp,nomPolPur,mode)
                        command='pat -f tempo2 -F -s template.ii %s*.%s.ii > %s'%(prefix,mode,tim1TxtFile)
                        print command
                        os.system(command)
                        
                        command="tempo2 -f %s %s |grep \"residual\" |awk '{print $11}' > %s"%(parFile,tim1TxtFile,rms1TxtFile)
                        print command
                        os.system(command)

                        tim1Data=readTimFile(tim1TxtFile)
                        tim1pkl=open(tim1PklFile,'wb')
                        pickle.dump(tim1Data,tim1pkl)
                        tim1pkl.close()
                    
                        #RMS values
                        fh=open('%s/%s'%(rootDir,rms1TxtFile),'r')
                        try: rms1=float(fh.readline())
                        #except ValueError: rms1=float('NaN')
                        except ValueError: rms1=-1.
                        fh.close()
                    else: rms1=-1
                    
                    #Matrix Template Matching
                    if timingModes['mtm'][mid] is True:
                        tim2TxtFile='tim2_sn%i_dj%.4f_pp%.4f.%s.txt'%(int(SN),deltaJAmp,nomPolPur,mode)
                        tim2PklFile='tim2_sn%i_dj%.4f_pp%.4f.%s.pkl'%(int(SN),deltaJAmp,nomPolPur,mode)
                        rms2TxtFile='rms2_sn%i_dj%.4f_pp%.4f.%s.txt'%(int(SN),deltaJAmp,nomPolPur,mode)

                        tim2str='FORMAT 1\n'
                        tempo_tim2str='FORMAT 1\n'
                        tempo_tim2TxtFile='tim2_sn%i_dj%.4f_pp%.4f.%s.tempo'%(int(SN),deltaJAmp,nomPolPur,mode)
                        solutionCtr=0
                        for ff in sorted(glob.glob('%s*.%s.fits'%(prefix,mode))):
                            command='pat -p -f tempo2 -F -s template.fits %s 2> mtm.log'%(ff)
                            cmdOutput=subprocess.check_output(command, shell=True)
                            lines=cmdOutput.split('\n')
                            if len(lines)==2: #no solution
                                tim2str+='no_solution 0 0 0 0 \n'
                                #TODO: add fits file to list of files to replace
                                #   run pat on new fits files
                                #   loop until all samples have a solution
                            else:
                                tim2str+=lines[1]+'\n'
                                tempo_tim2str+=lines[1]+'\n'
                                solutionCtr+=1
                        print '%s: MTM returned %i of %i (%f%%) solutions'%(mode,solutionCtr,nobs,100.*float(solutionCtr)/nobs)
                        #write MTM string to file
                        fh=open(tim2TxtFile,'w')
                        fh.write(tim2str)
                        fh.close()
                        #write to the file used to run tempo
                        fh=open(tempo_tim2TxtFile,'w')
                        fh.write(tempo_tim2str)
                        fh.close()
                        
                        #command="tempo2 -f %s %s |grep \"residual\" |awk '{print $11}' > %s"%(parFile,tim2TxtFile,rms2TxtFile)
                        command="tempo2 -f %s %s |grep \"residual\" |awk '{print $11}' > %s"%(parFile,tempo_tim2TxtFile,rms2TxtFile)
                        print command
                        os.system(command)
                        
                        tim2Data=readTimFile(tim2TxtFile)
                        tim2pkl=open(tim2PklFile,'wb')
                        pickle.dump(tim2Data,tim2pkl)
                        tim2pkl.close()
               
                        fh=open('%s/%s'%(rootDir,rms2TxtFile),'r')
                        try: rms2=float(fh.readline())
                        #except ValueError: rms2=float('NaN').
                        except ValueError: rms2=-1.
                        fh.close()
                    else: rms2=-1
                    
                    rmsDict[(mode,SN,deltaJAmp,nomPolPur)]=(rms0,rms1,rms2)
                    print 'RMS(%s): total intersity: %f \t invariant interval: %f \t matrix template matching: %f'%(mode,rms0,rms1,rms2)

                #save rmsDict to a pickle file after ever step of the simulation
                output=open('RMS.pkl','wb')
                pickle.dump(rmsDict,output)
                output.close()

                #clean up simulation files
                command='rm -f %s*'%(prefix)
                print command
                os.system(command)

                compt+=1
                perc=100.*compt/nsims
                print 'Precent complete: %f'%(perc)
                print 'Progress: %i of %i'%(compt,nsims)

	# Writing data in file
	fh=open('%s/RMS.dat'%rootDir,'w')
	fh.write('mode S/N deltaJ PolPur rms0 rms1 rms2\n')
	for k in rmsDict.keys():
		fh.write(str(k[0])+' '+str(k[1])+' '+str(k[2])+' '+str(k[3])+' '+str(rmsDict[k][0])+' '+str(rmsDict[k][1])+' '+str(rmsDict[k][2])+'\n')
	fh.close()

    #move output files to the results directory
    resultsDir='%s/results'%rootDir
    for fn in glob.glob(os.path.join(rootDir, 'tim*')): shutil.move(fn, resultsDir)
    for fn in glob.glob(os.path.join(rootDir, 'rms*')): shutil.move(fn, resultsDir)
    for fn in glob.glob(os.path.join(rootDir, 'RMS*')): shutil.move(fn, resultsDir)

if __name__ == "__main__":
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options]')
    o.set_description(__doc__)
    o.add_option('-d', '--dir', dest='directory', default=None,
        help='Working directory, no default')
    o.add_option('-m','--mjd',dest='mjd',default=None,
        help='MJD arrival times, no default')
    o.add_option('-p','--par',dest='par',default=None,
        help='pulsar PAR file, no default')
    o.add_option('--pd1', dest='pd1', default=0.005, type=float,
        help='Pulsar period in seconds, default=0.005')
    o.add_option('--beam2fits', dest='beam2fits', default=512, type=int,
        help='Select which beam2fits compile to use: 256, 512, 1024, 2048 default: 512')
    o.add_option('-s','--sim',dest='simParam',default=None,
        help='Use a simulation parameter file, default is to use the hardcoded values in the script')
    o.add_option('-M','--mode',dest='mode',default=None,
        help='Timing modes to use on calibrated or uncalibrated data, default: all i.e. ti_uncal,ti_cal,ii_uncal,ii_cal,mtm_uncal,mtm_cal, ')
    opts, args = o.parse_args(sys.argv[1:])

    startTime=time.time()
    
    if opts.directory==None:
        print 'Working directory not set'
        exit(1)
    elif opts.mjd==None:
        print 'MJD file not set'
        exit(1)
    elif opts.par==None:
        print 'PAR file not set'
        exit(1)
    directory=os.path.abspath(opts.directory)+'/'
    parFile=os.path.abspath(opts.par)
    mjdFile=os.path.abspath(opts.mjd)
    print 'Using %s as the working directory'%(directory)
    print 'Using PAR file %s'%(parFile)
    print 'Using MJD file %s'%(mjdFile)

    #Setup results directory
    if os.path.exists(os.path.abspath('.')+'/results'): shutil.rmtree(os.path.abspath('.')+'/results')
    os.makedirs(os.path.abspath('.')+'/results')

    command='rm -f template.fits template.ii'
    print command
    os.system(command)

    if opts.simParam is None:
        print 'Using hardcoded parameters'
        #Signal to Noise
        #setup such that tInt=1 is equivalent to an SNR=1000, SNR scales as the sqrt of the integration time
        tInt=np.array([0.001,0.01,0.1,1.])
        #tInt=np.array([0.1])
        print 'Integration Time:',tInt
        snRng=1000.*np.sqrt(tInt)
        print 'SNR Values:',snRng

        #Polarization purity
        # polPurRng=1/sqrt(IXR)
        IXRdbRng=[0.001,1.,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,35.,40.]
        IXRdbRng=np.array(IXRdbRng)
        IXRRng=10.**(IXRdbRng/10.)
        print 'IXR (dB):',IXRdbRng
        polPurRng=1/np.sqrt(IXRRng)
        print 'Polarization Purity:',polPurRng

        #Calibration Error
        dJRng=np.array([0.,1.,2.,3.,5.,7.,11.,15.,23.,30.])/100.
        #dJRng=np.array([7.])/100.
        print 'Calibration Error (%)',dJRng
    else:
        print 'Loading parameters from %s file'%opts.simParam
        fh=open(opts.simParam,'r')
        for line in fh.readlines():
            if line.startswith('#'): continue
            elif line.startswith('snr'): tInt=np.array(map(float,line.split(':')[-1].split(',')))
            elif line.startswith('ixr'): IXRdbRng=np.array(map(float,line.split(':')[-1].split(',')))
            elif line.startswith('dj'): dJRng=np.array(map(float,line.split(':')[-1].split(',')))/100.
        fh.close()
        #Signal to Noise
        #setup such that tInt=1 is equivalent to an SNR=1000, SNR scales as the sqrt of the integration time
        print 'Integration Time:',tInt
        snRng=1000.*np.sqrt(tInt)
        print 'SNR Values:',snRng

        #Polarization purity
        # polPurRng=1/sqrt(IXR)
        IXRRng=10.**(IXRdbRng/10.)
        print 'IXR (dB):',IXRdbRng
        polPurRng=1/np.sqrt(IXRRng)
        print 'Polarization Purity:',polPurRng

        #Calibration Error
        #dJRng=np.array([7.])/100.
        print 'Calibration Error (%)',dJRng

    print 'Running %i simulations'%(len(dJRng)*len(polPurRng)*len(snRng))

    if opts.mode is None:
        timingModes={'ti': [True, True], 'ii': [True, True], 'mtm': [True, True]}
    else:
        timingModes={'ti': [False, False], 'ii': [False, False], 'mtm': [False, False]}
        for mm in opts.mode.split(','):
            mpair=mm.split('_')
            if mpair[1].startswith('cal'): timingModes[mpair[0]][0]=True
            if mpair[1].startswith('uncal'): timingModes[mpair[0]][1]=True
    print 'Timing modes:', timingModes

    config={
        'snRng'     :   snRng,
        'polPurRng' :   polPurRng,
        'dJRng'     :   dJRng,
        'pd1'       :   opts.pd1,
        'parFile'   :   parFile,
        'mjdFile'   :   mjdFile,
        'beam2fits' :   opts.beam2fits,
        'timingModes' : timingModes,
    }

    #Run the simulation
    os.chdir(directory)
    TimingSimulation(config['snRng'],config['polPurRng'],config['dJRng'],config['pd1'],config['parFile'],config['mjdFile'],directory,config['beam2fits'],config['timingModes'])

    endTime=time.time()
    diffTime=endTime-startTime
    days=math.floor(diffTime/(24.*60.*60.))
    diffTime-=(24.*60.*60.)*days
    hours=math.floor(diffTime/(60.*60.))
    diffTime-=(60.*60.)*hours
    mins=math.floor(diffTime/60.)
    diffTime-=60.*mins
    secs=diffTime
    print  'Total Time: %03i:%02i:%02i:%2.4f'%(days,hours,mins,secs)

