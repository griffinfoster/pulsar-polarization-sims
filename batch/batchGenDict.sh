#!/bin/sh
TEMPLATEDIR=/home/griffin/pulsar/simulations/full
SCRIPTDIR=/home/griffin/pulsar/simulations/scripts
DICTSCRIPT=pklReduceDict.py

mkdir $TEMPLATEDIR/J0613-0200/results/analysis/
cd $TEMPLATEDIR/J0613-0200/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J0613.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1045-4509/results/analysis/
cd $TEMPLATEDIR/J1045-4509/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1045.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1824-2452/results/analysis/
cd $TEMPLATEDIR/J1824-2452/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1824.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1939+2134/results/analysis/
cd $TEMPLATEDIR/J1939+2134/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1939.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1600-3053/results/analysis/
cd $TEMPLATEDIR/J1600-3053/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1600.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1643-1224/results/analysis/
cd $TEMPLATEDIR/J1643-1224/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1643.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1732-5049/results/analysis/
cd $TEMPLATEDIR/J1732-5049/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1732.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1909-3744/results/analysis/
cd $TEMPLATEDIR/J1909-3744/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1909.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J2129-5721/results/analysis/
cd $TEMPLATEDIR/J2129-5721/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J2129.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J0437-4715/results/analysis/
cd $TEMPLATEDIR/J0437-4715/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J0437.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J0711-6830/results/analysis/
cd $TEMPLATEDIR/J0711-6830/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J0711.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1024-0719/results/analysis/
cd $TEMPLATEDIR/J1024-0719/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1024.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1603-7202/results/analysis/
cd $TEMPLATEDIR/J1603-7202/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1603.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1713+0747/results/analysis/
cd $TEMPLATEDIR/J1713+0747/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1713.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1730-2304/results/analysis/
cd $TEMPLATEDIR/J1730-2304/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1730.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1744-1134/results/analysis/
cd $TEMPLATEDIR/J1744-1134/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1744.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1857+0943/results/analysis/
cd $TEMPLATEDIR/J1857+0943/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1857.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J2124-3358/results/analysis/
cd $TEMPLATEDIR/J2124-3358/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J2124.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J1022+1001/results/analysis/
cd $TEMPLATEDIR/J1022+1001/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1022.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

mkdir $TEMPLATEDIR/J2145-0750/results/analysis/
cd $TEMPLATEDIR/J2145-0750/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J2145.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

