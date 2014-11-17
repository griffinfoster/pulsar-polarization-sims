#!/bin/sh
TEMPLATEDIR=/home/griffin/pulsar/simulations/full
SCRIPTDIR=/home/griffin/pulsar/simulations/scripts
DICTSCRIPT=pklReduceDict.py

mkdir $TEMPLATEDIR/J1824-2452/results/analysis/
cd $TEMPLATEDIR/J1824-2452/results/analysis/
$SCRIPTDIR/$DICTSCRIPT  --mjd=../../J1824.simulate.mjd ../tim*.pkl ../rms*.txt -o analysisDict.pkl

