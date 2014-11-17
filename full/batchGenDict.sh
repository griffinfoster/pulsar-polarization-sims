#!/bin/sh
TEMPLATEDIR=/home/foster/pulsar/simulations/full
SCRIPTDIR=/home/foster/pulsar/simulations/scripts
DICTSCRIPT=pklReduceDict.py

echo J0613-0200
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J0613-0200/J0613.simulate.mjd $TEMPLATEDIR/J0613-0200/results/tim*.pkl $TEMPLATEDIR/J0613-0200/results/rms*.txt -o analysisDict.J0613.pkl

echo J1045-4509
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1045-4509/J1045.simulate.mjd $TEMPLATEDIR/J1045-4509/results/tim*.pkl $TEMPLATEDIR/J1045-4509/results/rms*.txt -o analysisDict.J1045.pkl

echo J1824-2452
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1824-2452/J1824.simulate.mjd $TEMPLATEDIR/J1824-2452/results/tim*.pkl $TEMPLATEDIR/J1824-2452/results/rms*.txt -o analysisDict.J1824.pkl

echo J1939+2134
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1939+2134/J1939.simulate.mjd $TEMPLATEDIR/J1939+2134/results/tim*.pkl $TEMPLATEDIR/J1939+2134/results/rms*.txt -o analysisDict.J1939.pkl

echo J1600-3053
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1600-3053/J1600.simulate.mjd $TEMPLATEDIR/J1600-3053/results/tim*.pkl $TEMPLATEDIR/J1600-3053/results/rms*.txt -o analysisDict.J1600.pkl

echo J1643-1224
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1643-1224/J1643.simulate.mjd $TEMPLATEDIR/J1643-1224/results/tim*.pkl $TEMPLATEDIR/J1643-1224/results/rms*.txt -o analysisDict.J1643.pkl

echo J1732-5049
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1732-5049/J1732.simulate.mjd $TEMPLATEDIR/J1732-5049/results/tim*.pkl $TEMPLATEDIR/J1732-5049/results/rms*.txt -o analysisDict.J1732.pkl

echo J1909-3744
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1909-3744/J1909.simulate.mjd $TEMPLATEDIR/J1909-3744/results/tim*.pkl $TEMPLATEDIR/J1909-3744/results/rms*.txt -o analysisDict.J1909.pkl

echo J2129-5721
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J2129-5721/J2129.simulate.mjd $TEMPLATEDIR/J2129-5721/results/tim*.pkl $TEMPLATEDIR/J2129-5721/results/rms*.txt -o analysisDict.J2129.pkl

echo J0437-4715
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J0437-4715/J0437.simulate.mjd $TEMPLATEDIR/J0437-4715/results/tim*.pkl $TEMPLATEDIR/J0437-4715/results/rms*.txt -o analysisDict.J0437.pkl

echo J0711-6830
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J0711-6830/J0711.simulate.mjd $TEMPLATEDIR/J0711-6830/results/tim*.pkl $TEMPLATEDIR/J0711-6830/results/rms*.txt -o analysisDict.J0711.pkl

echo J1024-0719
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1024-0719/J1024.simulate.mjd $TEMPLATEDIR/J1024-0719/results/tim*.pkl $TEMPLATEDIR/J1024-0719/results/rms*.txt -o analysisDict.J1024.pkl

echo J1603-7202
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1603-7202/J1603.simulate.mjd $TEMPLATEDIR/J1603-7202/results/tim*.pkl $TEMPLATEDIR/J1603-7202/results/rms*.txt -o analysisDict.J1603.pkl

echo J1713+0747
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1713+0747/J1713.simulate.mjd $TEMPLATEDIR/J1713+0747/results/tim*.pkl $TEMPLATEDIR/J1713+0747/results/rms*.txt -o analysisDict.J1713.pkl

echo J1730-2304
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1730-2304/J1730.simulate.mjd $TEMPLATEDIR/J1730-2304/results/tim*.pkl $TEMPLATEDIR/J1730-2304/results/rms*.txt -o analysisDict.J1730.pkl

echo J1744-1134
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1744-1134/J1744.simulate.mjd $TEMPLATEDIR/J1744-1134/results/tim*.pkl $TEMPLATEDIR/J1744-1134/results/rms*.txt -o analysisDict.J1744.pkl

echo J1857+0943
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1857+0943/J1857.simulate.mjd $TEMPLATEDIR/J1857+0943/results/tim*.pkl $TEMPLATEDIR/J1857+0943/results/rms*.txt -o analysisDict.J1857.pkl

echo J2124-3358
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J2124-3358/J2124.simulate.mjd $TEMPLATEDIR/J2124-3358/results/tim*.pkl $TEMPLATEDIR/J2124-3358/results/rms*.txt -o analysisDict.J2124.pkl

echo J1022+1001
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J1022+1001/J1022.simulate.mjd $TEMPLATEDIR/J1022+1001/results/tim*.pkl $TEMPLATEDIR/J1022+1001/results/rms*.txt -o analysisDict.J1022.pkl

echo J2145-0750
$SCRIPTDIR/$DICTSCRIPT --mjd=$TEMPLATEDIR/J2145-0750/J2145.simulate.mjd $TEMPLATEDIR/J2145-0750/results/tim*.pkl $TEMPLATEDIR/J2145-0750/results/rms*.txt -o analysisDict.J2145.pkl

