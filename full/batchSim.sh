#!/bin/sh
TEMPLATEDIR=/home/foster/pulsar/simulations/full
SCRIPTDIR=/home/foster/pulsar/simulations/scripts
SIMSCRIPT=SimPulsarTiming.py
PARAMFILE=full.params

#256
cd $TEMPLATEDIR/J1824-2452
$SCRIPTDIR/$SIMSCRIPT -d . -m J1824.simulate.mjd -p J1824.par --beam2fits 256  --pd1=0.003054 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1824.log &

cd $TEMPLATEDIR/J1939+2134
$SCRIPTDIR/$SIMSCRIPT -d . -m J1939.simulate.mjd -p J1939.par --beam2fits 256  --pd1=0.001557 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1939.log &

#512
cd $TEMPLATEDIR/J0613-0200
$SCRIPTDIR/$SIMSCRIPT -d . -m J0613.simulate.mjd -p J0613.par --beam2fits 512  --pd1=0.00306  -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J0613.log &

cd $TEMPLATEDIR/J1045-4509
$SCRIPTDIR/$SIMSCRIPT -d . -m J1045.simulate.mjd -p J1045.par --beam2fits 512  --pd1=0.00747  -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1045.log &

cd $TEMPLATEDIR/J1600-3053
$SCRIPTDIR/$SIMSCRIPT -d . -m J1600.simulate.mjd -p J1600.par --beam2fits 512  --pd1=0.00359  -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1600.log &

cd $TEMPLATEDIR/J1643-1224
$SCRIPTDIR/$SIMSCRIPT -d . -m J1643.simulate.mjd -p J1643.par --beam2fits 512  --pd1=0.00462  -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1643.log &

cd $TEMPLATEDIR/J1732-5049
$SCRIPTDIR/$SIMSCRIPT -d . -m J1732.simulate.mjd -p J1732.par --beam2fits 512  --pd1=0.00531  -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1732.log &

cd $TEMPLATEDIR/J1909-3744
$SCRIPTDIR/$SIMSCRIPT -d . -m J1909.simulate.mjd -p J1909.par --beam2fits 512  --pd1=0.002947 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1909.log &

cd $TEMPLATEDIR/J2129-5721
$SCRIPTDIR/$SIMSCRIPT -d . -m J2129.simulate.mjd -p J2129.par --beam2fits 512  --pd1=0.003726 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J2129.log &

#1024
cd $TEMPLATEDIR/J0437-4715
$SCRIPTDIR/$SIMSCRIPT -d . -m J0437.simulate.mjd -p J0437.par --beam2fits 1024 --pd1=0.005757 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J0437.log &

cd $TEMPLATEDIR/J0711-6830
$SCRIPTDIR/$SIMSCRIPT -d . -m J0711.simulate.mjd -p J0711.par --beam2fits 1024 --pd1=0.005491 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J0711.log &

cd $TEMPLATEDIR/J1024-0719
$SCRIPTDIR/$SIMSCRIPT -d . -m J1024.simulate.mjd -p J1024.par --beam2fits 1024 --pd1=0.005162 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1024.log &

cd $TEMPLATEDIR/J1603-7202
$SCRIPTDIR/$SIMSCRIPT -d . -m J1603.simulate.mjd -p J1603.par --beam2fits 1024 --pd1=0.014841 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1603.log &

cd $TEMPLATEDIR/J1713+0747
$SCRIPTDIR/$SIMSCRIPT -d . -m J1713.simulate.mjd -p J1713.par --beam2fits 1024 --pd1=0.004570 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1713.log &

cd $TEMPLATEDIR/J1730-2304
$SCRIPTDIR/$SIMSCRIPT -d . -m J1730.simulate.mjd -p J1730.par --beam2fits 1024 --pd1=0.008122 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1730.log &

cd $TEMPLATEDIR/J1744-1134
$SCRIPTDIR/$SIMSCRIPT -d . -m J1744.simulate.mjd -p J1744.par --beam2fits 1024 --pd1=0.004074 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1744.log &

cd $TEMPLATEDIR/J1857+0943
$SCRIPTDIR/$SIMSCRIPT -d . -m J1857.simulate.mjd -p J1857.par --beam2fits 1024 --pd1=0.005362 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1857.log &

cd $TEMPLATEDIR/J2124-3358
$SCRIPTDIR/$SIMSCRIPT -d . -m J2124.simulate.mjd -p J2124.par --beam2fits 1024 --pd1=0.004931 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J2124.log &

#2048
cd $TEMPLATEDIR/J1022+1001
$SCRIPTDIR/$SIMSCRIPT -d . -m J1022.simulate.mjd -p J1022.par --beam2fits 2048 --pd1=0.016452 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J1022.log &

cd $TEMPLATEDIR/J2145-0750
$SCRIPTDIR/$SIMSCRIPT -d . -m J2145.simulate.mjd -p J2145.par --beam2fits 2048 --pd1=0.016052 -M ti_cal,ii_cal -s $TEMPLATEDIR/$PARAMFILE > pp_J2145.log &

