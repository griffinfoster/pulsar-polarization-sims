pulsar-polarization-sims
========================

Scripts for simulating the effects of polarization purity in pulsar timing in the paper *Intrinsic Instrumental Polarization and High-Precision Pulsar Timing* ([arXiv](http://arxiv.org/abs/1507.06839))

#### Required Software:

* [psrchive](http://psrchive.sourceforge.net/), follow http://psrchive.sourceforge.net/third/install.shtml

#### Useful stuff:

The profiles of the PTA MSPs are in their respective directories. To create a useable template, run:

```
pam -r 0.5 -e rr <filename>.SS
pdv -t <filename>.rr | sed '1 d' | awk '{print $4, $5, $6, $7}' > template.dat
```

plot fits:

```
psrplot -p flux <filename>.fits
/xs
```

simulated MJDs:

```
tempo2 -gr fake -ndobs 1 -nobsd 1 -randha y -start 50000 -end 50500 -rms auto -f <par file>
./parseTempo2fakeSim.py -f <simulate file>
```

beam2fits error: 256/512/1024/2048

```
vim PSRBeam/useful.inc
line 3: np=256 or 512 or 1024 or 2048
make beam2fits
```

generate par files:

```
psrcat -r <PULSAR NAME> > <par file>
```

tempo2 plot residuals:

```
tempo2 -gr plk -f J1603.par <tim file>
```

#### Profile Size:

* 256:
  * J1824-2452
  * J1939+2134
* 512:
  * J0613-0200
  * J1045-4509
  * J1600-3053
  * J1643-1224
  * J1732-5049
  * J1909-3744
  * J2129-5721
* 1024:
  * J0437-4715
  * J0711-6830
  * J1024-0719
  * J1603-7202
  * J1713+0747
  * J1730-2304
  * J1744-1134
  * J1857+0943
  * J2124-3358
* 2048:
  * J1022+1001
  * J2145-0750

