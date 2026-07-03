# Iceberg

A simple code to read in a list of defect diameters {d}, and return stats about the possible log-normal distributions which are compatible with that list. Deduces information about the hidden iceberg given the visible portion above the waterline.

Daniel Mason
(c) UKAEA 2026



## Compile 

```
    mkdir build ; cd build 
    cmake ../.
    make 
    ctest
```

## Input file format

Should be a simple one column list, eg
```
 2.10
 3.45
 0.86
 ...
```

## Simple run

Generate some test data with 
```
    cd build
    ./test/test_genData -h
    ./test/test_genData -f ../data/test.dat -mu 4.0 -sigma 2.0 
```
This generates 1000 points with lognormal distribution mean 4.0, stdev 2.0.
The results are random, but I get
```
test_genData
^^^^^^^^^^^^
 generate data with known lognormal distribution.

 generation function
LogNormal[mu,sigma =       3.57770876,      0.47238073]
 68% confidence interval      2.2366066917668501      :   4.0000000000000000      :   5.7229552460508506
 <mean>   :    4.0000000000000000
 <stdev>  :    2.0000000000000000
 visibility function
LogisticFunction unset

 n real   :         1000
 n seen   :         1000
 <d>      :    3.9890698311590067
 <d²>     :    19.966208288164008
 stdev    :    2.0133380665703071

 pass
```


Then find the result. There is lots of data, the important lines being
```
./bin/iceberg -f ../data/test.dat 


 chi-square value               14.524352946136730
 alpha value                   0.97543075856914863
 observed defect count          1000.0000000000000
 expected defect count          1027.4454891048622
 observed pointdefect count     65990.865911736750
 expected point defect count    66905.262250763772
 68% confidence interval <d>    2.1560622244710137      :   3.9385416520069074      :   5.6768935306798474
```
The confidence interval is right - so we do indeed have the right answer. But this was easy...
Now try adding the invisibility function, centre 2.0, width 1.0, and record the number of point defects by giving the volume per atom. 
```
    ./test/test_genData -f ../data/test.dat -mu 4.0 -sigma 2.0 -d0 2.0 -w 1.0 -omega0 1.0

    
 n real   :         1310
 n seen   :         1000
 npd real :    74858.630381334049
 npd obs  :    70671.476368651332
 <d>      :    4.2474852922383945
 <d²>     :    21.749893780289700
 stdev    :    1.9258147555017382
```

Now some of the defects are missing - I tried to place 1176 but only 1000 made it to the output file. 
Try to fit with no visibility function - this is done by setting the width to a -ve number, "-w -1".

```
./bin/iceberg -f ../data/test.dat -w -1 

 chi-square value               29.305325766912286
 alpha value                   0.54567024393954910
 observed defect count          1000.0000000000000
 expected defect count          1000.0000000000000
 observed pointdefect count     70671.476368651391
 expected point defect count    70295.729215202547
 68% confidence interval <d>    2.5247010055816621      :   4.2526160782935856      :   5.9490616762664006

```
The answer is now not so good- it is overestimating the average size. Now try to fit with a visibility function,
```
./bin/iceberg -f ../data/test.dat  

 chi-square value               28.984273223314666
 alpha value                   0.43313910212346152
 observed defect count          1000.0000000000000
 expected defect count          1297.9596121135100
 observed pointdefect count     70671.476368651391
 expected point defect count    79750.101731338393
 68% confidence interval <d>    2.4003330183428973      :   4.0573710318732701      :   5.6838009270568639
 ```
This gets a much better guess for the distribution.
We can also force the correct number of point defects, if we know it. Here I know the correct number is 74859,
```
./bin/iceberg -f ../data/test.dat -vol 1e8 -rho 0.000749 -omega0 1

 chi-square value               29.935767339047434
 alpha value                   0.43641703226706613
 observed defect count          1000.0000000000000
 expected defect count          1179.2473136855754
 observed pointdefect count     70671.476368651391
 expected point defect count    74900.000000000015
 68% confidence interval <d>    2.4032303048198380      :   4.0872141811268046      :   5.7393234093803907
```






## Theory


### Visibility

We make the assumption that not all defects are detectable, because some are too small.
The form of the visibility function is not critically important, as we always need to report what is actually seen as well as a best guess for what we miss.
The visibility function sets the probability of recording a defects with diameter d is a logistic function:

```
    p(d) = 1/(1 + Exp[ -(d-d0)/w ])
```
which tends to zero for small defects (d-d0)/w << 0 and tends to one for large defects (d-d0)/w >> 0.

### Lognormal distribution

The lognormal distribution looks power-law-like under some circumstances, and gaussian-like under others, so is a good general fit. It also does not have any weight for d<0. The lognormal has a few possible definitions, we use
```
    P(d) = 1/(sqrt(2 π) σ d ) exp[ - (ln d / μ)²/(2 σ²) ]
```
which means that the scale μ has diameter(length) units and the shape σ is dimensionless.


