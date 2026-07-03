# Iceberg

A simple code to read in a list of defect diameters {d}, and return stats about the possible log-normal distributions which are compatible with that list. Deduces information about the hidden iceberg given the visible portion above the waterline.

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


